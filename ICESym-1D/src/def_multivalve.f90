module def_multivalve

contains
  
  subroutine solve_multivalve(Ucyl, Uref, type_val, Area_val, Area_pipe, &
       RHS, ga, Upipe, Uthroat, solved_case, flag)
    !
    !  Solves the non-linear equation system for the multi-valve model.
    !  The solver is based on the selection of cases, each of them
    !  corresponding to a particular set of equations. The cylinders'
    !  states are assumed constant for this solver.
    !
    !  solve_multivalve is called by: solve_cylinder
    !  solve_multivalve calls the following subroutines and 
    !  functions: solve_multivalve_eqcases, get_multivalve_case
    !
    implicit none

    integer, intent(in) :: type_val
    logical, intent(in) :: flag
    real*8, intent(in) :: ga,Area_pipe
    real*8, dimension(2), intent(in) :: Area_val
    real*8, dimension(3), intent(in) :: RHS
    real*8, dimension(3,2), intent(in) :: Ucyl,Uref
    real*8, dimension(3), intent(inout) :: Upipe
    real*8, dimension(3,2), intent(inout) :: Uthroat
    integer, intent(inout) :: solved_case

    integer :: old_case, try_case, caso, max_iter, i, icase, iter, &
         n_P
    logical :: conv_flag, case_flag, do_cases, caso_ok
    integer, dimension(24) :: CASOS, CASOS_dbg
    real*8 :: rho_P, v_P, p_P, a_P, Up
    real*8, dimension(2) :: psi, rho_C, p_C, a_C
    real*8, dimension(3) :: U_P, UP0
    real*8, dimension(3,2) :: U_T

    n_P   = type_val

    U_P = Upipe
    UP0 = Upipe
    U_T = Uthroat
    psi(1) = Area_val(1)/Area_pipe
    psi(2) = Area_val(2)/Area_pipe

    do i=1,2
       p_C(i)   = Ucyl(2,i)
    end do

    call multivalve_cases(solved_case, CASOS, p_C)
    CASOS_dbg = CASOS
    if(solved_case.eq.0) solved_case = CASOS(1)
    old_case = solved_case
    try_case = old_case
    solved_case = 0
    max_iter = 24

    if(psi(1).gt.1d-6 .and. psi(2).gt.1d-6) then
       do_cases = .true.
       do i=1,2
          rho_C(i) = Ucyl(1,i)
          a_C(i)   = dsqrt(ga*p_C(i)/rho_C(i))
       end do
       conv_flag = .false.
       iter = 0
       do while (do_cases)
          iter = iter+1
          rho_P = Upipe(1)
          v_P   = Upipe(2)
          p_P   = Upipe(3)
          a_P   = dsqrt(ga*p_P/rho_P)

          if(v_P == 0) v_P = n_P
          Up = v_P*n_P/a_C(1)
          call solve_multivalve_eqcases(Up, Uref(:,2), Uref(:,1), RHS, Ucyl, &
               n_P, psi, ga, try_case, UP0, Upipe, Uthroat, case_flag)
          if((try_case.eq.3 .or. try_case.eq.9 .or. try_case.eq.15) &
               .and.case_flag) then
             do i=1,1
                call solve_multivalve_eqcases(Up, Uref(:,2), Uref(:,1), RHS, &
                     Ucyl, n_P, psi, ga, try_case, Upipe, Upipe, Uthroat, &
                     case_flag)
             end do
          end if

          do i=1, 24
             if(CASOS(i)==try_case) then
                CASOS(i)=0
                exit
             end if
          end do
          if(case_flag) then
             ! Here, we must to check the case
             ! call get_multivalve_case(Upipe, Uthroat, Ucyl, n_P, ga, caso)
             call check_multivalve_case(Upipe, Uthroat, Ucyl, n_P, ga, &
                  try_case, caso_ok)
             if(caso_ok) then
                solved_case = try_case
                conv_flag = .true.
             elseif(caso.ne.0) then
                icase = 0
                do i=1, 24
                   if(CASOS(i)==caso) then
                      icase=i
                      exit
                   end if
                end do
                if(icase.ne.0) then
                   try_case = CASOS(icase)
                else
                   icase = 0
                   do i=1, 24
                      if(CASOS(i).ne.0) then
                         icase=i
                         exit
                      end if
                   end do
                   if(icase.ne.0) then
                      try_case = CASOS(icase)
                   else
                      try_case = 0
                   end if
                end if
             else
                icase = 0
                do i=1, 24
                   if(CASOS(i).ne.0) then
                      icase=i
                      exit
                   end if
                end do
                if(icase.ne.0) then
                   try_case = CASOS(icase)
                else
                   try_case = 0
                end if
             end if
          else
             icase = 0
             do i=1, 24
                if(CASOS(i).ne.0) then
                   icase=i
                   exit
                end if
             end do
             if(icase.ne.0) then
                try_case = CASOS(icase)
             else
                try_case = 0
             end if
          end if

          if(conv_flag .or. try_case==0) &
               do_cases = .false.
       end do
    end if

  end subroutine solve_multivalve

  subroutine solve_multivalve_eqcases(Up, U_S, U_R, RHS, U_C, n_P, psi, ga, &
       caso, UP_old, U_P, U_T, ok_flag)
    !
    !  Solves the non-linear equation arising from the multivalve model for
    !  the specified case
    !
    !  solve_multivalve_eqcases is called by:
    !  solve_multivalve_eqcases calls the following subroutines and 
    !  functions:
    !
    implicit none

    integer, intent(in) :: caso, n_P
    real*8, intent(in) :: ga
    real*8, dimension(2), intent(in) :: psi
    real*8, dimension(3), intent(in) :: U_S, U_R, RHS, UP_old
    real*8, dimension(3,2), intent(in) :: U_C
    real*8, intent(inout) :: Up
    real*8, dimension(3), intent(out) :: U_P
    real*8, dimension(3,2), intent(out) :: U_T
    logical, intent(out) :: ok_flag

    integer :: max_ils,max_iter,iter,ils,caso_here
    real*8 :: rho_S,v_S,p_S,a_S,rho_R,v_R,p_R,a_R, &
         rho_P,v_P,p_P,a_P,rho_T1,v_T1,p_T1,a_T1, &
         rho_T2,v_T2,p_T2,a_T2
    real*8 :: Up_new,RHS1,RHS2,delta,F,dF,tol,dU,nu,relax,g1, &
         pC_tmp,rhoC_tmp,aC_tmp
    real*8, dimension(2) :: p_C,a_C,rho_C,psi_here
    real*8, dimension(3) :: UT_tmp
    logical :: do_newton,do_ls,swapped,F_isreal
    real*8 :: Fnew

    caso_here = caso
    g1 = ga-1.
    delta = 0.5*g1

    tol = 1e-9
    max_iter = 20
    max_ils  = 10
    nu = 2.0

    F = huge(0.0d0)

    ok_flag   = .false.
    do_newton = .false.

    swapped   = .false.
    if(caso_here==2.or.caso_here==5.or. &
         (caso_here>=13.and.caso_here<=18).or. &
         caso_here==20.or.caso_here==23) swapped = .true.

    rho_C = U_C(1,:)
    p_C   = U_C(2,:)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_R = U_R(1)
    v_R   = U_R(2)
    p_R   = U_R(3)

    psi_here = psi

    a_S = dsqrt(ga*p_S/rho_S)
    a_R = dsqrt(ga*p_R/rho_R)
    a_P = dsqrt(ga*UP_old(3)/UP_old(1))
    a_C(1) = dsqrt(ga*p_C(1)/rho_C(1))
    a_C(2) = dsqrt(ga*p_C(2)/rho_C(2))

    RHS1 = RHS(1)
    if(n_P < 0) then
       RHS2 = RHS(2)
    elseif(n_P > 0) then
       RHS2 = RHS(3)
    end if

    if(swapped) then
       Up = Up*a_C(1)/a_C(2)

       pC_tmp    = p_C(1)
       aC_tmp    = a_C(1)
       rhoC_tmp  = rho_C(1)
       p_C(1) = p_C(2)
       p_C(2) = pC_tmp
       rho_C(1) = rho_C(2)
       rho_C(2) = rhoC_tmp
       a_C(1) = a_C(2)
       a_C(2) = aC_tmp
       psi_here(1) = psi(2)
       psi_here(2) = psi(1)

       select case (caso_here)
       case (2)
          caso_here = 1
       case (5)
          caso_here = 4
       case (13)
          caso_here = 8
       case (14)
          caso_here = 7
       case (15)
          caso_here = 9
       case (16)
          caso_here = 11
       case (17)
          caso_here = 10
       case (18)
          caso_here = 12
       case (20)
          caso_here = 19
       case (23)
          caso_here = 22
       end select
    end if

    call check_initial_guess_multivalve(Up, v_S, p_S, a_S, p_R, a_R, &
         p_C, a_C, a_P, RHS1, RHS2, caso_here, n_P, ga)

    call compute_multivalve_res(Up, v_S, p_S, a_S, p_R, a_R, p_C, &
         a_C, a_P, RHS1, RHS2, psi_here, caso_here, n_P, ga, Fnew, F_isreal)

    if(F_isreal) then
       if(caso_here==4 .or. caso_here==6 .or. caso_here==10 .or. &
            caso_here==11 .or. caso_here==12 .or. &
            caso_here==22 .or. caso_here==24) then
          if(Up > 0.) do_newton=.true.
       elseif(caso_here==1 .or. caso_here==3 .or. caso_here==7 .or. &
            caso_here==8 .or. caso_here==9 .or. &
            caso_here==19 .or. caso_here==21) then
          if(Up < 0.) do_newton=.true.
       end if
       F = Fnew
    end if

    ! Newton loop
    iter = 0
    do while(do_newton)
       do_ls = .true.
       iter  = iter+1
       relax = 1.0
       call compute_multivalve_jaco(Up, v_S, p_S, a_S, p_R, a_R, &
            p_C, a_C, a_P, RHS1, RHS2, psi_here, caso_here, n_P, ga, dF)
       dU = -F/dF
       ! line-search loop
       ils = 0
       do while(do_ls)
          ils = ils+1
          Up_new = Up+relax*dU
          call compute_multivalve_res(Up_new, v_S, p_S, a_S, p_R, a_R, &
               p_C, a_C, a_P, RHS1, RHS2, psi_here, caso_here, n_P, ga, &
               Fnew, F_isreal)
          if(F_isreal .and. dabs(Fnew).lt.dabs(F)) then
             F  = Fnew
             Up = Up_new
             exit
          else
             relax = relax/nu
          end if
          if(ils.ge.max_ils) exit
       end do ! ends line-search loop
       if(dabs(F).lt.tol) then
          ok_flag =.true.
          exit
       end if
       if(ils.ge.max_ils .or. iter.gt.max_iter) exit
    end do ! ends newton loop

    !--------------------------------------------------------------------------
    if(ok_flag) then
       if(caso_here==1.or.caso_here==4.or.caso_here==7.or. &
            caso_here==10.or.caso_here==19.or.caso_here==22) then
          v_P = a_C(2)*Up*n_P
          p_P = -p_S*ga*a_C(2)/a_S*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       elseif(caso_here==3.or.caso_here==6.or.caso_here==8.or. &
            caso_here==9.or.caso_here==11.or.caso_here==12.or. &
            caso_here==21.or.caso_here==24) then
          v_P = a_C(1)*Up*n_P
          p_P = -p_S*ga*a_C(1)/a_S*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       end if
       if(caso_here==1.or.caso_here==4) then
          if(caso_here==1) then
             rho_P = ga*p_P/(a_C(2)**2*(1.-delta*Up**2))
          elseif(caso_here==4) then
             rho_P = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
             ! rho_P = ga*p_P/(a_C(2)**2*(1.-delta*Up**2))
          endif
          a_P    = dsqrt(ga*p_P/rho_P)
          p_T2   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -a_C(2)*dsqrt((1.-ga*p_T2/rho_T2/a_C(2)**2.)/delta)*n_P
          p_T1   = p_C(1)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1   = (rho_P*v_P-rho_T2*psi_here(2)*v_T2)/rho_T1/psi_here(1)
       elseif(caso_here==3) then
          p_T2   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -a_C(2)*dsqrt((1.-ga*p_T2/rho_T2/a_C(2)**2)/delta)*n_P
          p_T1   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T1 = rho_C(1)*(p_T1/p_C(1))**(1./ga)
          v_T1   = -a_C(1)*dsqrt((1.-ga*p_T1/rho_T1/a_C(1)**2)/delta)*n_P
          if(dabs(v_P) > 0.) then
             rho_P = 1./v_P*(rho_T1*psi_here(1)*v_T1+rho_T2*psi_here(2)*v_T2)
          else
             write(*,*) ' Warning! - Multivalve Case 3 '
          end if
       elseif(caso_here==9) then
          p_T2   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -a_C(2)*dsqrt((1.-ga*p_T2/rho_T2/a_C(2)**2)/delta)*n_P
          p_T1   = p_C(1)/(delta+1.)**(ga/g1)
          rho_T1 = rho_C(1)*(p_T1/p_C(1))**(1./ga)
          v_T1   = -a_C(1)/dsqrt(delta+1.)*n_P
          if(dabs(v_P) > 0.) then
             rho_P = 1./v_P*(rho_T1*psi_here(1)*v_T1+rho_T2*psi_here(2)*v_T2)
          else
             write(*,*) ' Warning! - Multivalve Case 9 '
          end if
       elseif(caso_here==6) then
          rho_P  = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          a_P    = dsqrt(ga*p_P/rho_P)
          p_T2   = p_C(2)
          rho_T2 = rho_P*(p_T2/p_P)**(1./ga)
          v_T2   = a_P*dsqrt((1.-ga*p_T2/rho_T2/a_P**2)/delta+(v_P/a_P)**2)*n_P
          p_T1   = p_C(1)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1 = 1./(rho_T1*psi_here(1))*(rho_P*v_P-rho_T2*psi_here(2)*v_T2)
       elseif(caso_here==8) then
          rho_P  = ga*p_P/(a_C(1)**2*(1.-delta*Up**2))
          p_T2   = p_C(2)
          rho_T2 = rho_P*(p_T2/p_P)**(1./ga)
          v_T2   = a_C(1)*dsqrt((1.-ga*p_T2/rho_T2/a_C(1)**2)/delta)*n_P
          p_T1   = p_C(1)/(delta+1.)**(ga/g1)
          rho_T1 = rho_C(1)*(p_T1/p_C(1))**(1./ga)
          v_T1   = -dsqrt(ga*p_T1/rho_T1)*n_P
       elseif(caso_here==7) then
          rho_P  = ga*p_P/(a_C(2)**2*(1.-delta*Up**2))
          a_P    = dsqrt(ga*p_P/rho_P)
          p_T2   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -a_C(2)*dsqrt((1.-ga*p_T2/rho_T2/a_C(2)**2)/delta)*n_P
          p_T1   = p_P*(a_C(2)**2*rho_P/(ga*(delta+1.)*p_P))**(ga/g1)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1   = dsqrt(ga*p_T1/rho_T1)*n_P
       elseif(caso_here==10) then
          rho_P  = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          a_P    = dsqrt(ga*p_P/rho_P)
          p_T2   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          a_T2   = dsqrt(ga*p_T2/rho_T2)
          v_T2   = -a_C(2)*dsqrt((1.-(a_T2/a_C(2))**2.)/delta)*n_P
          a_T1   = a_P*dsqrt((1.+delta*(v_P/a_P)**2.)/(1.+delta))
          p_T1   = p_P*(a_T1/a_P)**(2.*ga/g1)
          rho_T1 = ga*p_T1/a_T1**2.
          v_T1   = a_T1*n_P
       elseif(caso_here==11) then
          rho_P  = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          p_T1   = p_C(1)/(delta+1.)**(ga/g1)
          rho_T1 = rho_C(1)*(p_T1/p_C(1))**(1./ga)
          v_T1   = -a_C(1)/dsqrt(delta+1.)*n_P
          p_T2   = p_C(2)
          rho_T2 = rho_C(1)*(a_C(1)/a_R)**2.*dexp(-RHS1)* &
               (p_R/p_C(1))**(g1/ga)*(p_C(2)/p_C(1))**(1./ga)
          v_T2   = 1./(rho_T2*psi_here(2))*(rho_P*v_P-rho_T1*psi_here(1)*v_T1)
       elseif(caso_here==12) then
          rho_P  = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          a_P    = dsqrt(ga*p_P/rho_P)
          a_T1   = a_P*dsqrt(1./(delta+1.)*(1.+delta*(v_P/a_P)**2))
          a_T2   = a_R*(p_C(2)/p_R)**(g1/2./ga)*dexp(0.5*RHS1)
          p_T1   = p_P*(a_T1/a_P)**(2.*ga/g1)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1   = a_T1*n_P
          p_T2   = p_C(2)
          rho_T2 = rho_P*(p_T2/p_P)**(1./ga)
          v_T2   = dsqrt((a_T1**2*(delta+1.)-a_T2**2)/delta)*n_P
       elseif(caso_here==19) then
          rho_P  = ga*p_P/(a_C(2)**2*(1.-delta*Up**2))
          p_T2   = p_C(2)/(delta+1.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -dsqrt(ga*p_T2/rho_T2)*n_P
          p_T1   = p_P/((delta+1.)*(1.-delta*Up**2))**(ga/g1)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1   = dsqrt(ga*p_T1/rho_T1)*n_P
       elseif(caso_here==21) then
          p_T2   = p_C(2)/(delta+1.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -dsqrt(ga*p_T2/rho_T2)*n_P
          p_T1   = p_C(1)/(delta+1.)**(ga/g1)
          rho_T1 = rho_C(1)*(p_T1/p_C(1))**(1./ga)
          v_T1   = -dsqrt(ga*p_T1/rho_T1)*n_P
          if(dabs(v_P) > 0.) then
             rho_P = 1./v_P*(rho_T1*psi_here(1)*v_T1+rho_T2*psi_here(2)*v_T2)
          else
             write(*,*) ' Warning! - Multivalve Case 21 '
          end if
       elseif(caso_here==22) then
          rho_P = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          p_T2   = p_C(2)/(delta+1.)**(ga/g1)
          rho_T2 = rho_C(2)*(p_T2/p_C(2))**(1./ga)
          v_T2   = -dsqrt(ga*p_T2/rho_T2)*n_P
          rho_T1 = (rho_P**ga/ga/p_P* &
               (rho_P*v_P-rho_T2*psi(2)*v_T2)**2/psi(1)**2)**(1./(ga+1.))
          p_T1   = p_P*(rho_T1/rho_P)**ga
          v_T1   = dsqrt(ga*p_T1/rho_T1)*n_P
       elseif(caso_here==24) then
          rho_P  = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
          a_P    = dsqrt(ga*p_P/rho_P)
          a_T1   = a_P*dsqrt((1.+delta*(v_P/a_P)**2)/(1.+delta))
          p_T1   = p_P*(a_T1/a_P)**(ga/delta)
          rho_T1 = rho_P*(p_T1/p_P)**(1./ga)
          v_T1   = a_T1*n_P
          p_T2   = p_T1
          rho_T2 = rho_T1
          v_T2   = v_T1
       end if

       U_P(1) = rho_P
       U_P(2) = v_P
       U_P(3) = p_P

       U_T(1,1) = rho_T1
       U_T(2,1) = v_T1
       U_T(3,1) = p_T1

       U_T(1,2) = rho_T2
       U_T(2,2) = v_T2
       U_T(3,2) = p_T2
    end if

    if(swapped) then
       UT_tmp = U_T(:,1)
       U_T(:,1) = U_T(:,2)
       U_T(:,2) = UT_tmp
    end if

  end subroutine solve_multivalve_eqcases

  subroutine check_initial_guess_multivalve(Up, v_S, p_S, a_S ,p_R, a_R, p_C, a_C, a_P, &
       RHS1, RHS2, caso, n_P, ga)
    !
    !
    !  check_initial_guess_multivalve is called by:
    !  check_initial_guess_multivalve calls the following subroutines 
    !  and functions:
    !
    implicit none

    integer, intent(in) :: caso,n_P
    real*8, intent(in) :: ga,a_S,v_S,p_S,a_R,p_R,a_P, &
         RHS2,RHS1
    real*8, dimension(2), intent(in) :: p_C,a_C
    real*8, intent(inout) :: Up

    integer :: iter,max_iter
    real*8 :: C0,C1,C2,CR,K1,K2,a_T2,F,dF
    real*8 :: delta,g1,Up0,Up1,Up2,tol_ini,tol

    g1 = ga-1.
    delta = 0.5*g1

    tol_ini = 1e-8
    tol = 1e-10
    max_iter = 20

    select case (caso)
    case (1)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       C1 = -a_S/a_C(2)*p_C(1)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K1 = 1.-C1**(g1/ga)*(1.-delta*Up**2.)
       if(K1 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C1 = -ga*a_C(2)/a_S*p_S/p_C(1)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
             F  = 1.-delta*Up**2.-C1**(g1/ga)
             dF = -g1*Up+g1*a_C(2)/a_S*p_S/p_C(1)/C1**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = C2**(g1/ga)*(1.-delta*Up**2.)-1.
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
             F  = 1.-delta*Up**2.-C2**(g1/ga)
             dF = -g1*Up+g1*a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up + tol_ini
       end if
    case (3)
       C1 = -ga*a_C(1)/a_S*p_S/p_C(1)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K1 = 1.-C1**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)
       if(K1 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
             F  = 1.+delta*Up**2.*(a_C(1)/a_P)**2.-C1**(g1/ga)
             dF = g1*Up*(a_C(1)/a_P)**2.-g1*a_C(1)/a_S*p_S/p_C(1)*C1**((2.*ga-1.)/ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up*(1. - tol_ini)
       end if
       C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2 = 1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
             F  = 1.+delta*Up**2.*(a_C(1)/a_P)**2.-C2**(g1/ga)
             dF = g1*Up*(a_C(1)/a_P)**2.-g1*a_C(1)/a_S*p_S/p_C(2)*C2**((2.*ga-1.)/ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up*(1. - tol_ini)
       end if
    case (4)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       C0 = (a_R/a_C(2))**2.*dexp(RHS1)
       CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K1 = C0*CR**(g1/ga)+delta*Up**2.-C0*(p_C(1)/p_R)**(g1/ga)
       if(K1 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+ &
                  RHS2/a_C(2))
             F  = CR**(g1/ga)+delta*Up**2./C0-(p_C(1)/p_R)**(g1/ga)
             dF = g1*Up/C0-g1*a_C(2)/a_S*p_S/p_R/CR**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
       C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = 1.-C2**(g1/ga)-delta*Up**2./C0*(p_R/p_C(2))**(g1/ga)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+ &
                  RHS2/a_C(2))
             F  = C2**(g1/ga)+delta*Up**2./C0*(p_R/p_C(2))**(g1/ga)-1.
             dF = g1*Up/C0*(p_R/p_C(2))**(g1/ga)- &
                  g1*a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up + tol_ini
       end if
    case (6)
       C0  = (a_R/a_C(1))**2.*dexp(RHS1)
       C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       CR = C1*p_R/p_C(1)
       K1 = 1.+delta*Up**2./C0*CR**(g1/ga)-C1**(g1/ga)
       if(K1 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter  = iter+1
             C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+ &
                  RHS2/a_C(1))
             CR = C1*p_R/p_C(1)
             F  = 1.+delta*Up**2./C0*CR**(g1/ga)-C1**(g1/ga)
             dF = g1*Up/C0*CR**(g1/ga)+ &
                  g1*delta*Up**2./C0*a_C(1)/a_S*p_S/p_R*CR**((2.*ga-1)/ga)- &
                  g1*a_C(1)/a_S*p_S/p_C(1)*C1**((2.*ga-1.)/ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
       C2 = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       CR = C2*p_R/p_C(2)
       K2 = 1.+delta*Up**2./C0*CR**(g1/ga)-C2**(g1/ga)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter  = iter+1
             C2 = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+ &
                  RHS2/a_C(1))
             CR = C2*p_R/p_C(2)
             F  = 1.+delta*Up**2./C0*CR**(g1/ga)-C2**(g1/ga)
             dF = g1*Up/C0*CR**(g1/ga)+ &
                  g1*delta*Up**2./C0*a_C(1)/a_S*p_S/p_R*CR**((2.*ga-1)/ga)- &
                  g1*a_C(1)/a_S*p_S/p_C(2)*C2**((2.*ga-1.)/ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
    case (7)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = C2**(g1/ga)*(1.-delta*Up**2.)-1.
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
             F  = 1.-delta*Up**2.-C2**(g1/ga)
             dF = -g1*Up+g1*a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up*(1. - tol_ini)
       end if
    case (8)
       C2  = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2  = 1.-C2**(g1/ga)*(1.-delta*Up**2)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+ &
                  RHS2/a_C(1))
             F  = 1.-delta*Up**2.-C2**(g1/ga)
             dF = -g1*Up+g1*a_C(1)/a_S*p_S/p_C(2)/C2**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
    case (9)
       ! tol_ini = 1.0d-4
       if(Up.gt.0.) Up = -tol_ini
       C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2 = 1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+ &
                  RHS2/a_C(1))
             F  = 1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)
             dF = g1*a_C(1)/a_S*p_S/p_C(2)/C2**(1./ga)* &
                  (1.+delta*Up**2.*(a_C(1)/a_P)**2.)- &
                  g1*Up*(a_C(1)/a_P)**2*C2**(g1/ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up + tol_ini
       end if
    case (10)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       C0 = (a_R/a_C(2))**2.*dexp(RHS1)
       C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = 1.-C2**(g1/ga)-delta*Up**2./C0*(p_R/p_C(2))**(g1/ga)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             C2 = -ga*a_C(2)/a_S*p_S/p_C(2)*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+ &
                  RHS2/a_C(2))
             F  = C2**(g1/ga)+delta*Up**2./C0*(p_R/p_C(2))**(g1/ga)-1.
             dF = g1*Up/C0*(p_R/p_C(2))**(g1/ga)- &
                  g1*a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
    case (11)
       C0 = (a_R/a_C(1))**2.*dexp(RHS1)
       CR = -ga*a_C(1)/a_S*p_S/p_R*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2 = CR**(g1/ga)+delta*Up**2./C0-(p_C(2)/p_R)**(g1/ga)
       if(K2 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          ! the limits of the interval where the solution can be found
          iter = 0
          do while (.true.)
             iter = iter+1
             CR = -ga*a_C(1)/a_S*p_S/p_R*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+ &
                  RHS2/a_C(1))
             F  = CR**(g1/ga)+delta*Up**2./C0-(p_C(2)/p_R)**(g1/ga)
             dF = g1*Up/C0-g1*a_C(1)/a_S*p_S/p_R/CR**(1./ga)
             Up = Up-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up - tol_ini
       end if
    case (12)
       Up1 = a_S/a_C(1)/ga*(1.-p_C(1)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       Up2 = a_S/a_C(1)/ga*(1.-p_C(2)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       if(Up > Up1 .or. Up > Up2) Up = dmin1(Up1, Up2) - tol_ini
       C0   = a_R/a_C(1)*dexp(0.5*RHS1)
       a_T2 = a_R*(p_C(2)/p_R)**(g1/2./ga)*dexp(0.5*RHS1)
       Up0  = a_S/a_C(1)/ga+v_S/a_C(1)*n_P-RHS2/a_C(1)
       if(Up < 0.) Up = tol_ini
       CR = -a_S/a_C(1)*p_R/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = CR*p_C(2)/p_R
       K1 = C0**2/CR**(g1/ga)+delta*Up**2-(a_T2/a_C(1))**2
       if(K1 <= 0.) then
          ! Here, we solve a non-linear equation in order to get
          !the limits of the interval where the solution can be found
          Up1  = tol_ini
          iter = 0
          do while (.true.)
             iter  = iter+1
             F  = C0**2/(a_S/a_C(1)*p_R/p_S/ga)**(g1/ga)*(Up0-Up1)**(g1/ga)+ &
                  delta*Up1**2-(a_T2/a_C(1))**2
             dF = -g1/ga*C0**2/(a_S/a_C(1)*p_R/p_S/ga)**(g1/ga)*(Up0-Up1)**(-1./ga)+ &
                  2.*delta*Up1
             Up1 = Up1-F/dF
             if(dabs(F) < tol .or. iter > max_iter) exit
          end do
          Up = Up1 - tol_ini
       end if
    case (19)
       Up1 = a_S/a_C(2)/ga*(1.-p_C(1)/p_S)+v_S/a_C(2)*n_P-RHS2/a_C(2)
       Up2 = a_S/a_C(2)/ga*(1.-p_C(2)/p_S)+v_S/a_C(2)*n_P-RHS2/a_C(2)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       if(Up < Up2 .or. Up > Up1 .or. Up > 0.) &
            Up = dmax1(Up2+tol_ini, Up1-tol_ini)
    case (21)
       Up1 = a_S/a_C(1)/ga*(1.-p_C(1)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       Up2 = a_S/a_C(1)/ga*(1.-p_C(2)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       if(Up < dmax1(Up1,Up2) .or. Up > 0.) Up = dmax1(Up1,Up2) + tol_ini
    case (22)
       Up1 = a_S/a_C(2)/ga*(1.-p_C(1)/p_S)+v_S/a_C(2)*n_P-RHS2/a_C(2)
       Up2 = a_S/a_C(2)/ga*(1.-p_C(2)/p_S)+v_S/a_C(2)*n_P-RHS2/a_C(2)
       Up0 = a_S/a_C(2)/ga+v_S/a_C(2)*n_P-RHS2/a_C(2)
       ! v_P is allways non-dimensionalized by a_C(1),
       ! so we must redefine Up in this case
       Up = Up*a_C(1)/a_C(2)
       if(Up < 0. .or. Up > Up0) Up = Up0 - tol_ini
       if(Up < Up2 .or. Up > Up1) Up = dmin1(Up2+tol_ini, Up1-tol_ini)
    case (24)
       Up1 = a_S/a_C(1)/ga*(1.-p_C(1)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       Up2 = a_S/a_C(1)/ga*(1.-p_C(2)/p_S)+v_S/a_C(1)*n_P-RHS2/a_C(1)
       Up0 = a_S/a_C(1)/ga+v_S/a_C(1)*n_P-RHS2/a_C(1)
       if(Up < 0. .or. Up > Up0) Up = Up0 - tol_ini
       if(Up > Up1 .or. Up > Up2) Up = dmin1(Up2, Up1)-tol_ini
    end select

  end subroutine check_initial_guess_multivalve

  subroutine compute_multivalve_res(Up, v_S, p_S, a_S, p_R, a_R, p_C, a_C, a_P, &
       RHS1, RHS2, psi, caso, n_P, ga, F, isreal)
    !
    !
    !  compute_multivalve_res is called by:
    !  compute_multivalve_res calls the following subroutines 
    !  and functions: none
    !
    implicit none

    integer, intent(in) :: caso,n_P
    real*8, intent(in) :: Up,ga,a_S,v_S,p_S,a_R,p_R,a_P, &
         RHS2,RHS1
    real*8, dimension(2), intent(in) :: p_C,a_C,psi
    real*8, intent(out) :: F
    logical, intent(out) :: isreal

    real*8 :: C0,C1,C2,CR,K0,rho_T2,a_T2,rho_C1
    real*8 :: delta,g1
    real*8 :: K2,K1

    g1 = ga-1.
    delta = 0.5*g1

    F = huge(0.0d0)
    isreal = .false.

    select case (caso)
    case (1)
       C1 = -a_S/a_C(2)/ga*p_C(1)/p_S/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = C1*p_C(2)/p_C(1)
       K1 = (1.-(1.-delta*Up**2)*C1**(g1/ga))/delta
       K2 = (C2**(g1/ga)*(1.-delta*Up**2)-1.)/delta
       if(K1.gt.0. .and. K2.gt.0.) then
          F = Up-psi(1)*C1**(1./ga)*dsqrt(K1)+ &
               psi(2)*C2**(delta/ga)*(1.-delta*Up**2)**((ga-3.)/2./g1)*dsqrt(K2)
          isreal = .true.
       end if
    case (3)
       C1 = -ga*a_C(1)/a_S*p_S/p_C(1)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(1)/p_C(2)
       K1 = (1.-C1**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta
       K2 = (1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta
       if(K1.gt.0. .and. K2.gt.0.) then
          F  = Up*C1**(g1/ga)/(1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)+ &
               psi(1)*(1.-delta*Up**2.)*dsqrt(K1)+ &
               psi(2)*(p_C(2)/p_C(1))**(g1/ga)*a_C(1)/a_C(2)* &
               ((a_C(2)/a_C(1))**2.-delta*Up**2.)*dsqrt(K2)
          isreal = .true.
       end if
    case (4)
       C0 = (a_R/a_C(2))**2.*dexp(RHS1)
       C1 = -a_S/a_C(2)/ga*p_C(1)/p_S/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = p_C(1)/p_C(2)/C1
       CR = C1*p_R/p_C(1)
       K1 = (C0/CR**(g1/ga)+delta*Up**2.-C0*(p_C(1)/p_R)**(g1/ga))/delta
       K2 = (1.-C2**(g1/ga)-delta*Up**2./C0*(p_R/p_C(2))**(g1/ga))/delta
       if(K1.gt.0. .and. K2.gt.0.) then
          F = Up-psi(1)*C1**(1./ga)*dsqrt(K1)+ &
               psi(2)*C0*(p_C(2)/p_R)**(g1/ga)* &
               (1.+delta*Up**2./C0*CR**(g1/ga))**(1./g1)*dsqrt(K2)
          isreal = .true.
       end if
    case (6)
       C0 = a_C(1)/a_R*dexp(-0.5*RHS1)
       C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(2)/p_C(1)
       CR = C1*p_R/p_C(1)
       K1 = (1.+delta*Up**2*C0**2*CR**(g1/ga)-C1**(g1/ga))/delta
       K2 = (1.+delta*Up**2*C0**2*CR**(g1/ga)-C2**(g1/ga))/delta
       if(K1.gt.0. .and. K2.gt.0.) then
          F  = C0*Up*CR**(delta/ga)-psi(1)*C1**(1./ga)*dsqrt(K1)- &
               psi(2)*C2**(1./ga)*dsqrt(K2)
          isreal = .true.
       end if
    case (7)
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = (C2**(g1/ga)*(1.-delta*Up**2.)-1.)/delta
       if(K2.gt.0.) then
          F  = Up*(1.-delta*Up**2)**(1./g1)-psi(1)/(delta+1.)**((ga+1.)/2./g1)+ &
               psi(2)*C2**(delta/ga)*dsqrt(K2*(1.-delta*Up**2))
          isreal = .true.
       end if
    case (8)
       C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(2)/p_C(1)
       K2 = (1.-C2**(g1/ga)*(1.-delta*Up**2))/delta
       if(K2.gt.0.) then
          F  = Up+psi(1)/(delta+1.)**((ga+1.)/g1/2.)*(1.-delta*Up**2)*C1- &
               psi(2)*C2**(1./ga)*dsqrt(K2)
          isreal = .true.
       end if
    case (9)
       C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2 = (1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta
       if(K2.gt.0.) then
          F  = Up*C2+psi(1)/(delta+1.)**((ga+1.)/2./g1)*p_C(1)/p_C(2)*(1.-delta*Up**2.)+ &
               psi(2)*a_C(2)/a_C(1)*C2**(1./ga)* &
               (1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)* &
               (1.-delta*Up**2*(a_C(1)/a_C(2))**2.)*dsqrt(K2)
          isreal = .true.
       end if
    case (10)
       C0 = a_R/a_C(2)*dexp(0.5*RHS1)
       CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = CR*p_R/p_C(2)
       K2 = (1.-C2**(g1/ga)-delta*Up**2./C0**2.*(p_R/p_C(2))**(g1/ga))/delta
       K1 = 1.+delta*Up**2./C0**2./CR**(g1/ga)
       if(K2.gt.0.) then
          F  = Up-psi(1)/(1.+delta)**((ga+1.)/2./g1)*C0* &
               CR**(delta/ga)*K1**((ga+1.)/2./g1)+ &
               psi(2)*C0**2.*(p_C(2)/p_R)**(g1/ga)*K1**(1./g1)*dsqrt(K2)
          isreal = .true.
       end if
    case (11)
       C0 = (a_R/a_C(1))**2.*dexp(RHS1)
       C2 = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       CR = C2*p_R/p_C(2)
       K2 = (C0/CR**(g1/ga)+delta*Up**2.-C0*(p_C(2)/p_R)**(g1/ga))/delta
       if(K2.gt.0.) then
          F  = Up+psi(1)/(1.+delta)**((ga+1.)/2./g1)*C0*p_C(1)/p_R*CR**(1./ga)- &
               psi(2)*C2**(1./ga)*dsqrt(K2)
          isreal = .true.
       end if
    case (12)
       C0 = a_R/a_C(1)*dexp(0.5*RHS1)
       a_T2 = a_R*(p_C(2)/p_R)**(g1/2./ga)*dexp(0.5*RHS1)
       CR = -a_S/a_C(1)*p_R/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = CR*p_C(2)/p_R
       K0 = 1./(delta+1.)*(1.+delta/C0**2*CR**(g1/ga)*Up**2)
       K1 = (C0**2/CR**(g1/ga)+delta*Up**2-(a_T2/a_C(1))**2)/delta
       if(K1.gt.0.) then
          F  = Up-psi(1)*C0/CR**(g1/2./ga)*K0**((ga+1.)/2./g1)- &
               psi(2)*C2**(1./ga)*dsqrt(K1)
          isreal = .true.
       end if
    case (19)
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       F = Up*(1.-delta*Up**2)**(1./g1)-psi(1)/(delta+1.)**(0.5/g1)+ &
            C2*(1.-delta*Up**2)**(ga/g1)*psi(2)/(delta+1.)**((ga+1.)/2./g1)
       isreal = .true.
    case (21)
       C1 = -ga*a_C(1)/a_S*p_S/p_C(1)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       F = (delta+1.)**((ga+1.)/2./g1)*Up*C1+psi(1)- &
            delta*Up**2*(psi(1)+psi(2)*p_C(2)/p_C(1)*a_C(1)/a_C(2))+ &
            psi(2)*p_C(2)/p_C(1)*a_C(2)/a_C(1)
       isreal = .true.
    case (22)
       C0 = a_R/a_C(2)*dexp(0.5*RHS1)
       CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K0 = Up+psi(2)/(delta+1.)**((ga+1.)/2./g1)*C0**2*p_C(2)/p_R/CR**(1./ga)
       F = C0**2*Up*CR**(g1/ga)-Up*(1.-delta*Up**2)+ &
            K0*(1.-(delta+1.)*C0**(4./(ga+1.))*(CR*K0/psi(1))**(2.*g1/(ga+1.)))
       isreal = .true.
    case (24)
       C0 = a_R/a_C(1)*dexp(0.5*RHS1)
       CR = -a_S/a_C(1)*p_R/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K0 = 1.+delta*Up**2*CR**(g1/ga)/C0**2
       F = Up-(psi(1)+psi(2))*C0/CR**(g1/2./ga)*(K0/(delta+1.))**((ga+1.)/2./g1)
       isreal = .true.
    case default
       stop ' Wrong case passed to subroutine compute_multivalve_res '
    end select

  end subroutine compute_multivalve_res

  subroutine compute_multivalve_jaco(Up, v_S, p_S, a_S, p_R, a_R, p_C, a_C, a_P, &
       RHS1, RHS2, psi, caso, n_P, ga, dF)
    !
    !
    !  compute_multivalve_jaco is called by:
    !  compute_multivalve_jaco calls the following subroutines 
    !  and functions: none
    !
    implicit none

    integer, intent(in) :: caso,n_P
    real*8, intent(in) :: Up,ga,a_S,v_S,p_S,a_R,p_R,a_P, &
         RHS2,RHS1
    real*8, dimension(2), intent(in) :: p_C,a_C,psi
    real*8, intent(out) :: dF

    real*8 :: C0,C1,C2,CR,K0,K1,K2,rho_T2,a_T2,rho_C1
    real*8 :: delta,g1

    g1 = ga-1.
    delta = 0.5*g1

    select case (caso)
    case (1)
       C1 = -a_S/a_C(2)/ga*p_C(1)/p_S/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = C1*p_C(2)/p_C(1)
       K1 = dsqrt(Up**2+(1.-delta*Up**2)/delta*(1.-C1**(g1/ga)))
       K2 = dsqrt((C2**(g1/ga)*(1.-delta*Up**2)-1.)/delta)
       dF = 1.-psi(1)*p_S/p_C(1)*a_C(2)/a_S*C1**((ga+1.)/ga)*K1- &
            psi(1)*C1**(1./ga)/K1*(Up*C1**(g1/ga)- &
            p_S/p_C(1)*a_C(2)/a_S*(1.-delta*Up**2)*C1**((2.*ga-1.)/ga))+ &
            psi(2)*delta*p_S/p_C(2)*a_C(2)/a_S*C2**((3.*ga-1.)/2./ga)* &
            (1.-delta*Up**2)**((ga-3.)/2./g1)*K2- &
            psi(2)*(ga-3.)/2.*C2**(g1/2./ga)*Up/(1.-delta*Up**2)**((ga+1.)/2./g1)*K2- &
            psi(2)*C2**(g1/2./ga)*(1.-delta*Up**2)**((ga-3.)/2./g1)/K2* &
            (-p_S/p_C(2)*a_C(2)/a_S*(1.-delta*Up**2)*C2**((2.*ga-1)/ga)+Up*C2**(g1/ga))
    case (3)
       C1 = -ga*a_C(1)/a_S*p_S/p_C(1)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(1)/p_C(2)
       K1 = dsqrt((1.-C1**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta)
       K2 = dsqrt((1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta)
       dF = C1**(g1/ga)/(1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)- &
            g1*a_C(1)/a_S*p_S/p_C(1)*Up/C1**(1./ga)/ &
            (1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)-(a_C(1)/a_P)**2.*Up**2.* &
            C1**(g1/ga)/(1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(ga/g1)- &
            g1*psi(1)*Up*K1-psi(1)*(1.-delta*Up**2.)/K1*((a_C(1)/a_P)**2.*Up*C1**(g1/ga)- &
            a_C(1)/a_S*p_S/p_C(1)/C1**(1./ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))- &
            g1*psi(2)*(p_C(2)/p_C(1))**(g1/ga)*a_C(1)/a_C(2)*Up*K2- &
            psi(2)*(p_C(2)/p_C(1))**(g1/ga)*a_C(1)/a_C(2)* &
            ((a_C(2)/a_C(1))**2.-delta*Up**2.)/K2*((a_C(1)/a_P)**2.*Up*C2**(g1/ga)- &
            a_C(1)/a_S*p_S/p_C(2)/C2**(1./ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))
    case (4)
       C0 = (a_R/a_C(2))**2.*dexp(RHS1)
       C1 = -a_S/a_C(2)/ga*p_C(1)/p_S/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = p_C(1)/p_C(2)/C1
       CR = C1*p_R/p_C(1)
       K1 = dsqrt((C0/CR**(g1/ga)+delta*Up**2.-C0*(p_C(1)/p_R)**(g1/ga))/delta)
       K2 = dsqrt((1.-C2**(g1/ga)-delta*Up**2./C0*(p_R/p_C(2))**(g1/ga))/delta)
       dF = 1.-psi(1)*a_C(2)/a_S*p_S/p_C(1)*C1**((ga+1.)/ga)*K1- &
            psi(1)*C1**(1./ga)/K1*(Up-C0*a_C(2)/a_S*p_S/p_R*CR**(1./ga))+ &
            psi(2)*(p_C(2)/p_R)**(g1/ga)*K2* &
            (1.+delta*Up**2./C0*CR**(g1/ga))**((2.-ga)/g1)* &
            (Up*CR**(g1/ga)+delta*Up**2.*a_C(2)/a_S*p_S/p_R*CR**((2.*ga-1.)/ga))- &
            psi(2)*C0*(p_C(2)/p_R)**(g1/ga)* &
            (1.+delta*Up**2./C0*CR**(g1/ga))**(1./g1)/K2* &
            (Up/C0*(p_R/p_C(2))**(g1/ga)-a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga))
    case (6)
       C0 = a_C(1)/a_R*dexp(-0.5*RHS1)
       C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(2)/p_C(1)
       CR = C1*p_R/p_C(1)
       K1 = dsqrt((1.+delta*Up**2*C0**2*CR**(g1/ga)-C1**(g1/ga))/delta)
       K2 = dsqrt((1.+delta*Up**2*C0**2*CR**(g1/ga)-C2**(g1/ga))/delta)
       dF = C0*CR**(delta/ga)+delta*C0*a_C(1)/a_S*p_S/p_R*Up*CR**((3.*ga-1.)/2./ga)- &
            psi(1)*a_C(1)/a_S*p_S/p_C(1)*C1**((ga+1.)/ga)*K1- &
            psi(1)*C1**(1./ga)/K1*(C0**2*Up*CR**(g1/ga)+ &
            delta*C0**2*a_C(1)/a_S*p_S/p_R*Up**2*CR**((2.*ga-1.)/ga)- &
            a_C(1)/a_S*p_S/p_C(1)*C1**((2.*ga-1.)/ga))- &
            psi(2)*a_C(1)/a_S*p_S/p_C(2)*C2**((ga+1.)/ga)*K2- &
            psi(2)*C2**(1./ga)/K2*(C0**2*Up*CR**(g1/ga)+ &
            delta*C0**2*a_C(1)/a_S*p_S/p_R*Up**2*CR**((2.*ga-1.)/ga)- &
            a_C(1)/a_S*p_S/p_C(2)*C2**((2.*ga-1.)/ga))
    case (7)
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K2 = dsqrt((C2**(g1/ga)*(1.-delta*Up**2.)-1.)/delta)
       dF = (1.-delta*Up**2)**(1./g1)-Up**2.*(1.-delta*Up**2)**((2.-ga)/g1)+ &
            psi(2)*p_S/p_C(2)*a_C(2)/a_S*delta*C2**((3.*ga-1.)/2./ga)*dsqrt(1.-delta*Up**2)*K2- &
            psi(2)*C2**(3.*g1/2./ga)*dsqrt(1.-delta*Up**2)/K2* &
            (Up-a_C(2)/a_S*p_S/p_C(2)*C2*(1.-delta*Up**2))- &
            psi(2)*delta*C2**(delta/ga)*Up/dsqrt(1.-delta*Up**2)*K2
    case (8)
       C1 = -a_S/a_C(1)*p_C(1)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = C1*p_C(2)/p_C(1)
       K2 = dsqrt((1.-C2**(g1/ga)*(1.-delta*Up**2))/delta)
       dF = 1.+psi(1)/(delta+1.)**((ga+1.)/g1/2.)*(-g1*Up*C1+ &
            ga*p_S/p_C(1)*a_C(1)/a_S*(1.-delta*Up**2)*C1**2.)- &
            psi(2)*p_S/p_C(2)*a_C(1)/a_S*C2**((ga+1.)/ga)*K2- &
            psi(2)/K2*(C2*Up-p_S/p_C(2)*a_C(1)/a_S*C2**2.*(1.-delta*Up**2))
    case (9)
       C2 = -ga*a_C(1)/a_S*p_S/p_C(2)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K2 = dsqrt((1.-C2**(g1/ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.))/delta)
       dF = C2-ga*a_C(1)/a_S*p_S/p_C(2)*Up-g1*psi(1)/(delta+1.)**((ga+1.)/2./g1)* &
            p_C(1)/p_C(2)*Up-psi(2)*a_C(2)/a_S*p_S/p_C(2)/C2**(g1/ga)* &
            (1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)* &
            (1.-delta*Up**2.*(a_C(1)/a_C(2))**2.)*K2+psi(2)*a_C(2)/a_C(1)*C2**(1./ga)* &
            (1.+delta*Up**2.*(a_C(1)/a_P)**2.)**((2.-ga)/g1)*Up*(a_C(1)/a_P)**2.* &
            (1.-delta*Up**2.*(a_C(1)/a_C(2))**2.)*K2-g1*psi(2)*a_C(1)/a_C(2)* &
            C2**(1./ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)*Up*K2- &
            psi(2)*a_C(2)/a_C(1)*C2**(1./ga)*(1.+delta*Up**2.*(a_C(1)/a_P)**2.)**(1./g1)* &
            (1.-delta*Up**2.*(a_C(1)/a_C(2))**2.)/K2* &
            (Up*(a_C(1)/a_P)**2.*C2**(g1/ga)-a_C(1)/a_S*p_S/p_C(2)/C2**(1./ga)* &
            (1.+delta*Up**2.*(a_C(1)/a_P)**2.))
    case (10)
       C0 = a_R/a_C(2)*dexp(0.5*RHS1)
       CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       C2 = CR*p_R/p_C(2)
       K2 = dsqrt((1.-C2**(g1/ga)-delta*Up**2./C0**2.*(p_R/p_C(2))**(g1/ga))/delta)
       K1 = 1.+delta*Up**2./C0**2./CR**(g1/ga)
       dF = 1.-psi(1)/(1.+delta)**((ga+1.)/2./g1)*C0* &
            (-delta*a_C(2)/a_S*p_S/p_R/CR**((ga+1.)/2./ga)*K1**((ga+1.)/2./g1)+ &
            0.5*(ga+1.)*CR**(delta/ga)*K1**((3.-ga)/2./g1)*(Up/C0**2./CR**(g1/ga)+ &
            delta*Up**2./C0**2.*a_C(2)/a_S*p_S/p_R/CR**((2.*ga-1.)/ga)))+ &
            psi(2)*C0**2.*(p_C(2)/p_R)**(g1/ga)*K1**((2.-ga)/g1)* &
            (Up/C0**2./CR**(g1/ga)+ &
            delta*Up**2./C0**2.*a_C(2)/a_S*p_S/p_R/CR**((2.*ga-1.)/ga))*K2- &
            psi(2)*C0**2.*(p_C(2)/p_R)**(g1/ga)*K1**(1./g1)/K2* &
            (Up/C0**2.*(p_R/p_C(2))**(g1/ga)-a_C(2)/a_S*p_S/p_C(2)/C2**(1./ga))
    case (11)
       C0 = (a_R/a_C(1))**2.*dexp(RHS1)
       C2 = -a_S/a_C(1)*p_C(2)/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       CR = C2*p_R/p_C(2)
       K2 = dsqrt((C0/CR**(g1/ga)+delta*Up**2.-C0*(p_C(2)/p_R)**(g1/ga))/delta)
       dF = 1.+psi(1)/(1.+delta)**((ga+1.)/2./g1)*C0*p_C(1)/p_R* &
            a_C(1)/a_S*p_S/p_R*CR**((ga+1.)/ga)- &
            psi(2)*a_C(1)/a_S*p_S/p_C(2)*C2**((ga+1.)/ga)*K2- &
            psi(2)*C2**(1./ga)/K2*(Up-C0*a_C(1)/a_S*p_S/p_R*CR**(1./ga))
    case (12)
       C0 = a_R/a_C(1)*dexp(0.5*RHS1)
       a_T2 = a_R*(p_C(2)/p_R)**(g1/2./ga)*dexp(0.5*RHS1)
       CR = -a_S/a_C(1)*p_R/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       C2 = CR*p_C(2)/p_R
       K0 = 1./(delta+1.)*(1.+delta/C0**2*CR**(g1/ga)*Up**2)
       K1 = dsqrt((C0**2/CR**(g1/ga)+delta*Up**2-(a_T2/a_C(1))**2)/delta)
       dF = 1.+psi(1)*C0*delta*p_S/p_R*a_C(1)/a_S*CR**((ga+1.)/2./ga)*K0**((ga+1.)/2./g1)- &
            psi(1)/C0*CR**(g1/2./ga)*K0**((3.-ga)/2./g1)*Up*(1.+delta*p_S/p_R*a_C(1)/a_S*CR*Up)- &
            psi(2)*p_S/p_C(2)*a_C(1)/a_S*C2**((ga+1.)/ga)*K1- &
            psi(2)*C2**(1./ga)*(Up-C0**2*p_S/p_R*a_C(1)/a_S*CR**(1./ga))/K1
    case (19)
       C2 = -a_S/a_C(2)*p_C(2)/p_S/ga/(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       dF = (1.-delta*Up**2)**(1./g1)-Up**2*(1.-delta*Up**2)**((2.-ga)/g1)- &
            psi(2)/(delta+1.)**((ga+1.)/2./g1)*C2*(1.-delta*Up**2)**(1./g1)* &
            ga*(Up-C2*(1.-delta*Up**2)*p_S/p_C(2)*a_C(2)/a_S)
    case (21)
       C1 = -ga*a_C(1)/a_S*p_S/p_C(1)*(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       dF = (delta+1.)**((ga+1.)/2./g1)*(C1-ga*p_S/p_C(1)*a_C(1)/a_S*Up)- &
            2*delta*Up*(psi(1)+psi(2)*p_C(2)/p_C(1)*a_C(1)/a_C(2))
    case (22)
       C0 = a_R/a_C(2)*dexp(0.5*RHS1)
       CR = -ga*a_C(2)/a_S*p_S/p_R*(Up-a_S/a_C(2)/ga-n_P*v_S/a_C(2)+RHS2/a_C(2))
       K0 = Up+psi(2)/(delta+1.)**((ga+1.)/2./g1)*C0**2*p_C(2)/p_R/CR**(1./ga)
       dF = C0**2*CR**(-1./ga)*(CR-g1*Up*p_S/p_R*a_C(2)/a_S)+3.*delta*Up**2- &
            (delta+1.)*C0**(4./(ga+1.))*CR**((ga-3.)/(ga+1.))*(K0/psi(1))**(2.*g1/(ga+1.))* &
            (-2.*g1*ga/(ga+1.)*p_S/p_R*a_C(2)/a_S*K0+ &
            (3.*ga-1)/(ga+1.)*CR*(1+psi(2)/(delta+1.)**((ga+1.)/2./g1)*C0**2* &
            p_C(2)/p_R*p_S/p_R*a_C(2)/a_S*CR**(-(ga+1.)/ga)))+ &
            psi(2)/(delta+1.)**((ga+1.)/2./g1)*C0**2*p_C(2)/p_R* &
            p_S/p_R*a_C(2)/a_S*CR**(-(ga+1.)/ga)
    case (24)
       C0 = a_R/a_C(1)*dexp(0.5*RHS1)
       CR = -a_S/a_C(1)*p_R/p_S/ga/(Up-a_S/a_C(1)/ga-n_P*v_S/a_C(1)+RHS2/a_C(1))
       K0 = 1.+delta*Up**2*CR**(g1/ga)/C0**2
       dF = 1.-(psi(1)+psi(2))/(delta+1.)**((ga+1.)/2./g1)*C0* &
            (-ga**2/delta*p_S/p_R*a_C(1)/a_S*CR**((ga+1.)/2./ga)*K0**((ga+1.)/2./g1)+ &
            (ga+1.)/2./C0**2*Up*(K0/CR)**((3.-ga)/2./g1)*(1.+ &
            delta*p_S/p_R*a_C(1)/a_S*Up*CR))
    end select

  end subroutine compute_multivalve_jaco

  subroutine get_multivalve_case(U_P, U_T, U_C, n_P, ga, caso)
    !
    !  Given the states at the pipe end (U_P) and the ones at the 
    !  nozzles' throats (U_T), this subroutine computes the corresponding 
    !  multivalve case
    !
    !  get_multivalve_case is called by: solve_multivalve
    !  get_multivalve_case calls the following subroutines and functions:
    !  none
    !
    implicit none

    real*8, dimension(3), intent(in) :: U_P
    real*8, dimension(3,2), intent(in) :: U_T, U_C
    real*8, intent(in) :: ga
    integer, intent(out) :: caso, n_P

    real*8 :: delta,rho_P,v_P,p_P,a_P,p_0P,tol,tol_vel
    real*8, dimension(2) :: rho_T,v_T,p_T,a_T,p_C
    logical :: subsonic_T1, subsonic_T2

    tol     = 0*1d-4
    tol_vel = 0*1d-2
    delta = 0.5*(ga-1.)
    caso  = 0

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_T = U_T(1,:)
    v_T   = U_T(2,:)
    p_T   = U_T(3,:)

    p_C = U_C(2,:)

    a_P = dsqrt(ga*p_P/rho_P)
    a_T(1) = dsqrt(ga*p_T(1)/rho_T(1))
    a_T(2) = dsqrt(ga*p_T(2)/rho_T(2))

    p_0P = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/(ga-1.))

    if(-n_P*v_T(1) > 0.) then
       subsonic_T1 = p_C(1)/p_0P < (delta+1.)**(ga/(ga-1.))
    else
       subsonic_T1 = p_0P/p_C(1) < (delta+1.)**(ga/(ga-1.))
    end if
    if(-n_P*v_T(2) > 0.) then
       subsonic_T2 = p_C(2)/p_0P < (delta+1.)**(ga/(ga-1.))
    else
       subsonic_T2 = p_0P/p_C(2) < (delta+1.)**(ga/(ga-1.))
    end if

    ! n_P is the external unit normal to the pipe
    if(n_P*v_P > -tol_vel) then
       if(-n_P*v_T(1) > 0.) then
          if(-n_P*v_T(2) > 0.) then
             stop ' IMPOSIBLE CASE: uP* > 0 & uT1* > 0 & uT2* > 0 '
          else ! u_T2*n_T2 < 0
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 5
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 17
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 11
             else
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 23
             end if
          end if
       else ! u_T1*n_T1 < 0
          if(-n_P*v_T(2) > 0.) then
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 4
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 16
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 10
             else
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 22
             end if
          else ! u_T2*n_T2 < 0
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 6
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 18
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 12
             else
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 24
             end if
          end if
       end if
    else ! u_P*n_P < 0
       if(-n_P*v_T(1) > 0.) then
          if(-n_P*v_T(2) > 0.) then
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 3
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 15
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 9
             else
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 21
             end if
          else ! u_T2*n_T2 < 0
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 2
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 14
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 8
             else
                if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                     caso = 20
             end if
          end if
       else ! u_T1*n_T1 < 0
          if(-n_P*v_T(2) > 0.) then
             if(subsonic_T1.and.subsonic_T2) then ! subsonic flow at the throats
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 1
             elseif(subsonic_T1) then ! sonic flow at the throat 2
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 13
             elseif(subsonic_T2) then ! sonic flow at the throat 1
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 7
             else
                if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                     caso = 19
             end if
          else ! u_T2*n_T2 < 0
             stop ' IMPOSIBLE CASE: uP* < 0 & uT1* < 0 & uT2* < 0 '
          end if
       end if
    end if

  end subroutine get_multivalve_case

  subroutine multivalve_cases(caso, CASOS, p_C)
    !
    !  From the list of 24 possible cases, selects the allowed cases
    !  taking into acount the case solved at the last time step
    !  
    !  multivalve_cases is called by: solve_multivalve
    !  multivalve_cases calls the following subroutines and functions:
    !  none
    !
    implicit none

    integer, intent(in) :: caso
    real*8, dimension(2), intent(in) :: p_C
    integer, dimension(24), intent(out) :: CASOS

    real*8 :: tol

    tol = 0*1d-4

    select case(caso)
    case(0)
       CASOS = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
            19,20,21,22,23,24/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,0,0,0,0,0,0,0,0/)
       end if
    case(1)
       CASOS = (/1,4,3,6,2,5,7,10,13,16,15,19,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/3,6,2,5,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/1,4,3,6,7,10,13,16,15,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(2)
       CASOS = (/2,5,3,6,1,4,8,11,9,14,17,20,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/2,5,3,6,8,11,9,14,17,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/3,6,1,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(3)
       CASOS = (/3,2,1,5,4,6,9,8,15,13,21,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/3,2,5,6,9,8,15,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/3,1,4,6,9,15,13,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(4)
       CASOS = (/4,1,6,3,5,2,10,7,12,16,13,22,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/6,3,5,2,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/4,1,6,3,10,7,12,16,13,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(5)
       CASOS = (/5,2,6,3,4,1,11,8,17,14,18,23,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/5,2,6,3,11,8,17,14,18,23,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/6,3,4,1,18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(6)
       CASOS = (/6,4,5,1,2,3,12,10,18,17,24,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/6,5,2,3,12,18,17,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       elseif(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/6,4,1,3,12,10,18,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(7)
       CASOS = (/7,1,10,4,13,16,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(8)
       CASOS = (/8,2,9,5,3,14,17,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/9,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(9)
       CASOS = (/9,3,8,11,2,15,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/9,3,15,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(10)
       CASOS = (/10,4,7,12,1,6,16,13,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/12,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(11)
       CASOS = (/11,5,9,2,17,14,23,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(12)
       CASOS = (/12,6,10,4,18,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/12,6,18,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(13)
       CASOS = (/13,1,15,4,3,7,10,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/15,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(14)
       CASOS = (/14,2,17,5,8,11,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(15)
       CASOS = (/15,3,13,1,9,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/15,3,9,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(16)
       CASOS = (/16,4,1,10,7,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(17)
       CASOS = (/17,5,14,18,2,6,11,8,23,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/18,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(18)
       CASOS = (/18,6,17,5,12,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/18,6,12,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(19)
       CASOS = (/19,13,7,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(20)
       CASOS = (/20,14,8,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(21)
       CASOS = (/21,15,9,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    case(22)
       CASOS = (/22,16,10,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(1)/p_C(2).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(23)
       CASOS = (/23,17,11,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       if(p_C(2)/p_C(1).gt.(1.-tol)) then
          CASOS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       end if
    case(24)
       CASOS = (/24,18,12,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    end select

  end subroutine multivalve_cases

  subroutine check_multivalve_case(U_P, U_T, U_C, n_P, ga, caso, caso_ok)
    !
    !  Given the states at the pipe end (U_P) and the ones at the 
    !  nozzles' throats (U_T), this subroutine checks the corresponding 
    !  multivalve case
    !
    !  check_multivalve_case is called by: solve_multivalve
    !  check_multivalve_case calls the following subroutines and functions:
    !  none
    !
    implicit none

    integer, intent(in) :: caso, n_P
    real*8, dimension(3), intent(in) :: U_P
    real*8, dimension(3,2), intent(in) :: U_T, U_C
    real*8, intent(in) :: ga
    logical, intent(out) :: caso_ok

    real*8 :: delta,rho_P,v_P,p_P,a_P,p_0P,tol,tol_vel
    real*8, dimension(2) :: rho_T,v_T,p_T,a_T,p_C
    logical :: subsonic_T1, subsonic_T2

    tol     = 1.0d-4
    tol_vel = 0*1d-2
    delta = 0.5*(ga-1.)

    caso_ok = .false.

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_T = U_T(1,:)
    v_T   = U_T(2,:)
    p_T   = U_T(3,:)

    p_C = U_C(2,:)

    a_P = dsqrt(ga*p_P/rho_P)
    a_T(1) = dsqrt(ga*p_T(1)/rho_T(1))
    a_T(2) = dsqrt(ga*p_T(2)/rho_T(2))

    p_0P = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/(ga-1.))

    if(-n_P*v_T(1) > 0.) then
       subsonic_T1 = p_C(1)/p_0P < (delta+1.)**(ga/(ga-1.))
    else
       subsonic_T1 = p_0P/p_C(1) < (delta+1.)**(ga/(ga-1.))
    end if
    if(-n_P*v_T(2) > 0.) then
       subsonic_T2 = p_C(2)/p_0P < (delta+1.)**(ga/(ga-1.))
    else
       subsonic_T2 = p_0P/p_C(2) < (delta+1.)**(ga/(ga-1.))
    end if

    ! n_P is the external unit normal to the pipe
    select case(caso)
    case(1)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                  caso_ok = .true.
          end if
       end if
    case(2)
       if((n_P*v_P < 1.0d-3).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                  caso_ok = .true.
          end if
       end if
    case(3)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                  caso_ok = .true.
          end if
       end if
    case(4)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) <= (1.+tol)) &
                  caso_ok = .true.
          end if
       end if
    case(5)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) <= (1.+tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                  caso_ok = .true.
          end if
       end if
    case(6)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol) .and. p_0P/p_C(2) >= (1.-tol)) &
                  caso_ok = .true.
          end if
       end if
    case(7)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) <= (1.+tol)) caso_ok = .true.
          end if
       end if
    case(8)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(9)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) <= (1.+tol)) caso_ok = .true.
          end if
       end if
    case(10)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) <= (1.+tol)) caso_ok = .true.
          end if
       end if
    case(11)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(12)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and.subsonic_T2) then
             if(p_0P/p_C(2) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(13)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(14)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             if(p_0P/p_C(1) <= (1.+tol)) caso_ok = .true.
          end if
       end if
    case(15)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             ! if(p_0P/p_C(1) <= (1.+tol)) 
             caso_ok = .true.
          end if
       end if
    case(16)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(17)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             if(p_0P/p_C(1) <= (1.+tol)) caso_ok = .true.
          end if
       end if
    case(18)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(subsonic_T1.and. .not.subsonic_T2) then
             if(p_0P/p_C(1) >= (1.-tol)) caso_ok = .true.
          end if
       end if
    case(19)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    case(20)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    case(21)
       if((n_P*v_P < tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    case(22)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) > 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    case(23)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) > 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    case(24)
       if((n_P*v_P > -tol_vel).and.(-n_P*v_T(1) < 0.).and. &
            (-n_P*v_T(2) < 0.)) then
          if(.not.subsonic_T1.and. .not.subsonic_T2) caso_ok = .true.
       end if
    end select

  end subroutine check_multivalve_case

  subroutine solve_multivalve_explicit(Ucyl, Uref, type_val, Area_val, Area_pipe, &
       RHS, ga, Upipe, Uthroat, solved_case)
    !
    !  Solves the non-linear equation system for the multi-valve model.
    !  The solver is based on the selection of cases, each of them
    !  corresponding to a particular set of equations. The cylinders'
    !  states are assumed constant for this solver.
    !
    !  solve_multivalve is called by: solve_cylinder
    !  solve_multivalve calls the following subroutines and 
    !  functions: solve_multivalve_eqcases, get_multivalve_case
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: ga,Area_pipe
    real*8, dimension(2), intent(in) :: Area_val
    real*8, dimension(3), intent(in) :: RHS
    real*8, dimension(3,2), intent(in) :: Ucyl,Uref
    real*8, dimension(3), intent(inout) :: Upipe
    real*8, dimension(3,2), intent(inout) :: Uthroat
    integer, intent(inout) :: solved_case

    integer :: try_case, max_iter, iter, n_P, i
    logical :: conv_flag, case_flag, do_cases
    real*8 :: g1, delta, rp_cr, rp, Up
    real*8 :: rho_P, v_P, p_P, a_P, rho_0P, p_0P, a_0P, &
         rho_S, p_S, a_S
    real*8 :: rho_T1, v_T1, p_T1, a_T1, rho_T2, v_T2, p_T2, a_T2
    real*8, dimension(2) :: psi, rho_C, p_C, a_C, mdot, Edot
    real*8, dimension(3) :: U_P
    real*8, dimension(3,2) :: U_T

    n_P = type_val

    U_T = Uthroat
    psi(1) = Area_val(1)/Area_pipe
    psi(2) = Area_val(2)/Area_pipe

    g1 = ga-1.
    delta = g1/2.
    rp_cr = 1./(1.+delta)**(ga/g1)

    do i=1,2
       rho_C(i) = Ucyl(1,i)
       p_C(i)   = Ucyl(2,i)
       a_C(i)   = dsqrt(ga*p_C(i)/rho_C(i))
    end do

    do i=1,1
       U_P = Upipe
       
       rho_P = U_P(1)
       v_P   = U_P(2)
       p_P   = U_P(3)
       
       a_P = dsqrt(ga*p_P/rho_P)
       a_0P   = a_P*dsqrt(1.+delta*(v_P/a_P)**2.)
       p_0P   = p_P*(a_0P/a_P)**(ga/delta)
       rho_0P = rho_P*(a_0P/a_P)**(2./g1)
       
       mdot = (/ 0., 0./)
       Edot = (/ 0., 0./)
       
       ! Computes the state at throat of nozzle 1
       if(psi(1).gt.1d-6) then
          rp = p_C(1)/p_0P
          if(rp.le.rp_cr) then
             ! Sonic inflow to chamber 1
             p_T1   = p_0P*rp_cr
             rho_T1 = rho_0P/(1.+delta)**(1./g1)
             v_T1   = a_0P/dsqrt(1.+delta)*n_P
             a_T1   = dabs(v_T1)
          elseif(rp.le.1.) then
             ! Subsonic inflow to chamber 1
             p_T1   = p_C(1)
             rho_T1 = rho_0P*rp**(1./ga)
             a_T1   = a_0P*rp**(delta/ga)
             v_T1   = dsqrt((a_0P**2.-a_T1**2.)/delta)*n_P
          elseif(rp.lt.1./rp_cr) then
             ! Subsonic outflow from chamber 1
             p_T1   = p_0P
             rho_T1 = rho_C(1)/rp**(1./ga)
             a_T1   = a_C(1)/rp**(delta/ga)
             v_T1   = -dsqrt((a_C(1)**2.-a_T1**2.)/delta)*n_P
          else
             ! Sonic outflow from chamber 1
             p_T1   = p_C(1)*rp_cr
             rho_T1 = rho_C(1)/(1.+delta)**(1./g1)
             v_T1   = -a_C(1)/dsqrt(1.+delta)*n_P
             a_T1   = dabs(v_T1)
          end if
          mdot(1) = -rho_T1*psi(1)*v_T1*n_P
          Edot(1) = mdot(1)*(a_T1**2.+delta*v_T1**2.)
       end if
       
       ! Computes the state at throat of nozzle 2
       if(psi(2).gt.1d-6) then
          rp = p_C(2)/p_0P
          if(rp.le.rp_cr) then
             ! Sonic inflow to chamber 2
             p_T2   = p_0P*rp_cr
             rho_T2 = rho_0P/(1.+delta)**(1./g1)
             v_T2   = a_0P/dsqrt(1.+delta)*n_P
             a_T2   = dabs(v_T2)
          elseif(rp.le.1.) then
             ! Subsonic inflow to chamber 2
             p_T2   = p_C(2)
             rho_T2 = rho_0P*rp**(1./ga)
             a_T2   = a_0P*rp**(delta/ga)
             v_T2   = dsqrt((a_0P**2.-a_T2**2.)/delta)*n_P
          elseif(rp.lt.1./rp_cr) then
             ! Subsonic outflow from chamber 2
             p_T2   = p_0P
             rho_T2 = rho_C(2)/rp**(1./ga)
             a_T2   = a_C(2)/rp**(delta/ga)
             v_T2   = -dsqrt((a_C(2)**2.-a_T2**2.)/delta)*n_P
          else
             ! Sonic outflow from chamber 2
             p_T2   = p_C(2)*rp_cr
             rho_T2 = rho_C(2)/(1.+delta)**(1./g1)
             v_T2   = -a_C(2)/dsqrt(1.+delta)*n_P
             a_T2   = dabs(v_T2)
          end if
          mdot(2) = -rho_T2*psi(2)*v_T2*n_P
          Edot(2) = mdot(2)*(a_T2**2.+delta*v_T2**2.)
       end if
       
       Uthroat(1,1) = rho_T1
       Uthroat(2,1) = v_T1
       Uthroat(3,1) = p_T1
       
       Uthroat(1,2) = rho_T2
       Uthroat(2,2) = v_T2
       Uthroat(3,2) = p_T2
       
       ! Computes the state at the point P
       rho_S = Uref(1,2)
       p_S   = Uref(3,2)
       a_S    = dsqrt(ga*p_S/rho_S)
       
       Up = v_P*n_P/a_S
       
       do_cases = .true.
       conv_flag = .false.
       iter = 0
       do while (do_cases)
          iter = iter+1
          if(Up.le.0.) then
             try_case = 1
          else
             try_case = 2
          end if
          call solve_xmultivalve_eqcases(Up, Uref(:,2), Uref(:,1), RHS, Ucyl, &
               n_P, psi, ga, try_case, mdot, Edot, Upipe, case_flag)
          if(case_flag) then
             conv_flag = .true.
             solved_case = try_case
          else
             Up = -v_P*n_P/a_S
          end if
          if(conv_flag .or. iter==2) do_cases = .false.
       end do

       Upipe = U_P
       if(.not.conv_flag) then
          write(*,*) 'WARNING! - no solution founded at multivalve'
          Upipe = U_P
       end if
    end do
    
  end subroutine solve_multivalve_explicit

  subroutine solve_xmultivalve_eqcases(Up, U_S, U_R, RHS, U_C, n_P, psi, ga, &
       caso, mdot, Edot, U_P, ok_flag)
    !
    !  Solves the non-linear equation arising from the multivalve model for
    !  the specified case
    !
    !  solve_multivalve_eqcases is called by:
    !  solve_multivalve_eqcases calls the following subroutines and 
    !  functions: none
    !
    implicit none

    integer, intent(in) :: caso, n_P
    real*8, intent(in) :: ga
    real*8, dimension(2), intent(in) :: psi, mdot, Edot
    real*8, dimension(3), intent(in) :: U_S, U_R, RHS
    real*8, dimension(3,2), intent(in) :: U_C
    real*8, intent(inout) :: Up
    real*8, dimension(3), intent(out) :: U_P
    logical, intent(out) :: ok_flag

    integer :: max_ils,max_iter,iter,ils
    real*8 :: rho_S,v_S,p_S,a_S,rho_R,v_R,p_R,a_R, &
         rho_P,v_P,p_P
    real*8 :: Up_new,RHS1,RHS2,delta,F,dF,tol,dU,nu,relax,g1
    real*8 :: A,B,C,D,tol_ini
    real*8, dimension(2) :: p_C,a_C,rho_C
    logical :: do_newton,do_ls,F_isreal
    real*8 :: Fnew

    g1 = ga-1.
    delta = 0.5*g1

    tol = 1d-8
    max_iter = 20
    max_ils  = 10
    nu = 2.0

    tol_ini = 1d-6

    F = huge(0.0d0)

    ok_flag   = .false.
    do_newton = .false.

    rho_C = U_C(1,:)
    p_C   = U_C(2,:)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_R = U_R(1)
    v_R   = U_R(2)
    p_R   = U_R(3)

    a_S = dsqrt(ga*p_S/rho_S)
    a_R = dsqrt(ga*p_R/rho_R)
    a_C(1) = dsqrt(ga*p_C(1)/rho_C(1))
    a_C(2) = dsqrt(ga*p_C(2)/rho_C(2))

    RHS1 = RHS(1)
    if(n_P < 0) then
       RHS2 = RHS(2)
    elseif(n_P > 0) then
       RHS2 = RHS(3)
    end if

    if(.true.) then
    select case(caso)
    case(1)
       A = ga*a_S+delta*sum(mdot)/rho_S
       B = -a_S-ga*v_S*n_P+ga*RHS2
       C = -sum(Edot)/ga/p_S

       D = B**2.-4.*A*C
       if(D.ge.0.) then
          Up = (-B-dsqrt(D))/2./A
          if(Up.le.0.) then
             ok_flag = .true.
          else
             Up = (-B+dsqrt(D))/2./A
             if(Up.le.0.) ok_flag = .true.
          end if
       end if
    case(2)
       A = 1./ga+n_P*v_S/a_S-RHS2/a_S
       if(Up.ge.A) Up = Up-tol_ini
       if(Up > 0.) do_newton=.true.

       ! Newton loop
       iter = 0
       do while(do_newton)
          do_ls = .true.
          iter  = iter+1
          relax = 1.0
          C = -ga*p_S/p_R*(Up-1./ga-n_P*v_S/a_S+RHS2/a_S)
          
          F  = ga*p_R*a_S/a_R**2.*C**(1./ga)*Up+sum(mdot)
          dF = ga*p_R*a_S/a_R**2.*C**(1./ga)- &
               ga*p_S*a_S/a_R**2./C**(g1/ga)*Up

          dU = -F/dF
          ! line-search loop
          ils = 0
          do while(do_ls)
             ils = ils+1
             Up_new = Up+relax*dU
             C = -ga*p_S/p_R*(Up_new-1./ga-n_P*v_S/a_S+RHS2/a_S)
             if(C.gt.0.) then
                F_isreal = .true.
                Fnew = ga*p_R*a_S/a_R**2.*C**(1./ga)*Up_new+sum(mdot)
             else
                F_isreal = .false.
                Fnew = huge(0.0d0)
             end if
             if(F_isreal .and. dabs(Fnew).lt.dabs(F)) then
                F  = Fnew
                Up = Up_new
                exit
             else
                relax = relax/nu
             end if
             if(ils.ge.max_ils) exit
          end do ! ends line-search loop

          if(dabs(F).lt.tol .and. Up.gt.0.) then
             ok_flag =.true.
             exit
          end if
          if(ils.ge.max_ils .or. iter.gt.max_iter) exit
       end do ! ends newton loop
    end select
    end if

    if(ok_flag) then
       v_P = a_S*Up*n_P
       p_P = -p_S*ga*(Up-1./ga-n_P*v_S/a_S+RHS2/a_S)
       
       if(caso==1) then
          rho_P = -sum(mdot)/v_P*n_P
       elseif(caso==2) then
          rho_P = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS1)
       endif

       U_P(1) = rho_P
       U_P(2) = v_P
       U_P(3) = p_P
    end if

  end subroutine solve_xmultivalve_eqcases

end module def_multivalve
