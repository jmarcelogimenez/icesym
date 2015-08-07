module def_valve

  use def_simulator
  use utilities
  use gasdyn_utils, only : compute_cp

  type :: valve
     integer :: Nval, type_dat, tube, valve_model
     integer :: solved_case
     real*8, dimension(:,:), pointer  :: Cd, Lv
     real*8 :: angle_VO, angle_VC, Dv, Lvmax, &
          Area, dArea, Area_max
     real*8 :: Area_tube, twall_tube, dAreax_tube, dx_tube
     real*8 :: mass_res_port, mass_flow_factor
     real*8, dimension(3) :: state, state_old, state_pipe_old
     real*8, dimension(3,2) :: state_ref
  end type valve

contains

  subroutine solve_valve(Ucyl, Uref, type_val, Area_val, Area_pipe, &
       ga, R_gas, Upipe, Uthroat)
    !
    !
    !  solve_valve is called by: solve_cylinder
    !  solve_valve calls the following subroutines and functions: 
    !    solve_valve_equation
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: Area_val,Area_pipe,ga,R_gas
    real*8, dimension(3), intent(in) :: Ucyl,Uref
    real*8, dimension(3), intent(out) :: Upipe,Uthroat

    real*8 :: tol,cp,g1,delta,psi,rp,rp_cr
    real*8 :: rho_P,v_P,p_P,T_P,a_P,rho_T,v_T,p_T,a_T
    real*8 :: rho0,v0,p0,a0,T0,a_ref,p_C
    real*8 :: A,A_0,At,U,Ut,rA,v1
    real*8, dimension(3) :: param
    logical :: noflow,inflow

    tol = 1.0d-10

    g1    = ga-1.
    delta = .5*g1

    p_P = Uref(3)
    p_C = Ucyl(2)
    rp  = (p_P/p_C)**(delta/ga)

    psi = Area_val/Area_pipe
    if(psi.gt.1.) &
         write(*,*) ' *** WARNING *** Area_throat > Area_pipe the nozzle flow convergence could fail '

    noflow = .false.
    if(rp.lt.1.0d0-tol) then
       inflow = .false.
    elseif(rp.gt.1.0d0+tol) then
       inflow = .true.
    else
       noflow = .true.
    endif
    if(psi.lt.tol) noflow = .true.

    param = 0.
    
    if(noflow) then
       rho0 = Uref(1)
       v0   = Uref(2)
       p0   = Uref(3)
       a0   = dsqrt(ga*p0/rho0)
       T0   = p0/rho0/R_gas

       cp    = compute_cp(T0, 0, R_gas)
       T_P   = T0+0.5*v0**2/cp
       p_P   = p0*(T0/T_P)**(-ga/g1)
       rho_P = p_P/R_gas/T_P
       v_P   = 0.

       rho_T = rho_P
       v_T   = v_P
       p_T   = p_P
    else
       if(inflow) then
          ! inflow (pipe to cylinder)
          rho0 = Uref(1)
          v0   = Uref(2)
          p0   = Uref(3)
          p0   = p0 + 0.5*rho0*v0**2  ! stagnation pressure
          a0   = dsqrt(ga*p0/rho0)    ! stagnation speed of sound

          a_ref = a0*(p_C/p0)**(delta/ga)
          A_0   = a0/a_ref
          
          param(1) = A_0
          U        = 0.5
          call solve_valve_equation(1, U, psi, rp, ga, param)

          A  = dsqrt(A_0**2-delta*U**2)
          At = 1.
          Ut = U*(A/At)**(1./delta)/psi

          if(Ut.ge.1.) then
             rA = 1.1
             call solve_valve_equation(2, rA, psi, rp, ga, param)

             v1  = psi**2/rA**((ga+1.)/delta)
             U   = dsqrt(v1*A0**2/(1.+delta*v1))
             A   = dsqrt(A_0**2-delta*U**2)
             At  = A/rA
             Ut  = At
          endif

          a_P = A*a_ref
          a_T = At*a_ref

          p_P   = Uref(3)
          rho_P = ga*p_P/a_P**2
          v_P   = type_val*a_ref*U

          v_T   = type_val*Ut*a_ref
          p_T   = p0*(a_T/a0)**(ga/delta)
          rho_T = ga*p_T/a_T**2
       else
          ! outflow (cylinder to pipe)
          rho0 = Ucyl(1)
          p0   = Ucyl(2)
          a0   = dsqrt(ga*p0/rho0)

          U = 0.1
          call solve_valve_equation(3, U, psi, rp, ga, param)

          a_P = a0*dsqrt(1.-delta*U**2)
          a_T = rp*a0

          p_P   = Uref(3)
          rho_P = ga*p_P/a_P**2
          v_P   = -type_val*U*a0

          p_T   = p_P
          rho_T = ga*p_T/a_T**2
          v_T   = -type_val*(1./psi*(a_T/a_P)**2*a0*U)

          if(dabs(v_T)/a_T .ge. 1.) then
             ! Sonic outflow from cylinder
             rp_cr = dsqrt(2.0d0/(ga+1.0d0))

             if(rp/rp_cr.gt.1.2) &
                  stop ' RES_16 | OUTFLOW | rp/rp_cr > 1.2 '
             if(rp/rp_cr .lt. (1.-tol)) then
                ! Sonic flow through valve throat
                v1 = rp**(ga/delta)/psi/((2.0/(1.0+ga))**((ga+1.0)/2.0/g1))
                U  = (-v1+dsqrt(v1**2+2.0*g1))/g1

                a_P = a0*dsqrt(1.0d0-delta*U**2)

                rho_P = ga*p_P/a_P**2
                v_P   = -type_val*U*a0

                p_T   = p0*dsqrt(2.0d0/(ga+1.0d0))
                rho_T = rho0*(p_T/p0)**(1./ga)
                v_T   = 1./psi*(rho_P/rho_T)*v_P
             else
                ! Sonic boundary
                U = 0.5
                call solve_valve_equation(4, U, psi, rp, ga, param)

                a_P = a0*dsqrt(1.0d0-delta*U**2)
                a_T = rp*a0

                rho_P = ga*p_P/a_P**2
                v_P   = -type_val*U*a0
                
                rho_T = ga*p_T/a_T**2
                v_T   = -type_val*(1./psi*(a_T/a_P)**2*a0*U)
             endif

             if(dabs(v_P)/a_P .ge. 1.) then
                U  = dsqrt(2.0d0/(ga+1.0d0))

                v_P   = -type_val*U*a0
                a_P   = v_P
                rho_P = ga*p_P/a_P**2

                ! Computing rho_T
                rho_T    = rho0
                param(1) = 1.0/rho0
                ! param(2) = rho_P*v_P/psi  ! CHEQUEAR
                param(2) = rho_P*v_P*psi
                param(3) = a0
                call solve_valve_equation(5, rho_T, psi, rp, ga, param)

                ! v_T = rho_P*v_P/(psi*rho_T)  ! CHEQUEAR
                v_T = rho_P*v_P*psi/rho_T
                a_T = a0*(rho_T/rho0)**delta
                p_T = a_T**2*rho_T/ga
             endif
          endif
       endif  ! outflow
    endif  ! noflow

    Upipe(1) = rho_P
    Upipe(2) = v_P
    Upipe(3) = p_P

    Uthroat(1) = rho_T
    Uthroat(2) = v_T
    Uthroat(3) = p_T

  end subroutine solve_valve

  subroutine solve_valve_equation(case, U, psi, rp, ga, param)
    !
    !
    !  solve_valve_equation is called by: solve_valve
    !  solve_valve_equation calls the following subroutines and functions: 
    !    valve_res_jac
    !
    implicit none

    integer, intent(in) :: case
    real*8, intent(in) :: psi,rp,ga
    real*8, dimension(3), intent(in) :: param
    real*8, intent(inout) :: U

    integer :: niter,niter_max
    logical :: continue_flag
    real*8 :: dU,R,J,relax,normres,tol
    real*8 :: normres0,dU0

    niter_max = 100
    tol       = 1.0d-6
    relax     = 1.

    niter = 0
    continue_flag = .true.
    do while (continue_flag)
       R = 0.0d0
       J = 0.0d0
       relax = tanh(0.05*(niter+1))
       call valve_res_jac(case, U, psi, rp, ga, param, R, J)
       dU = -relax*R/J
       ! if(isnan(dU)) then
       !    stop
       ! end if
       if(niter.eq.0) then
          normres0 = dabs(R)
          dU0      = dabs(dU)
       end if
       U = U+dU

       normres = dabs(R)
       niter   = niter+1
       if((normres/normres0.lt.tol).or.(niter.gt.niter_max) &
            .or.(dabs(dU)/dU0).lt.tol) &
            continue_flag = .false.
    end do

  end subroutine solve_valve_equation

  subroutine valve_res_jac(case, U, psi, rp, ga, param, res, jac)
    !
    !
    !  valve_res_jac is called by: solve_valve_equation
    !  valve_res_jac calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: case
    real*8, intent(in) :: U,psi,rp,ga
    real*8, dimension(3), intent(in) :: param
    real*8, intent(out) :: res,jac

    real*8 :: g1,delta,v1,v2,a0

    g1    = ga-1.
    delta = .5*g1

    if(case.eq.1) then
       ! ec 3.15
       a0 = param(1)
       v2 = (a0**2-delta*U**2-1.)/delta
       v1 = (a0**2-delta*U**2)**(1./delta)/psi**2
       res = U**2*(v1-1.)-v2
       jac = 2.*U*v1-2.*U**3/psi**2*(a0**2-delta*U**2)**((1.-delta)/delta)
    elseif(case.eq.2) then
       ! ec 3.21
       v1 = (ga+1.)/g1-U**2/delta
       v2 = U**(2./delta)
       res = psi**2-v1*v2
       jac = 2./delta*(U*v2-v1*U**((2.-delta)/delta))
    elseif(case.eq.3) then
       ! ec 4.12
       res = psi/rp*dsqrt((1./rp**2-1.)/delta) - &
            U/(1.-delta*U**2)
       jac = -(1.+delta*U**2)/(1.-delta*U**2)**2
    elseif(case.eq.4) then
       ! ec 4.16
       res = psi-dsqrt(2./(1.+ga))*U/(1.-delta*U**2)
       jac = -dsqrt(2./(1.+ga))*(1.+delta*U**2)/(1.-delta*U**2)**2
    elseif(case.eq.5) then
       ! ec 4.21
       v1 = param(1)
       v2 = param(2)
       a0 = param(3)
       res = a0**2*(1.-(v1*U)**g1)-delta*(v2/U)**2
       jac = 2.*delta*(-a0**2*v1**(2.*delta)*U**(2.*delta-1.)+ &
            v2**2/U**3)
    else
       stop ' Case is wrong @ valve_residue '
    endif

  end subroutine valve_res_jac

  subroutine solve_valve_implicit(Ucyl, Uref, type_val, Area_val, Area_pipe, &
       RHS, ga, Upipe, Uthroat, solved_case)
    !
    !  Solves the non-linear equation system for the valve model.
    !  The solver is based on the selection of cases, each of them
    !  corresponding to a particular set of equations. The cylinder
    !  state is assumed constant for this solver.
    !
    !  solve_valve_implicit is called by: solve_cylinder, solve_tank
    !  solve_valve_implicit calls the following subroutines and 
    !  functions: solve_valve_eqcases, get_valve_case
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: ga, Area_val, Area_pipe
    real*8, dimension(3), intent(in) :: RHS, Ucyl
    real*8, dimension(3,2), intent(in) :: Uref
    real*8, dimension(3), intent(inout) :: Upipe
    real*8, dimension(3), intent(inout) :: Uthroat
    integer, intent(out) :: solved_case

    integer :: old_case, try_case, caso, max_iter, i, icase, iter
    logical :: conv_flag, case_flag, do_cases
    integer, dimension(4) :: CASOS
    real*8 :: delta, n_P, psi, rho_S, v_S, p_S, a_S, &
         rho_R, v_R, p_R, rho_C, T_C, p_C, a_C, a_P, &
         rho_P, v_P, p_P, Up
    real*8, dimension(3) :: U_S, U_R
    real*8, dimension(3) :: U_P, U_T

    delta = 0.5*(ga-1.)
    n_P   = type_val

    U_R = Uref(:,1)
    U_S = Uref(:,2)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_R = U_R(1)
    v_R   = U_R(2)
    p_R   = U_R(3)

    U_P = Upipe
    U_T = Uthroat
    psi = Area_val/Area_pipe

    old_case = 1
    try_case = old_case
    CASOS = (/1,2,3,4/)

    solved_case = 0

    max_iter = 3

    if(psi.gt.1d-6) then
       do_cases = .true.
       rho_C = Ucyl(1)
       p_C   = Ucyl(2)
       T_C   = Ucyl(3)

       a_C = dsqrt(ga*p_C/rho_C)

       conv_flag = .false.
       iter = 0
       do while (do_cases)
          iter = iter+1
          rho_P = U_P(1)
          v_P   = U_P(2)
          p_P   = U_P(3)
          a_P   = dsqrt(ga*p_P/rho_P)

          if (v_P == 0) v_P = n_P
          Up = v_P*n_P/a_C

          call solve_valve_eqcases(Up, U_S, U_R, RHS, Ucyl, n_P, psi, ga, &
               try_case, Upipe, Uthroat, case_flag)

          do i=1, 4
             if (CASOS(i)==try_case) then
                CASOS(i)=0
                exit
             end if
          end do

          if(case_flag) then
             ! Here, we must to check the case
             call get_valve_case(Upipe, Uthroat, Ucyl, n_P, ga, caso)

             if (try_case==caso) then
                conv_flag = .true.
                solved_case = try_case
             else
                icase = 0
                do i=1, 4
                   if (CASOS(i)==caso) then
                      icase=i
                      exit
                   end if
                end do
                if (icase.ne.0) then
                   try_case=icase
                else
                   do i=1, 4
                      if (CASOS(i).ne.0) then
                         icase=i
                         exit
                      end if
                   end do
                   if (icase.ne.0) then
                      try_case = icase
                   else
                      try_case = 0
                   end if
                end if
             end if
          else
             icase = 0
             do i=1, 4
                if (CASOS(i).ne.0) then
                   icase=i
                   exit
                end if
             end do
             if (icase.ne.0) then
                try_case = icase
             else
                try_case = 0
             end if
          end if

          if(conv_flag .or. try_case==0) &
               do_cases = .false.
       end do
       ! For debugging
       ! if(.not.conv_flag .or. any(isnan(Uthroat))) then
       !    write(*,*) 'VALVE SOLVER'
       !    write(*,*) 'n_P: ', n_P, 'psi: ', psi
       !    write(*,*) 'U_S: ', Uref
       !    write(*,*) 'U_C: ', Ucyl
       !    write(*,*) 'RHS: ', RHS
       !    ! stop
       ! end if
    else
       a_S = dsqrt(ga*p_S/rho_S)
       v_P = 0.
       if (type_val.gt.0) then
          p_P = p_S+rho_S*a_S*(v_S-RHS(2))
       else
          p_P = p_S-rho_S*a_S*(v_S+RHS(3))
       end if
       rho_P = rho_S*(p_P/p_S)**(1./ga)*dexp(-RHS(1))

       Upipe(1) = rho_P
       Upipe(2) = v_P
       Upipe(3) = p_P

       Uthroat = Upipe
    end if

  end subroutine solve_valve_implicit

  subroutine valve_flow(nval, type_val, ga, F_v, UT, UP, UC, &
       mdot, hdot)
    !
    !
    !  valve_flow is called by: cylinder_solver
    !  valve_flow calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: nval, type_val
    real*8, intent(in) :: ga
    real*8, dimension(3), intent(in) :: UC
    real*8, dimension(nval), intent(in) :: F_v
    real*8, dimension(3,nval), intent(in) :: UT, UP
    real*8, dimension(nval), intent(out) :: mdot, hdot

    integer :: i
    real*8 :: g1,rho,vel,pre

    g1 = ga-1.

    do i=1,nval
       rho = UT(1,i)
       vel = UT(2,i)
       pre = UT(3,i)

       mdot(i) = type_val*rho*vel*F_v(i)

       if(vel*type_val.gt.0.) then
          ! Outflow from end pipe

          rho = UP(1,i)
          vel = UP(2,i)
          pre = UP(3,i)
          
          hdot(i) = mdot(i)*(ga*pre/rho/g1+0.5*vel**2)
       else
          ! Inflow to end pipe
          hdot(i) = mdot(i)*ga*UC(2)/UC(1)/g1
       endif
    end do

  end subroutine valve_flow

  subroutine area_valve(theta, Dv, VO, VC, Lvmax, Nval, type_dat, &
       Lv_array, Cd, theta_cycle, F_v)
    !
    !  Computes the passage area through valves
    !
    !  area_valve is called by: solve_cylinder
    !  area_valve calls the following subroutines and functions: interpolant
    !
    implicit none

    integer, intent(in) :: type_dat,Nval
    real*8, intent(in) :: Dv,VO,VC,Lvmax,theta,theta_cycle
    real*8, dimension(:,:), intent(in) :: Lv_array
    real*8, dimension(:,:), intent(in) :: Cd
    real*8, intent(out) :: F_v

    real*8 :: dangle,angle,alpha,Lv,Cd_int
    logical :: is_open

    is_open = .false.
    if(modulo(theta-VO,theta_cycle) >= 0.0 .and. &
         modulo(theta-VO,theta_cycle) <= &
         modulo(VC-VO,theta_cycle)) &
         is_open = .true.
    
    if (is_open) then
       if(type_dat.eq.0) then
          dangle = VC-VO
          angle  = theta-VO
          if(angle.lt.0.)  angle  = theta_cycle+angle
          if(dangle.lt.0.) dangle = theta_cycle+dangle
          alpha = angle/dangle*pi
          Lv    = Lvmax*dsin(alpha)**2
       elseif(type_dat.eq.-1) then
          dangle = VC-VO
          angle  = theta-VO
          if(angle.lt.0.)  angle  = theta_cycle+angle
          if(dangle.lt.0.) dangle = theta_cycle+dangle
          alpha = angle/dangle
          if(alpha.le.0.5) then
             Lv = Lvmax*(1.-dexp(-50.*alpha))
          else
             Lv = Lvmax*(1.-dexp(-50.*(1.-alpha)))
          endif
       else
          call interpolant(Lv_array(:,1), Lv_array(:,2), &
               theta*180./pi, Lv)
       endif
       
       call interpolant(Cd(:,1), Cd(:,2), Lv, Cd_int)

       F_v = Cd_int * Nval * pi * Dv * Lv
    else
       F_v = 0.0d0
    endif

  end subroutine area_valve

  subroutine solve_valve_eqcases(Up, U_S, U_R, RHS, U_C, n_P, psi, ga, &
       caso, U_P, U_T, ok_flag)
    !
    !  Solves the non-linear equation arising from the valve model for
    !  the specified case
    !
    !  solve_valve_eqcases is called by: solve_valve_implicit
    !  solve_valve_eqcases calls the following subroutines and functions:
    !  none
    !
    implicit none

    integer, intent(in) :: caso
    real *8, intent(in) :: psi, ga, n_P
    real *8, intent(inout) :: Up
    real *8, dimension(3), intent(in) :: U_S, U_R, RHS, U_C
    real *8, dimension(3), intent(out) :: U_P, U_T
    logical, intent(out) :: ok_flag

    real*8 :: C0,C1,C2,a_S,a_R,a_C,p_C,p_S,v_S,v_R,RHS2nd, &
         delta,F,tol,a_T,rho_S,rho_C
    real*8 ::A,B,C,B1,B2,B3,C_new,Up_new, &
         v_P,v_T,p_P,p_T,rho_T,rho_P,dF,dU,nu,relax, &
         rho_R,tol_ini,Up1,Up_check,p_R,T_C
    integer:: max_ils,max_iter,iter,ils
    logical:: solve_case, do_newton, do_ls
    complex*16 :: Fnew

    tol     = 1.0d-10
    tol_ini = 1.0d-10

    delta = 0.5*(ga-1.)

    max_iter = 30
    max_ils  = 10

    ok_flag    = .false.
    solve_case = .false.
    do_newton  = .true.

    rho_C = U_C(1)
    p_C   = U_C(2)
    T_C   = U_C(3)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_R = U_R(1)
    v_R   = U_R(2)
    p_R   = U_R(3)

    a_S = dsqrt(ga*p_S/rho_S)
    a_R = dsqrt(ga*p_R/rho_R)
    a_C = dsqrt(ga*p_C/rho_C)

    ! esto es para obtener RHS
    if(n_P.gt.0) then
       RHS2nd = RHS(3)/a_C ! RHS_2^+
    else
       RHS2nd = RHS(2)/a_C ! RHS_2^-
    end if

    Up1 = 1./ga*a_S/a_C*(1.-p_C/p_S)+v_S/a_C*n_P-RHS2nd
    if(caso==1) then
       if(Up < Up1.or.Up > 0.) Up = Up1+tol_ini
       if(Up < 0.) solve_case = .true.
    elseif(caso==2) then
       if(Up > Up1.or.Up < 0.) Up = Up1-tol_ini
       if(Up > 0.) solve_case = .true.
       a_T = a_R*(p_C/p_R)**(delta/ga)*dexp(0.5*RHS(1))
    elseif(caso==3) then
       solve_case = .true.
    elseif(caso==4) then
       Up1 = 1./ga*a_S/a_C+v_S/a_C*n_P-RHS2nd
       if(Up > Up1.or.Up < 0.) Up = Up1-tol_ini
       if(Up > 0.) solve_case = .true.
    else
       stop ' Wrong case pased to subroutine solve_valve_eqcases '
    end if
    !
    if(solve_case) then
       iter = 0
       nu   = 2.
       !  estos son terminos que se repiten en C1 para todos los casos
       A = -1./ga*(a_S/a_C)*(p_C/p_S)
       B = n_P*v_S/a_C
       do while (do_newton)
          do_ls = .true.
          iter = iter+1
          C = Up-1./ga*(a_S/a_C)

          ! de a cuerdo al caso que se trata selecciona
          ! las ecuaciones correspondientes para resolver Newton
          if (caso==1) then
             !
             ! up*np<0,uT<aT
             !
             C1 = A/(C-B+RHS2nd)
             C2 = (C1**((ga-1.)/ga)-1.)/delta
             if(C2.lt.0.) exit ! do_newton = .false.
             F = Up/(1.-delta*Up**2)+psi*C1**(delta/ga)*dsqrt(C2)
             dF = (1.+delta*Up**2)/(1.-delta*Up**2)**2+psi*(a_C/a_S)*(p_S/p_C)*&
                  C1**((3.*ga-1.)/(2.*ga))/dsqrt(C2)*(2.*C1**((ga-1.)/ga)-1.)
          elseif(caso==2) then
             !
             ! up*np>0,uT<aT
             !
             C1 = (C-B+RHS2nd)/A
             C2 = (C1**((ga-1.)/ga)-1.)/delta
             if(C2.lt.0.) exit ! do_newton = .false.
             F = Up*dsqrt(C1**(2./ga)-psi**2)-psi*a_T/a_C*dsqrt(C2)
             dF = dsqrt(C1**(2./ga)-psi**2)-a_C/a_S*p_S/p_C*(Up*C1**((2.-ga)/ga)/ &
                  dsqrt(C1**(2./ga)-psi**2)-psi*dsqrt(delta)*a_T/a_C/C1**(1./ga)/ &
                  dsqrt(C1**((ga-1.)/ga)-1.))
          elseif(caso==3) then
             !
             ! up*np<0,uT=aT
             !
             B1 = -A*p_S/p_C+B-RHS2nd
             B2 = psi*(delta+1.)**(-(ga+1.)/(2.*(ga-1.)))*A
             B3 = B1**2-4.*B2*(1.-delta*B2)

             if(B3 >= 0.) then
                Up = (B1-dsqrt(B3))/(2.*(1.-delta*B2))
                if(Up < 0.) then
                   Up_check = (B1+dsqrt(B3))/(2.*(1.-delta*B2))
                   if(Up_check < 0.) then
                      ! Two possible solutions for this case
                      write(*,*) 'WARNING - valve case: 3 '
                   end if
                   ok_flag = .true.
                else
                   Up = (B1+dsqrt(B3))/(2.*(1.-delta*B2))
                   if(Up < 0.) then
                      ok_flag = .true.
                   end if
                end if
                do_newton = .false.
                do_ls     = .false.
             end if
          elseif(caso==4) then
             !
             ! up*np>0,uT=aT
             !
             C0 = a_C/a_R*dexp(-0.5*RHS(1))
             C1 = (A*p_R/p_C)/(C-B+RHS2nd)
             C2 = 1./(delta+1.)*(1.+delta*C0**2*Up**2*C1**((ga-1.)/ga))
             F  = C0*Up*C1**((ga-1.)/(2.*ga))-psi*C2**((ga+1.)/(2.*(ga-1.)))
             dF = C0*(C1**((ga-1.)/(2.*ga))+delta*a_C/a_S*p_S/p_R*Up* &
                  C1**((3.*ga-1.)/(2.*ga)))-psi*C0**2*C2**((3.-ga)/(2.*(ga-1.)))* &
                  (Up*C1**((ga-1.)/ga)+delta*a_C/a_S*p_S/p_R*Up**2*C1**((2.*ga-1.)/ga))
          end if
          !-------------------------------------------------------------------
          dU  = -F/dF
          ils = 0
          !
          ! Line search
          !
          relax = 1.
          do while(do_ls)
             ils = ils+1
             Up_new = Up+relax*dU
             C_new  = Up_new-1./ga*(a_S/a_C)
             !
             if(caso==1) then
                C1 = A/(C_new-B+RHS2nd)  
                C2 = (C1**((ga-1.)/ga)-1.)/delta
                Fnew = Up_new/(1.-delta*Up_new**2)+psi*C1**(delta/ga)*dsqrt(C2)
             elseif(caso==2) then
                C1 = (C_new-B+RHS2nd)/A
                Fnew = Up_new*dsqrt(C1**(2./ga)-psi**2)-psi/dsqrt(delta)* &
                     a_T/a_C*dsqrt(C1**((ga-1.)/ga)-1.)
             elseif(caso==4) then
                C1 = (A*p_R/p_C)/(C_new-B+RHS2nd)
                C2 = 1./(delta+1.)*(1.+delta*C0**2*Up_new**2*C1**((ga-1.)/ga))
                Fnew = C0*Up_new*C1**((ga-1.)/(2.*ga))-psi*C2**((ga+1.)/(2.*(ga-1.)))
             end if

             if((aimag(Fnew)==0 .and. dabs(dreal(Fnew)) < dabs(F)) .or. ils>=max_ils) then
                F  = dreal(Fnew)
                Up = Up_new
                do_ls = .false.
             else
                relax = relax/nu
             end if
          end do  ! Line-seach ends here
          
          if(dabs(F) < tol) then
             ok_flag   = .true.
             do_newton = .false.
          end if

          if(iter > max_iter) then
             do_newton = .false.
          end if

       end do ! Newton loop ends here
    end if !fin de entrar

    !
    ! Variables computation
    !
    if(ok_flag) then
       C   = Up-1./ga*(a_S/a_C)
       v_P = a_C*Up*n_P
       p_P = -p_S*ga*(a_C/a_S)*(C-B+RHS2nd)
       if(caso==1) then
          p_T   = p_P
          rho_P = ga*p_P/(a_C**2*(1.-delta*Up**2))
          rho_T = rho_C*(p_T/p_C)**(1./ga)
          v_T   = rho_P*v_P/(rho_T*psi)
       elseif(caso==2) then
          p_T   = p_C
          rho_T = ga*p_T/a_T**2
          rho_P = rho_T*(p_P/p_T)**(1./ga)
          v_T   = rho_P*v_P/(rho_T*psi)
       elseif(caso==3) then
          p_T   = p_C*(delta+1.)**(-ga/(ga-1.))
          rho_T = rho_C*(p_T/p_C)**(1./ga)
          v_T   = -n_P*a_C*(p_T/p_C)**((ga-1.)/(2.*ga))
          rho_P = rho_T*v_T*psi/v_P
       elseif(caso==4) then
          rho_P = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS(1))
          v_T   = dsqrt(1./(delta+1.)*(ga*p_P/rho_P+delta*v_P**2))
          rho_T = rho_P*v_P/(v_T*psi)
          p_T   = p_P*(rho_T/rho_P)**ga
       end if
       U_P(1) = rho_P
       U_P(2) = v_P
       U_P(3) = p_P

       U_T(1) = rho_T
       U_T(2) = v_T
       U_T(3) = p_T
    end if

  end subroutine solve_valve_eqcases

  subroutine get_valve_case(U_P, U_T, U_C, n_P, ga, caso)
    !
    !  Given the states at the pipe end (U_P) and the one at the 
    !  nozzle throat (U_T), this subroutine computes the corresponding 
    !  valve case
    !
    !  get_valve_case is called by: solve_valve_implicit
    !  get_valve_case calls the following subroutines and functions:
    !  none
    !
    implicit none

    real*8, dimension(3), intent(in) :: U_P,U_T,U_C
    real*8, intent(in) :: n_P,ga
    integer, intent(out) :: caso

    real*8 :: delta,rho_P,v_P,p_P,rho_T,v_T,p_T,a_P,a_T,tol,p_C
    logical :: subsonic_T

    tol   = 1d-10
    delta = 0.5*(ga-1.)

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_T = U_T(1)
    v_T   = U_T(2)
    p_T   = U_T(3)

    p_C = U_C(2)

    a_P = dsqrt(ga*p_P/rho_P)
    a_T = dsqrt(ga*p_T/rho_T)

    ! n_P is the external unit normal to the pipe
    if(-n_P*v_T > 0.) then
       subsonic_T = p_C/p_P < (delta+1.)**(ga/(ga-1.))
    else
       subsonic_T = p_P/p_C < ((delta+1.)/(1.+delta*(v_P/a_P)**2))**(ga/(ga-1.))
    end if
    ! subsonic_T = dabs(v_T) < a_T-tol
    if(n_P*v_P > 0.) then
       if(subsonic_T) then ! subsonic flow at the throat
          caso = 2
       else ! sonic flow at the throat
          caso = 4
       end if
    else
       if(subsonic_T) then ! subsonic flow at the throat
          caso = 1
       else ! sonic flow at the throat
          caso = 3
       end if
    end if

  end subroutine get_valve_case

  subroutine solve_valve_eqcases2(Up, U_S, U_R, RHS, U_C, n_P, psi, ga, &
       caso, U_P, U_T, ok_flag)
    !
    !  Solves the non-linear equation arising from the valve model for
    !  the specified case
    !
    !  solve_valve_eqcases is called by: solve_valve_implicit
    !  solve_valve_eqcases calls the following subroutines and functions:
    !  none
    !
    implicit none

    integer, intent(in) :: caso
    real *8, intent(in) :: psi, ga, n_P
    real *8, intent(inout) :: Up
    real *8, dimension(3), intent(in) :: U_S, U_R, RHS, U_C
    real *8, dimension(3), intent(out) :: U_P, U_T
    logical, intent(out) :: ok_flag

    real*8 :: C0,C1,C2,a_S,a_R,a_C,p_C,p_S,v_S,v_R,RHS2nd, &
         delta,F,tol,a_T,rho_S,rho_C
    real*8 ::A,B,C,B1,B2,B3,C_new,Up_new, &
         v_P,v_T,p_P,p_T,rho_T,rho_P,dF,dU,nu,relax, &
         rho_R,tol_ini,Up1,Up_check,p_R,T_C,a_P
    integer:: max_ils,max_iter,iter,ils
    logical:: solve_case, do_newton, do_ls
    complex*16 :: Fnew

    tol     = 1.0d-10
    tol_ini = 1.0d-10

    delta = 0.5*(ga-1.)

    max_iter = 30
    max_ils  = 10

    ok_flag    = .false.
    solve_case = .false.
    do_newton  = .true.

    rho_C = U_C(1)
    p_C   = U_C(2)
    T_C   = U_C(3)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_R = U_R(1)
    v_R   = U_R(2)
    p_R   = U_R(3)

    a_S = dsqrt(ga*p_S/rho_S)
    a_R = dsqrt(ga*p_R/rho_R)
    a_C = dsqrt(ga*p_C/rho_C)

    ! esto es para obtener RHS
    if(n_P.gt.0) then
       RHS2nd = RHS(3)/a_C ! RHS_2^+
    else
       RHS2nd = RHS(2)/a_C ! RHS_2^-
    end if

    Up1 = 1./ga*a_S/a_C*(1.-p_C/p_S)+v_S/a_C*n_P-RHS2nd
    if(caso==1) then
       if(Up < Up1.or.Up > 0.) Up = Up1+tol_ini
       if(Up < 0.) solve_case = .true.
    elseif(caso==2) then
       if(Up > Up1.or.Up < 0.) Up = Up1-tol_ini
       if(Up > 0.) solve_case = .true.
       a_T = a_R*(p_C/p_R)**(delta/ga)*dexp(0.5*RHS(1))
    elseif(caso==3) then
       solve_case = .true.
    elseif(caso==4) then
       Up1 = 1./ga*a_S/a_C+v_S/a_C*n_P-RHS2nd
       if(Up > Up1.or.Up < 0.) Up = Up1-tol_ini
       if(Up > 0.) solve_case = .true.
    else
       stop ' Wrong case pased to subroutine solve_valve_eqcases '
    end if
    !
    if(solve_case) then
       iter = 0
       nu   = 2.
       !  estos son terminos que se repiten en C1 para todos los casos
       A = -1./ga*(a_S/a_C)*(p_C/p_S)
       B = n_P*v_S/a_C
       do while (do_newton)
          do_ls = .true.
          iter = iter+1
          C = Up-1./ga*(a_S/a_C)

          ! de a cuerdo al caso que se trata selecciona
          ! las ecuaciones correspondientes para resolver Newton
          if (caso==1) then
             !
             ! up*np<0,uT<aT
             !
             C1 = A/(C-B+RHS2nd)
             C2 = (C1**((ga-1.)/ga)*(1.-delta*Up**2)-1.)/delta
             if(C2.lt.0.) exit ! do_newton = .false.
             F = Up+ &
                  psi*C1**(delta/ga)*(1.-delta*Up**2.)**((ga-3.)/2./(ga-1.))* &
                  dsqrt(C2)
             dF = 1.+psi*delta*p_S/p_C*a_C/a_S*C1**((3.*ga-1.)/2./ga)* &
                  (1.-delta*Up**2.)**((ga-3.)/2./(ga-1.))*dsqrt(C2)- &
                  psi*(ga-3.)/2.*C1**(delta/ga)* &
                  Up/(1.-delta*Up**2)**((ga+1.)/2./(ga-1.))*dsqrt(C2)- &
                  psi*C1**(delta/ga)* &
                  (1.-delta*Up**2.)**((ga-3.)/2./(ga-1.))/dsqrt(C2)* &
                  (-p_S/p_C*a_C/a_S*(1.-delta*Up**2.)*C1**((2.*ga-1)/ga)+ &
                  Up*C1**((ga-1.)/ga))
          elseif(caso==2) then
             !
             ! up*np>0,uT<aT
             !
             C1 = (C-B+RHS2nd)/A
             C2 = (C1**((ga-1.)/ga)-1.)/delta
             if(C2.lt.0.) exit ! do_newton = .false.
             F = Up*dsqrt(C1**(2./ga)-psi**2.)-psi*a_T/a_C*dsqrt(C2)
             dF = dsqrt(C1**(2./ga)-psi**2.)- &
                  a_C/a_S*p_S/p_C*(Up*C1**((2.-ga)/ga)/ &
                  dsqrt(C1**(2./ga)-psi**2.)- &
                  psi*dsqrt(delta)*a_T/a_C/C1**(1./ga)/ &
                  dsqrt(C1**((ga-1.)/ga)-1.))
          elseif(caso==3) then
             !
             ! up*np<0,uT=aT
             !
             B1 = -A*p_S/p_C+B-RHS2nd
             B2 = psi*(delta+1.)**(-(ga+1.)/(2.*(ga-1.)))*A
             B3 = B1**2-4.*B2*(1.-delta*B2)

             if(B3 >= 0.) then
                Up = (B1-dsqrt(B3))/(2.*(1.-delta*B2))
                if(Up < 0.) then
                   Up_check = (B1+dsqrt(B3))/(2.*(1.-delta*B2))
                   if(Up_check < 0.) then
                      ! Two possible solutions for this case
                      write(*,*) 'WARNING - valve case: 3 '
                   end if
                   ok_flag = .true.
                else
                   Up = (B1+dsqrt(B3))/(2.*(1.-delta*B2))
                   if(Up < 0.) then
                      ok_flag = .true.
                   end if
                end if
                do_newton = .false.
                do_ls     = .false.
             end if
          elseif(caso==4) then
             !
             ! up*np>0,uT=aT
             !
             C0 = a_C/a_R*dexp(-0.5*RHS(1))
             C1 = (A*p_R/p_C)/(C-B+RHS2nd)
             C2 = 1./(delta+1.)*(1.+delta*C0**2*Up**2*C1**((ga-1.)/ga))
             F  = C0*Up*C1**((ga-1.)/(2.*ga))-psi*C2**((ga+1.)/(2.*(ga-1.)))
             dF = C0*(C1**((ga-1.)/(2.*ga))+delta*a_C/a_S*p_S/p_R*Up* &
                  C1**((3.*ga-1.)/(2.*ga)))- &
                  psi*C0**2*C2**((3.-ga)/(2.*(ga-1.)))* &
                  (Up*C1**((ga-1.)/ga)+ &
                  delta*a_C/a_S*p_S/p_R*Up**2*C1**((2.*ga-1.)/ga))
          end if
          !--------------------------------------------------------------------
          dU  = -F/dF
          ils = 0
          !
          ! Line search
          !
          relax = 1.
          do while(do_ls)
             ils = ils+1
             Up_new = Up+relax*dU
             C_new  = Up_new-1./ga*(a_S/a_C)
             !
             if(caso==1) then
                C1 = A/(C_new-B+RHS2nd)
                C2 = (C1**((ga-1.)/ga)*(1.-delta*Up_new**2)-1.)/delta
                Fnew = Up_new+psi*C1**(delta/ga)* &
                     (1.-delta*Up_new**2.)**((ga-3.)/2./(ga-1.))*dsqrt(C2)
             elseif(caso==2) then
                C1 = (C_new-B+RHS2nd)/A
                Fnew = Up_new*dsqrt(C1**(2./ga)-psi**2.)-psi/dsqrt(delta)* &
                     a_T/a_C*dsqrt(C1**((ga-1.)/ga)-1.)
             elseif(caso==4) then
                C1 = (A*p_R/p_C)/(C_new-B+RHS2nd)
                C2 = 1./(delta+1.)*(1.+delta*C0**2.*Up_new**2*C1**((ga-1.)/ga))
                Fnew = C0*Up_new*C1**((ga-1.)/(2.*ga))- &
                     psi*C2**((ga+1.)/(2.*(ga-1.)))
             end if

             if((aimag(Fnew)==0 .and. dabs(dreal(Fnew)) < dabs(F)) .or. ils>=max_ils) then
                F  = dreal(Fnew)
                Up = Up_new
                do_ls = .false.
             else
                relax = relax/nu
             end if
          end do  ! Line-seach ends here
          
          if(dabs(F) < tol) then
             ok_flag   = .true.
             do_newton = .false.
          end if

          if(iter > max_iter) then
             do_newton = .false.
          end if

       end do ! Newton loop ends here
    end if !fin de entrar

    !
    ! Variables computation
    !
    if(ok_flag) then
       C   = Up-1./ga*(a_S/a_C)
       v_P = a_C*Up*n_P
       p_P = -p_S*ga*(a_C/a_S)*(C-B+RHS2nd)
       if(caso==1) then
          rho_P = ga*p_P/(a_C**2.*(1.-delta*Up**2.))
          a_P   = dsqrt(ga*p_P/rho_P)
          p_T   = p_P*(1.+delta*(v_P/a_P)**2.)**(ga/(ga-1.))
          rho_T = rho_C*(p_T/p_C)**(1./ga)
          v_T   = rho_P*v_P/(rho_T*psi)
       elseif(caso==2) then
          p_T   = p_C
          rho_T = ga*p_T/a_T**2
          rho_P = rho_T*(p_P/p_T)**(1./ga)
          v_T   = rho_P*v_P/(rho_T*psi)
       elseif(caso==3) then
          p_T   = p_C*(delta+1.)**(-ga/(ga-1.))
          rho_T = rho_C*(p_T/p_C)**(1./ga)
          v_T   = -n_P*a_C*(p_T/p_C)**((ga-1.)/(2.*ga))
          rho_P = rho_T*v_T*psi/v_P
       elseif(caso==4) then
          rho_P = rho_R*(p_P/p_R)**(1./ga)*dexp(-RHS(1))
          v_T   = dsqrt(1./(delta+1.)*(ga*p_P/rho_P+delta*v_P**2))
          rho_T = rho_P*v_P/(v_T*psi)
          p_T   = p_P*(rho_T/rho_P)**ga
       end if
       U_P(1) = rho_P
       U_P(2) = v_P
       U_P(3) = p_P

       U_T(1) = rho_T
       U_T(2) = v_T
       U_T(3) = p_T
    end if

  end subroutine solve_valve_eqcases2

end module def_valve
