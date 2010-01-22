module def_valve

  use def_simulator
  use utilities
  use gasdyn_utils, only : compute_cp

  type :: valve
     integer :: Nval, type_dat, tube, valve_model
     real*8, dimension(:,:), pointer  :: Cd, Lv
     real*8 :: angle_VO, angle_VC, Dv, Lvmax
     real*8 :: Area_tube, twall_tube, dAreax_tube, dx_tube
  end type valve

contains

  subroutine solve_valve(globalData, Ucyl, Uref, type_val, Area_val, Area_pipe, &
       Upipe, Uthroat)
    !
    !
    !  solve_valve is called by: solve_cylinder
    !  solve_valve calls the following subroutines and functions: 
    !    solve_valve_equation
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: Area_val,Area_pipe
    real*8, dimension(3), intent(in) :: Ucyl,Uref
    type(DataSim), intent(in) :: globalData
    real*8, dimension(3), intent(out) :: Upipe,Uthroat

    real*8 :: tol,ga,R_gas,cp,g1,delta,psi,rp,rp_cr
    real*8 :: rho_P,v_P,p_P,T_P,a_P,rho_T,v_T,p_T,a_T
    real*8 :: rho0,v0,p0,a0,T0,a_ref,p_C
    real*8 :: A,A_0,At,U,Ut,rA,v1
    real*8, dimension(3) :: param
    logical :: noflow,inflow

    tol = 1.0d-10

    ga    = globalData%ga
    R_gas = globalData%R_gas

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
       RHS, ga, Upipe, Uthroat)
    !
    !
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: ga,Area_val,Area_pipe
    real*8, dimension(3), intent(in) :: RHS,Ucyl,Uref
    real*8, dimension(3), intent(inout) :: Upipe
    real*8, dimension(3), intent(out) :: Uthroat

    integer :: niter,niter_max
    logical :: continue_flag
    real*8 :: normres,tol,relax
    real*8, dimension(6) :: R,dU
    real*8, dimension(6,6) :: J
    real*8, dimension(6,9) :: mat

    niter_max = 50
    tol       = 1.0d-6
    relax     = 1.

    ! Upipe = Uref
    call state_valve_predictor(Upipe, Ucyl, Uref, RHS, type_val, Area_val, Area_pipe, &
         ga, Uthroat)

    niter = 0
    continue_flag = .true.
    do while (continue_flag)
       R = 0.0d0
       J = 0.0d0
       call valve_residue(Uthroat, Upipe, Ucyl, Uref, ga, Area_val, Area_pipe, type_val, &
            RHS, R)
       call valve_analytical_jaco(Uthroat, Upipe, Ucyl, Uref, ga, Area_val, Area_pipe, type_val, &
            RHS, mat)
       J = mat(:,1:6)
       call linear_solver(R, J, dU, 6)
       relax = tanh(0.05*(niter+1))
       Upipe   = Upipe-relax*dU(1:3)
       Uthroat = Uthroat-relax*dU(4:6)

       normres = normvec(R, 6, 2)
       niter   = niter+1
       if((normres.lt.tol).or.(niter.gt.niter_max)) continue_flag = .false.
    end do

  end subroutine solve_valve_implicit

  subroutine valve_residue(U_T,U_P,U_C,Uref,ga,Area_T,Area_P,type_end, &
       RHS,res)
    !
    !  Computes the residue on valve model
    !
    !  valve_residue is called by:
    !  valve_residue calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: type_end
    real*8, intent(in) :: ga,Area_T,Area_P
    real*8, dimension(3), intent(in) :: U_T,U_P,U_C,Uref,RHS
    real*8, dimension(6), intent(out) :: res

    real*8 :: delta,g1,p_0P,psi
    real*8 :: Area_min,sigma,alpha
    real*8 :: rho_T,v_T,p_T,a_T
    real*8 :: rho_P,v_P,p_P,a_P
    real*8 :: rho_C,T_C,p_C,a_C
    real*8 :: rho_S,v_S,p_S,a_S
    real*8 :: r_sonic,r_subsonic

    Area_min = 1d-10

    g1    = ga-1.
    delta = .5*g1

    rho_T = U_T(1)
    v_T   = U_T(2)
    p_T   = U_T(3)

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_C = U_C(1)
    p_C   = U_C(2)
    T_C   = U_C(3)

    rho_S = Uref(1)
    v_S   = Uref(2)
    p_S   = Uref(3)

    res = 0d0
    
    if (Area_T.gt.Area_min) then
       a_C = dsqrt(ga*p_C/rho_C)
       a_T = dsqrt(ga*p_T/rho_T)
       a_P = dsqrt(ga*p_P/rho_P)
       a_S = dsqrt(ga*p_S/rho_S)

       p_0P = p_P*(1.+delta*(v_P/a_P)**2)**(ga/g1)

       ! Compatibility along Mach lines (lambda minus)
       if (type_end.gt.0) then
          res(1) = (p_P-p_S)/(rho_S*a_S)-v_S+RHS(3)+v_P
       else
          res(1) = (p_P-p_S)/(rho_S*a_S)+v_S+RHS(2)-v_P
       end if
       ! Continuity throat-pipe
       psi    = Area_T/Area_P
       res(2) = rho_T*v_T*psi-rho_P*v_P

       if (p_C.gt.p_0P) then
          ! From cylinder to pipe through valve (IMBOCCO)

          ! Energy conservation cylinder to pipe
          res(3) = a_C**2-(a_P**2+delta*v_P**2)
          ! Isentropic cylinder to throat
          res(4) = p_C/p_T-(rho_C/rho_T)**ga
          ! Energy conservation cylinder throat
          res(5) = a_C**2-(a_T**2+delta*v_T**2)
          ! Isobaric throat-cylinder
          sigma      = .025*a_T
          alpha      = .5*(1.+dtanh((a_T-v_T)/sigma))
          r_sonic    = a_T-v_T
          r_subsonic = p_T-p_P
          res(6)     = (1.-alpha)*r_sonic+alpha*r_subsonic
       else
          ! From pipe to cylinder through valve (SBOCCO)

          ! Compatibility along path line (u)
          res(3) = (p_P/p_S)**(1./ga)-(rho_P/rho_S)*dexp(RHS(1))
          ! Isoentropic pipe-throat
          res(4) = p_P/p_T-(rho_P/rho_T)**ga
          ! Energy conservation pipe throat
          res(5) = a_P**2+delta*v_P**2-(a_T**2+delta*v_T**2)
          ! Isobaric throat-cylinder
          sigma      = .025*a_T
          alpha      = .5*(1.+dtanh((a_T-v_T)/sigma))
          r_sonic    = a_T-v_T
          r_subsonic = p_T-p_C
          res(6)     = (1.-alpha)*r_sonic+alpha*r_subsonic
       end if
    else
       ! Compatibility along Mach lines (lambda minus)
       a_S = dsqrt(ga*p_S/rho_S)
       if (type_end.gt.0) then
          res(1) = (p_P-p_S)/(rho_S*a_S)-v_S+RHS(2)+v_P
       else
          res(1) = (p_P-p_S)/(rho_S*a_S)+v_S+RHS(3)-v_P
       end if
       res(2) = v_P
       res(3) = (p_P/p_S)**(1./ga)-(rho_P/rho_S)*dexp(RHS(1))
       res(4) = v_T
       res(5) = rho_T-rho_P
       res(6) = p_T-p_P
    end if

  end subroutine valve_residue

  subroutine valve_analytical_jaco(U_T,U_P,U_C,Uref,ga,Area_T,Area_P,type_end, &
       RHS,mat)
    !
    !  Computes the analytical jacobian on valve model
    !
    !  valve_analytical_jaco is called by:
    !  valve_analytical_jaco calls the following subroutines and 
    !    functions: none
    !
    implicit none

    integer, intent(in) :: type_end
    real*8, intent(in) :: ga,Area_T,Area_P
    real*8, dimension(3), intent(in) :: U_T,U_P,U_C,Uref,RHS
    real*8, dimension(6,9), intent(out) :: mat

    real*8 :: delta,g1,p_0P,psi
    real*8 :: Area_min,sigma,alpha,dalphadpG,dalphadrhoG,dalphaduG
    real*8 :: rho_T,v_T,p_T,a_T
    real*8 :: rho_P,v_P,p_P,a_P
    real*8 :: rho_C,T_C,p_C,a_C
    real*8 :: rho_S,v_S,p_S,a_S

    Area_min = 1d-10

    g1    = ga-1.
    delta = .5*g1

    rho_T = U_T(1)
    v_T   = U_T(2)
    p_T   = U_T(3)

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_C = U_C(1)
    p_C   = U_C(2)
    T_C   = U_C(3)

    rho_S = Uref(1)
    v_S   = Uref(2)
    p_S   = Uref(3)

    mat = 0d0

    if(Area_T.gt.Area_min) then
       a_C = dsqrt(ga*p_C/rho_C)
       a_T = dsqrt(ga*p_T/rho_T)
       a_P = dsqrt(ga*p_P/rho_P)
       a_S = dsqrt(ga*p_S/rho_S)
			
       p_0P = p_P*(1.+delta*(v_P/a_P)**2)**(ga/g1)

       ! Compatibility along Mach lines (lambda minus)
       if(type_end.gt.0) then
          mat(1,2) = 1.
       else
          mat(1,2) = -1.
       endif
       mat(1,3) = 1./(rho_S*a_S)
       ! Continuity throat-pipe
       psi      = Area_T/Area_P
       mat(2,1) = -v_P
       mat(2,2) = -rho_P
       mat(2,4) = v_T*psi
       mat(2,5) = rho_T*psi
       
       if(p_C.gt.p_0P) then
          ! From cylinder to pipe through valve (IMBOCCO)
          
          ! Energy conservation cylinder to pipe
          mat(3,1) = ga*p_P/rho_P**2
          mat(3,2) = -2.*delta*v_P
          mat(3,3) = -ga/rho_P
          mat(3,7) = -ga*p_C/rho_C**2
          mat(3,8) = ga/rho_C
          ! Isoentropic cylinder to throat
          mat(4,4) = ga/rho_T*(rho_C/rho_T)**ga
          mat(4,6) = -p_C/p_T**2
          mat(4,7) = -ga/rho_C*(rho_C/rho_T)**ga
          mat(4,8) = 1./p_T
          ! Energy conservation cylinder throat
          mat(5,4) = ga*p_T/rho_T**2
          mat(5,5) = -2.*delta*v_T
          mat(5,6) = -ga/rho_T
          mat(5,7) = -ga*p_C/rho_C**2
          mat(5,8) = ga/rho_C
          ! Isobaric throat-cylinder
          sigma = .025*a_T
          alpha = .5*(1.+dtanh((a_T-v_T)/sigma))
          dalphadrhoG = -v_T/(4.*sigma*rho_T)*(1.-dtanh((a_T-v_T)/sigma)**2)
          dalphaduG   = -1./(2.*sigma)*(1.-dtanh((a_T-v_T)/sigma)**2)
          dalphadpG   = v_T/(4.*sigma*p_T)*(1.-dtanh((a_T-v_T)/sigma)**2)
          mat(6,3) = -alpha
          mat(6,4) = -a_T/(2.*rho_T)*(1.-alpha)+dalphadrhoG*((p_T-p_P)-(a_T-v_T))
          mat(6,5) = -1.+alpha+dalphaduG*((p_T-p_P)-(a_T-v_T))
          mat(6,6) = alpha+a_T/(2.*p_T)*(1.-alpha)+dalphadpG*((p_T-p_P)-(a_T-v_T))
       else
          ! From pipe to cylinder through valve (SBOCCO)
          
          ! Compatibility along path line (u)
          mat(3,1) = -dexp(RHS(1))/rho_S
          mat(3,3) = 1./(ga*p_P)*(p_P/p_S)**(1./ga)
          ! Isoentropic pipe-throat
          mat(4,1) = -ga/rho_P*(rho_P/rho_T)**ga
          mat(4,3) = 1./p_T
          mat(4,4) = ga/rho_T*(rho_P/rho_T)**ga
          mat(4,5) = -p_P/p_T**2
          ! Energy conservation pipe throat
          mat(5,1) = -ga*p_P/rho_P**2
          mat(5,2) = 2.*delta*v_P
          mat(5,3) = ga/rho_P
          mat(5,4) = ga*p_T/rho_T**2
          mat(5,5) = -2.*delta*v_T
          mat(5,6) = -ga/rho_T
          ! Isobaric throat-cylinder
          sigma = .025*a_T
          alpha = .5*(1.+dtanh((a_T-v_T)/sigma))
          dalphadrhoG = -v_T/(4.*sigma*rho_T)*(1.-dtanh((a_T-v_T)/sigma)**2)
          dalphaduG   = -1./(2.*sigma)*(1.-dtanh((a_T-v_T)/sigma)**2)
          dalphadpG   = v_T/(4.*sigma*p_T)*(1.-dtanh((a_T-v_T)/sigma)**2)
          mat(6,4) = -a_T/(2.*rho_T)*(1.-alpha)+dalphadrhoG*((p_T-p_C)-(a_T-v_T))
          mat(6,5) = -1.+alpha+dalphaduG*((p_T-p_C)-(a_T-v_T))
          mat(6,6) = alpha+a_T/(2.*p_T)*(1.-alpha)+dalphadpG*((p_T-p_C)-(a_T-v_T))
          mat(6,8) = -alpha
       endif
    else
       ! Compatibility along Mach lines (lambda minus)
       a_S  = dsqrt(ga*p_S/rho_S)
       if(type_end.gt.0) then
          mat(1,2) = 1.
       else
          mat(1,2) = -1.
       endif
       mat(1,3) = 1./(rho_S*a_S)
       mat(2,2) = 1.
       mat(3,1) = -dexp(RHS(1))/rho_S
       mat(3,3) = 1./(ga*p_P)*(p_P/p_S)**(1./ga)
       mat(4,5) = 1.
       mat(5,4) = 1.
       mat(5,1) = -1.
       mat(6,6) = 1.
       mat(6,3) = -1.
    endif

  end subroutine valve_analytical_jaco

  subroutine state_valve_predictor(U_P, U_C, U_S, RHS, type_val, Area_T, Area_P, ga, U_T)
    !
    !
    !
    implicit none

    integer, intent(in) :: type_val
    real*8, intent(in) :: ga,Area_T,Area_P
    real*8, dimension(3), intent(in) :: U_S,U_C,RHS
    real*8, dimension(3), intent(inout) :: U_P
    real*8, dimension(3), intent(out) :: U_T

    real*8 :: delta,g1,psi,Area_min,beta,beta_cr
    real*8 :: rho_0P,p_0P,a_0P
    real*8 :: rho_P,v_P,p_P,a_P,rho_C,p_C,a_C
    real*8 :: rho_T,v_T,p_T,rho_S,v_S,p_S,a_S

    Area_min = 1d-10

    g1    = ga-1.
    delta = .5*g1

    beta_cr = (2./(ga+1.))**(ga/g1)

    rho_P = U_P(1)
    v_P   = U_P(2)
    p_P   = U_P(3)

    rho_S = U_S(1)
    v_S   = U_S(2)
    p_S   = U_S(3)

    rho_C = U_C(1)
    p_C   = U_C(2)

    if(Area_T.gt.Area_min) then
       psi = Area_T/Area_P
       a_C = dsqrt(ga*p_C/rho_C)
       a_S = dsqrt(ga*p_S/rho_S)

       if(type_val.gt.0) then
          v_P = -(p_P-p_S)/(rho_S*a_S)+v_S-RHS(3)
       else
          v_P = (p_P-p_S)/(rho_S*a_S)+v_S+RHS(2)
       end if

       if(v_P*type_val.gt.0.) then
          a_P   = dsqrt(a_C**2-delta*v_P**2)
          rho_P = ga*p_P/a_P**2
       else
          rho_P = rho_S*(p_P/p_S)**(1./ga)*dexp(-RHS(1))
          a_P   = dsqrt(ga*p_P/rho_P)
       end if

       p_0P   = p_P*(1.+delta*(v_P/a_P)**2)**(ga/g1)
       rho_0P = rho_P*(p_0P/p_P)**(1./ga)
       a_0P   = dsqrt(ga*p_0P/rho_0P)
       beta   = p_0P/p_C

       if(beta.gt.1./beta_cr) then
          rho_T = rho_0P/(1.+delta)**(1./g1)
          p_T   = rho_T/rho_0P*p_0P/(1.+delta)
          v_T   = rho_0P/rho_T*a_0P*beta_cr**((ga+1.)/(2.*ga))
       else if(beta.gt.1.) then
          rho_T = rho_P*(p_C/p_P)**(1./ga)
          p_T   = p_C
          v_T   = rho_0P/rho_T*a_0P*dsqrt(1./delta*((p_C/p_0P)**(2./ga)- &
               (p_C/p_0P)**((ga+1.)/ga)))
       else if(beta.lt.beta_cr) then
          rho_T = rho_C/(1.+delta)**(1./g1)
          p_T   = rho_T/rho_C*p_C/(1.+delta)
          v_T   = -rho_C/rho_T*a_C*beta_cr**((ga+1.)/(2.*ga))
       else
          rho_T = rho_C*(p_P/p_C)**(1./ga)
          p_T   = p_P
          v_T   = -rho_C/rho_T*a_C*dsqrt(1./delta*((p_P/p_C)**(2./ga)- &
               (p_P/p_C)**((ga+1.)/ga)))
       end if
    else
       rho_T = rho_P
       v_T   = 0.0d0
       p_T   = p_P
    end if

    U_T(1) = rho_T
    U_T(2) = type_val*v_T
    U_T(3) = p_T

    U_P(1) = rho_P
    U_P(2) = type_val*sign(v_P, v_T)
    U_P(3) = p_P

  end subroutine state_valve_predictor

  subroutine valve_flow(nval, type_val, ga, F_v, UT, UP, UC, &
       nmdot, nhdot)
    !
    !
    !  valve_flow is called by: cylinder_solver
    !  valve_flow calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: nval,type_val
    real*8, intent(in) :: ga
    real*8, dimension(3), intent(in) :: UC
    real*8, dimension(nval), intent(in) :: F_v
    real*8, dimension(3,nval), intent(in) :: UT,UP
    real*8, intent(out) :: nmdot,nhdot

    integer :: i
    real*8 :: g1,rho,vel,pre
    real*8, dimension(nval) :: mdot,hdot

    g1 = ga-1.

    do i=1,nval
       rho = UT(1,i)
       vel = UT(2,i)
       pre = UT(3,i)

       mdot(i) = type_val*rho*vel*F_v(i)

       rho = UP(1,i)
       vel = UP(2,i)
       pre = UP(3,i)

       if(vel*type_val.gt.0.) then
          ! Outflow from end pipe
          hdot(i) = mdot(i)*(ga*pre/rho/g1+0.5*vel**2)
       else
          ! Inflow to end pipe
          hdot(i) = mdot(i)*ga*UC(2)/UC(1)/g1
       endif
    end do

    nmdot = sum(mdot)
    nhdot = sum(hdot)

  end subroutine valve_flow

  subroutine area_valve(theta, Dv, VO, VC, Lvmax, Nval, type_dat, &
       Lv_array, Cd, nstroke, F_v)
    !
    !  Computes the passage area through valves
    !
    !  area_valve is called by: solve_cylinder
    !  area_valve calls the following subroutines and functions: interpolant
    !
    implicit none

    integer, intent(in) :: type_dat,Nval,nstroke
    real*8, intent(in) :: Dv,VO,VC,Lvmax,theta
    real*8, dimension(:,:), intent(in) :: Lv_array
    real*8, dimension(:,:), intent(in) :: Cd
    real*8, intent(out) :: F_v

    real*8 :: dangle,angle,alpha,Lv,Cd_int
    logical :: is_open

    is_open = .false.
    if(modulo(theta-VO,nstroke*pi) >= 0.0 .and. &
         modulo(theta-VO,nstroke*pi) <= &
         modulo(VC-VO,nstroke*pi)) &
         is_open = .true.
    
    if (is_open) then
       if(type_dat.eq.0) then
          dangle = VC-VO
          angle  = theta-VO
          if(angle.lt.0.)  angle  = nstroke*pi+angle
          if(dangle.lt.0.) dangle = nstroke*pi+dangle
          alpha = angle/dangle*pi
          Lv    = Lvmax*dsin(alpha)**2
       elseif(type_dat.eq.-1) then
          dangle = VC-VO
          angle  = theta-VO
          if(angle.lt.0.)  angle  = nstroke*pi+angle
          if(dangle.lt.0.) dangle = nstroke*pi+dangle
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

end module def_valve
