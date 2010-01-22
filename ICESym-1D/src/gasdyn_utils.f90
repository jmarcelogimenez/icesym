module gasdyn_utils

  use utilities, only : pi

contains

  subroutine rhschar(U, Area, dAreax, Twall, ga, R_gas, dt, &
       viscous_flow, heat_flow, RHS)
    !
    !  Computes the RHS of compatibility equations
    !
    !
    implicit none

    integer, intent(in) :: viscous_flow,heat_flow
    real*8, intent(in) :: Area,dAreax,Twall,dt,ga,R_gas
    real*8, dimension(3), intent(in) :: U
    real*8, dimension(3), intent(out) :: RHS

    real*8 :: g1,rho,vel,pre,T,a
    real*8 :: RHS1,RHS2m,RHS2p,alpha,G,qptot

    g1 = ga-1.

    rho = U(1)
    vel = U(2)
    pre = U(3)

    T = pre/(R_gas*rho)
    a = dsqrt(ga*pre/rho)

    call source_char(U, dAreax, area, Twall, ga, R_gas, viscous_flow, &
         heat_flow,alpha, G, qptot)

    ! RHS Path line
    RHS1 = g1/(ga*pre)*(qptot+rho*vel*G)*dt

    ! RHS Mach line minus
    RHS2m = (g1*(qptot+rho*vel*G)+rho*a*G-rho*vel*a**2*alpha)*dt
    RHS2m = -RHS2m/(rho*a)

    ! RHS Mach line plus
    RHS2p = (g1*(qptot+rho*vel*G)-rho*a*G-rho*vel*a**2*alpha)*dt
    RHS2p = -RHS2p/(rho*a)

    RHS(1) = RHS1
    RHS(2) = RHS2m
    RHS(3) = RHS2p

  end subroutine rhschar
  
  subroutine source_char(U, dAreax, area, Twall, ga, R_gas, &
       viscous_flow, heat_flow, alpha, G, qptot)
    !
    !  Computes the source term for 
    !    variable section,
    !    friction flow, and
    !    heat transfer
    !
    !  source is called by: rhschar
    !  source calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: viscous_flow,heat_flow
    real*8, intent(in) :: dAreax,Area,Twall,ga,R_gas
    real*8, dimension(3), intent(in) :: U
    real*8, intent(out) :: alpha,G,qptot

    real*8 :: g1,cp,Pr
    real*8 :: cf_codo,kappa
    real*8 :: pres,temp,rho,vel,dia,mu,Re,ff,cf_rect,cf

    g1 = ga-1.
    cp = R_gas*ga/g1

    G     = 0.0d0
    qptot = 0.0d0

    rho  = U(1)
    vel  = U(2)
    pres = U(3)

    alpha = dareax/area
    
    Temp = pres/rho/R_gas

    dia = dsqrt(4.0d0*area/pi)
    ! dynamic viscosity  in [pascal sec]
    mu = compute_visco(Temp)
    Re = dabs(rho*vel) * dia / mu
    Re = dmax1(Re,1.0d-10)

    if(viscous_flow.eq.1) then
       cf_codo = 0.0d0
       ff      = 0.2373 / Re**0.25
       cf_rect = 4.*ff /dia
       cf      = cf_codo+cf_rect
       
       ! friction source
       G = cf/2.0*vel*dabs(vel)
    endif

    if(heat_flow.eq.1) then
       ! Pruebo la correlacion que puse en la tesis
       Pr    = 0.71
       kappa = mu*cp/Pr
       !  heat transfer source
       qptot = pi*0.0297*kappa*Re**0.75*Pr**(1./3.)*(Twall-Temp)
    endif

  end subroutine source_char

  function compute_cp(T,isp,R_gas)
    !
    !  Computes the Cp for isp specie at T Kelvin degrees
    !    (see pp. 130 Heywood for details)
    !
    !  We begin only with one specie (Nitrogen)
    !
    implicit none

    integer, intent(in) :: isp
    real*8, intent(in) :: T,R_gas
    real*8 :: compute_cp

    real*8, dimension(5) :: C

    if(T.lt.1000) then
       C(1) =  3.6748
       C(2) = -1.2082e-3
       C(3) =  2.3240e-6
       C(4) = -6.3218e-10
       C(5) = -2.2577e-13 
    else
       C(1) =  2.8963
       C(2) =  1.5155e-3
       C(3) = -5.7235e-7
       C(4) =  9.9807e-11
       C(5) = -6.5224e-15
    endif

    compute_cp = R_gas * (C(1)+C(2)*T+ &
         C(3)*T**2+C(4)*T**3+C(5)*T**4)

  end function compute_cp

  function compute_visco(T)
    !
    !  Computes the viscosity coefficient taken from 1st homework 
    !    of Abraham course
    !
    implicit none

    real*8, intent(in) :: T
    real*8 :: compute_visco

    real*8, dimension(2) :: Co_mu

    Co_mu(1) = 1.457d-6
    Co_mu(2) = 110.0d0

    ! dynamic viscosity  [Pascal sec]
    compute_visco = Co_mu(1)*T**(3./2.)/(Co_mu(2)+T)

  end function compute_visco

end module gasdyn_utils
