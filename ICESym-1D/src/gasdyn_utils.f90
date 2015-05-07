module gasdyn_utils

  use utilities, only : pi

  type :: gas_properties
     real*8 :: R_gas,cp,cv,gamma
     real*8, dimension(11) :: composition
  end type gas_properties

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

  subroutine compute_composition(burned_gas_model,phi,alpha, &
       beta,gamma,delta,mass_fuel,mass_air,mass_res,X)
    !
    !
    !  species CO2 H2O N2 O2 CO H2 H O OH and NO
    !
    implicit none

    integer, intent(in) :: burned_gas_model
    real*8, intent(in) :: alpha,beta,gamma,delta
    real*8, intent(in) :: phi,mass_air,mass_fuel,mass_res
    real*8, dimension(11), intent(out) :: X

    integer :: i
    real*8 :: nb,mRP,Mb,n_res,nT,psi,zetha,phi_c,z,eps
    real*8 :: A,B,C,K,nCO,Mw_fuel,y
    real*8, dimension(11) :: nmol

    do i=1,11
       X(i)    = 0.0d0
       nmol(i) = 0.0d0
    end do

    select case (burned_gas_model)
    case(0) ! composition for low-temperature burned gas
            ! see Heywood, pp. 102-107
       y = beta/alpha
       eps = 4.0d0/(4.0d0+y)
       z = gamma/alpha
       zetha = 2.0d0/(2.0d0-eps*z*(1.0d0-phi))
       psi   = (1.0d0-eps*z/2.0d0)*zetha*3.773d0
       phi_c = zetha*phi
       ! burned gas
       if(phi_c.le.1.0d0) then
          nb  = (1.0d0-eps)*phi_c+1.0d0+psi
          mRP = 32.0d0+4.0d0*phi_c*(1.0d0+2.0d0*eps)+28.16d0*psi
          Mb  = mRP/nb
          n_res = mass_res/Mb
          nmol(1) = n_res/nb*eps*phi_c
          nmol(2) = n_res/nb*2.0d0*phi_c*(1.0d0-eps)
          nmol(3) = n_res/nb*psi
          nmol(4) = n_res/nb*(1.0d0-phi_c)
       else
          nb  = (2.0d0-eps)*phi_c+psi
          mRP = 32.0d0+4.0d0*phi_c*(1.0d0+2.0d0*eps)+28.16d0*psi
          Mb  = mRP/nb
          n_res = mass_res/Mb

          K = 3.5d0
          A = K-1.0d0
          B = -(K*(2.0d0*(phi_c-1.0d0)+eps*phi_c)+2.0d0*(1.0d0-eps*phi_c))
          C = 2.0d0*K*eps*phi_c*(phi_c-1.0d0)
          nCO = (-B-dsqrt(B**2-4.0d0*A*C))/2.0d0/A

          nmol(1) = n_res/nb*(eps*phi_c-nCO)
          nmol(2) = n_res/nb*(2.0d0*(1.0d0-eps*phi_c)+nCO)
          nmol(3) = n_res/nb*psi
          nmol(5) = n_res/nb*nCO
          nmol(6) = n_res/nb*(2.0d0*(phi_c-1.0d0)-nCO)
       end if
       ! air
       nmol(3) = nmol(3)+mass_air/28.962d0*0.79d0
       nmol(4) = nmol(4)+mass_air/28.962d0*0.21d0
       ! fuel
       Mw_fuel = 12.01d0*alpha+1.008d0*beta+16.0d0*gamma+14.01d0*delta
       nmol(11) = mass_fuel/Mw_fuel

       nT = sum(nmol)
       do i=1,11
          X(i) = nmol(i)/nT
       end do
    case(1) ! chemical equilibrium
       stop ' Not implemented yet '
    case(2) ! user-defined
       stop ' Not implemented '
    case default
       stop ' Wrong burned gas model option'
    end select

  end subroutine compute_composition

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

    if(isp.ne.0) write(*,*) ' WARNING! In cp computation, isp!=0 '

    if(T.lt.1000.) then
       C(1) =  3.6748d0
       C(2) = -1.2082d-3
       C(3) =  2.3240d-6
       C(4) = -6.3218d-10
       C(5) = -2.2577d-13 
    else
       C(1) =  2.8963d0
       C(2) =  1.5155d-3
       C(3) = -5.7235d-7
       C(4) =  9.9807d-11
       C(5) = -6.5224d-15
    endif

    compute_cp = R_gas * (C(1)+C(2)*T+ &
         C(3)*T**2+C(4)*T**3+C(5)*T**4)

  end function compute_cp

  function compute_Rgas(X,Xfuel,Mfuel)
    !
    !  Computes the particular gas constant for an ideal gas mixture
    !  with molar composition X and fuel vapor with Xfuel molar 
    !  concentration
    !
    !  species CO2 H2O N2 O2 CO H2 H O OH and NO
    !
    implicit none

    real*8, intent(in) :: Xfuel,Mfuel
    real*8, intent(in), dimension(10) :: X
    real*8 :: compute_Rgas

    integer :: i
    real*8 :: MW,Ru
    real*8, dimension(10) :: M

    Ru = 8314.34 ! J/kmol*K

    M = (/44.01, 18.02, 28.008, 32.00, 28.01, 2.018, 1.009, &
         16.00, 17.009, 30.004/)

    MW = Xfuel*Mfuel
    do i=1,10
       MW = MW+X(i)*M(i)
    end do

    compute_Rgas = Ru/MW

  end function compute_Rgas

  function compute_cp_mix(T,R_gas,X,Xfuel,fuel_coef_cp)
    !
    !  Computes the specific heat at constant pressure for an ideal
    !  gas mixture with molar composition X and fuel vapor with 
    !  Xfuel molar concentration, at T Kelvin degrees
    !
    !  species CO2 H2O N2 O2 CO H2 H O OH and NO
    !
    implicit none

    real*8, intent(in) :: T,R_gas,Xfuel
    real*8, intent(in), dimension(5) :: fuel_coef_cp
    real*8, intent(in), dimension(10) :: X
    real*8 :: compute_cp_mix

    integer :: i
    real*8 :: cp,cp_fuel
    real*8, dimension(10,5) :: C

    if(T.lt.1000.) then
       ! CO2
       C(1,:) = (/2.4007797d0, 0.87350957d-2, -0.66070878d-5, 0.20021861d-8, &
            0.63274039d-15/)
       ! H20
       C(2,:) = (/4.0701275d0, -0.11084499d-2, 0.41521180d-5, -0.29637404d-8, &
            0.80702103d-12/)
       ! N2
       C(3,:) = (/3.6748d0, -1.2082d-3, 2.3240d-6, -6.3218d-10, -2.2577d-13/)
       ! O2
       C(4,:) = (/3.6255985d0, -0.18782184d-2, 0.70554544d-5, -0.67635137d-8, &
            0.21555993d-11/)
       ! CO
       C(5,:) = (/3.7100928d0, -0.16190964d-2, 0.36923594d-5, -0.20319674d-8, &
            0.23953344d-12/)
       ! H2
       C(6,:) = (/3.0574451d0, 0.26765200d-2, -0.58099162d-5, 0.55210391d-8, &
            -0.18122739d-11/)
       ! H
       C(7,:) = (/2.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
       ! O
       C(8,:) = (/2.9464287d0, -0.16381665d-2, 0.24210316d-5, -0.16028432d-8, &
            0.38906964d-12/)
       ! OH
       C(9,:) = (/3.8375943d0, -0.10778858d-2, 0.96830378d-6, 0.18713972d-9, &
            -0.22571094d-12/)
       ! NO
       C(10,:) = (/4.0459521d0, -0.34181783d-2, 0.79819190d-5, -0.61139316d-8, &
            0.15919076d-11/)
    else
       ! CO2
       C(1,:) = (/4.4608041d0, 0.30981719d-2, -0.12392571d-5, 0.22741325d-9, &
            -0.15525954d-13/)
       ! H20
       C(2,:) = (/2.7167633d0, 0.29451374d-2, -0.80224374d-6, 0.10226682d-9, &
            -0.48472145d-14/)
       ! N2
       C(3,:) = (/2.8963194d0, 0.15154866d-2, -0.57235277d-6, 0.99807393d-10, &
            -0.65223555d-14/)
       ! O2
       C(4,:) = (/3.6219535d0, 0.73618264d-3, -0.19652228d-6, 0.36201558d-10, &
            -0.28945627d-14/)
       ! CO
       C(5,:) = (/2.9840696d0, 0.14891390d-2, -0.57899684d-6, 0.10364577d-9, &
            -0.69353550d-14/)
       ! H2
       C(6,:) = (/3.1001901d0, 0.51119464d-3, 0.52644210d-7, -0.34909973d-10, &
            0.36945345d-14/)
       ! H
       C(7,:) = (/2.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
       ! O
       C(8,:) = (/2.5420596d0, -0.27550619d-4, -0.31028033d-8, 0.45510674d-11, &
            -0.43680515d-15/)
       ! OH
       C(9,:) = (/2.9106427d0, 0.95931650d-3, -0.19441702d-6, 0.13756646d-10, &
            0.14224542d-15/)
       ! NO
       C(10,:) = (/3.189d0, 0.13382281d-2, -0.52899318d-6, 0.95919332d-10, &
            -0.64847932d-14/)
    endif

    cp = 0.0d0
    do i=1,10
       cp = cp+X(i)*(C(i,1)+C(i,2)*T+ &
            C(i,3)*T**2+C(i,4)*T**3+C(i,5)*T**4)
    end do

    cp_fuel = fuel_coef_cp(1)+fuel_coef_cp(2)*T+fuel_coef_cp(3)*T**2+ &
         fuel_coef_cp(4)*T**3+fuel_coef_cp(5)/T**2

    compute_cp_mix = R_gas*(cp+Xfuel*cp_fuel)

  end function compute_cp_mix

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
