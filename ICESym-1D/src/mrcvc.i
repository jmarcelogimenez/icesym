!
!  Includes the geometry, heat transfer correlations and ideal 
!  cycle states for the MRCVC engine
!

subroutine geometry_mrcvc(myData, globalData, Vol, Area, Vdot)
  !
  !  Computes the chamber volume and wall area of
  !    the MRCVC engine, considering "sharp" geometry
  !
  !  geometry_mrcvc is called by: geometry
  !  geometry_mrcvc calls the following subroutines and functions:
  !    region1_geo, region2_geo, region3_geo, region4_geo, region5_geo
  !
  implicit none

  type(this), intent(in) :: myData
  type(dataSim), intent(in) :: globalData
  real*8, intent(out) :: Vol,Area,Vdot

  integer :: n
  real*8 :: R,r_,h,R_i,R_e,alpha,beta,delta
  real*8 :: theta_cycle,rpm,theta_g,theta,dVol
  real*8 :: dalpha,dbeta,dgamma,epsi,epsr,l1,l2,l3,s
  real*8 :: ddalpha,ddbeta,depsi,depsr,dl1,dl2,dl3,ds

  theta_cycle = globalData%theta_cycle
  rpm         = globalData%rpm
  theta_g     = globalData%theta

  R  = myData%major_radius  ! Radius of the vanes' center (R)
  r_ = myData%minor_radius  ! Radius of the vanes (r)
  h  = myData%chamber_heigh ! Chamber heigh (h)
  n  = myData%nvanes        ! Number of vanes (n)

  theta = modulo(theta_g+myData%theta_0-(n+2.)/(2.*n)*pi, &
       theta_cycle)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)
  delta = pi/n-alpha

  if(theta.ge.0. .and. theta.lt.(n-2.)/(2.*n)*pi-alpha) then
     call region1_geo(theta, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(n-2.)/(2.*n)*pi) then
     call region2_geo(theta, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(n-2.)/(2.*n)*pi+alpha) then
     call region3_geo(theta, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(n+2.)/(2.*n)*pi-alpha) then
     call region4_geo(theta, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(n+2.)/(2.*n)*pi+alpha) then
     call region5_geo(theta, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(n+6.)/(2.*n)*pi-alpha) then
     call region4_geo((n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(n+6.)/(2.*n)*pi) then
     call region3_geo((n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(n+6.)/(2.*n)*pi+alpha) then
     call region2_geo((n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(n+2.)/n*pi) then
     call region1_geo((n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(3.*n+2.)/(2.*n)*pi-alpha) then
     call region1_geo(theta-(n+2.)/n*pi, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(3.*n+2.)/(2.*n)*pi) then
     call region2_geo(theta-(n+2.)/n*pi, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(3.*n+2.)/(2.*n)*pi+alpha) then
     call region3_geo(theta-(n+2.)/n*pi, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.3.*(n+2.)/(2.*n)*pi-alpha) then
     call region4_geo(theta-(n+2.)/n*pi, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.3.*(n+2.)/(2.*n)*pi+alpha) then
     call region5_geo(theta-(n+2.)/n*pi, n, R, r_, h, Vol, Area, dVol)
  elseif(theta.lt.(3.*n+10.)/(2.*n)*pi-alpha) then
     call region4_geo(2.*(n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(3.*n+10.)/(2.*n)*pi) then
     call region3_geo(2.*(n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.lt.(3.*n+10.)/(2.*n)*pi+alpha) then
     call region2_geo(2.*(n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  elseif(theta.le.2.*(n+2.)/n*pi) then
     call region1_geo(2.*(n+2.)/n*pi-theta, n, R, r_, h, Vol, Area, dVol)
     dVol = -1.0d0 * dVol
  endif

  Area = Area + 2.*Vol/h

  Vdot = rpm*pi/30.*dVol

end subroutine geometry_mrcvc

subroutine region1_geo(theta, n, R, r_, h, Vol, Area, dVol)
  !
  !  Computes the chamber volume and wall surface area for
  !  region 1
  !
  !  region1_geo is called by: geometry_mrcvc, ideal_cycle_mrcvc
  !  region1_geo calls the following subroutines and functions:
  !    none
  !
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: theta,R,r_,h
  real*8, intent(out) :: Vol,Area,dVol

  real*8 :: R_i,R_e,alpha,beta,delta,pi

  pi = 4.*datan(1.d0)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)
  delta = pi/n-alpha

  Area = 2.*h*(2.*(alpha*R+alpha*r_+pi*R/n)+delta*R_i+beta*R_e)

  Vol  = h*(4.*R*r_*dsin(pi/n)*dcos(theta)+R_i*(2.*r_-3.*alpha*R_i)- &
       4.*pi/n*r_**2+R_e**2*dacos((R**2+r_**2)/(R*R_e)))

  dVol = -4.*h*R*r_*dsin(pi/n)*dsin(theta)
 
end subroutine region1_geo

subroutine region2_geo(theta, n, R, r_, h, Vol, Area, dVol)
  !
  !  Computes the chamber volume and wall surface area for
  !  region 2
  !
  !  region2_geo is called by: geometry_mrcvc
  !  region2_geo calls the following subroutines and functions:
  !    none
  !
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: theta,R,r_,h
  real*8, intent(out) :: Vol,Area,dVol

  real*8 :: R_i,R_e,alpha,beta,delta,pi
  real*8 :: dalpha,dbeta,dgamma,epsi,epsr,l1,l2,l3,s
  real*8 :: ddalpha,ddbeta,depsi,depsr,dl1,dl2,dl3,ds

  pi = 4.*datan(1.d0)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)
  delta = pi/n-alpha

  Area = h*((R+r_)*((n-2.)*pi/(2.*n)+3.*alpha-theta)+ &
       2.*R*((n+2.)*pi/(2.*n)-alpha-theta)+2.*delta*R_i- &
       R_e*(alpha+theta-2.*beta-(n-2.)*pi/(2.*n)))

  dgamma = (n-2.)/(2.*n)*pi-theta
  s      = dsqrt(2.*R*(R+R_i*dcos(dgamma))-r_**2)
  epsi   = dasin(dsqrt(4.*R**2*s**2-(s**2+r_**2)**2)/(2.*R_i*s))
  epsr   = dasin(dsqrt(4.*R**2*s**2-(s**2+r_**2)**2)/(2.*R*s))
  dalpha = 2.*epsi-dgamma
  dbeta  = dgamma-2.*epsr
  l1     = 2.*R*dsin(.5*(alpha+dbeta))
  l2     = 2.*r_*dsin(.5*(alpha+dgamma))
  l3     = 2.*R_i*dsin(.5*(alpha-dalpha))
  s      = .5*(l1+l2+l3)

  Vol  = h*(2.*r_*(R_i+R*dsin(pi/n-theta))+R_e**2/2.*datan(2.*r_/R_i)- &
       2.*r_**2*((n+2.)/(2.*n)*pi-theta)+.5*(R_e**2*dacos((R**2+r_**2)/(R*R_e))- &
       3.*alpha*R_i**2)+dsqrt(s*(s-l1)*(s-l2)*(s-l3))+ &
       .5*(r_**2*(alpha+dgamma-dsin(alpha+dgamma))-R**2*(alpha+dbeta-dsin(alpha+dbeta))- &
       R_i**2*(alpha-dalpha-dsin(alpha-dalpha)))-R_i**2*(epsi-.5*dsin(2.*epsi))- &
       R**2*(epsr-.5*dsin(2.*epsr)))

  s       = dsqrt(2.*R*(R+R_i*dcos(dgamma))-r_**2)
  ds      = R*R_i*dsin((n-2.)/(2.*n)*pi-theta)/s
  depsi   = ds*(r_**4-s**4)/(2.*R_i*s**2*dsqrt(4.*R**2*s**2-(s**2+r_**2)**2)*dcos(epsi))
  depsr   = ds*(r_**4-s**4)/(2.*R*s**2*dsqrt(4.*R**2*s**2-(s**2+r_**2)**2)*dcos(epsr))
  ddalpha = 2.*depsi+1.
  ddbeta  = -1.-2.*depsr
  dl1     = R*ddbeta*dcos(.5*(alpha+dbeta))
  dl2     = -r_*dcos(.5*(alpha+dgamma))
  dl3     = -R_i*ddalpha*dcos(.5*(alpha-dalpha))
  s       = .5*(l1+l2+l3)
  ds      = .5*(dl1+dl2+dl3)

  dVol = h*(-2.*r_*R*dcos(pi/n-theta)+2.*r_**2+ &
       .5*dsqrt(s*(s-l1)*(s-l2)*(s-l3))*(ds/s+(ds-dl1)/(s-l1)+(ds-dl2)/(s-l2)+ &
       (ds-dl3)/(s-l3))-.5*r_**2*(1.-dcos(alpha+dgamma))- &
       .5*R**2*ddbeta*(1.-dcos(alpha+dbeta))+.5*R_i**2*ddalpha*(1.-dcos(alpha-dalpha))- &
       R_i**2*depsi*(1.-dcos(2.*epsi))-R**2*depsr*(1.-dcos(2.*epsr)))
  
end subroutine region2_geo

subroutine region3_geo(theta, n, R, r_, h, Vol, Area, dVol)
  !
  !  Computes the chamber volume and wall surface area for
  !  region 3
  !
  !  region3_geo is called by: geometry_mrcvc
  !  region3_geo calls the following subroutines and functions:
  !    none
  !
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: theta,R,r_,h
  real*8, intent(out) :: Vol,Area,dVol

  real*8 :: R_i,R_e,alpha,beta,delta,pi
  real*8 :: dalpha,l1,l2,l3,s
  real*8 :: dl1,dl2,dl3,ds

  pi = 4.*datan(1.d0)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)
  delta = pi/n-alpha

  Area = h*((R+r_)*((n-2.)*pi/(2.*n)+3.*alpha-theta)+ &
       2.*R*((n+2.)*pi/(2.*n)-alpha-theta)+2.*delta*R_i- &
       R_e*(alpha+theta-2.*beta-(n-2.)*pi/(2.*n)))
  
  dalpha = (n-2.)/(2.*n)*pi+alpha-theta
  l1     = 2.*R*dsin(.5*dalpha)
  l2     = 2.*r_*dsin(.5*dalpha)
  l3     = 2.*R_i*dsin(.5*dalpha)
  s      = .5*(l1+l2+l3)
  
  Vol  = h*(2.*r_*(R_i+R*dsin(pi/n-theta))+R_e**2/2.*datan(2.*r_/R_i)- &
       2.*r_**2*((n+2.)/(2.*n)*pi-theta)+.5*(R_e**2*dacos((R**2+r_**2)/(R*R_e))- &
       3.*alpha*R_i**2)+dsqrt(s*(s-l1)*(s-l2)*(s-l3))-R_i**2*(dalpha-dsin(dalpha)))
  
  dl1 = -R*dcos(.5*dalpha)
  dl2 = -r_*dcos(.5*dalpha)
  dl3 = -R_i*dcos(.5*dalpha)
  ds  = .5*(dl1+dl2+dl3)
  
  dVol = h*(-2.*r_*R*dcos(pi/n-theta)+2.*r_**2+ &
       .5*dsqrt(s*(s-l1)*(s-l2)*(s-l3))*(ds/s+(ds-dl1)/(s-l1)+(ds-dl2)/(s-l2)+ &
       (ds-dl3)/(s-l3))+R_i**2*(1.-dcos(dalpha)))

end subroutine region3_geo

subroutine region4_geo(theta, n, R, r_, h, Vol, Area, dVol)
  !
  !  Computes the chamber volume and wall surface area for
  !  region 4
  !
  !  region4_geo is called by: geometry_mrcvc
  !  region4_geo calls the following subroutines and functions:
  !    none
  !
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: theta,R,r_,h
  real*8, intent(out) :: Vol,Area,dVol

  real*8 :: R_i,R_e,alpha,beta,pi

  pi = 4.*datan(1.d0)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)

  Area = h*(2.*alpha*(R+r_)+ (2.*R+R_i)*((n+2.)*pi/(2.*n)-alpha-theta)- &
       R_e*(alpha+theta-2.*beta-(n-2.)*pi/(2.*n)))

  Vol  = h*(2.*r_*(R_i+R*dsin(pi/n-theta))+R_e**2/2.*datan(2.*r_/R_i)- &
       2.*r_**2*((n+2.)/(2.*n)*pi-theta)+.5*(R_e**2*dacos((R**2+r_**2)/(R*R_e))- &
       3.*alpha*R_i**2))

  dVol = h*(-2.*r_*R*dcos(pi/n-theta)+2.*r_**2)

end subroutine region4_geo

subroutine region5_geo(theta, n, R, r_, h, Vol, Area, dVol)
  !
  !  Computes the chamber volume and wall surface area for
  !  region 5
  !
  !  region5_geo is called by: geometry_mrcvc, ideal_cycle_mrcvc
  !  region5_geo calls the following subroutines and functions:
  !    none
  !
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: theta,R,r_,h
  real*8, intent(out) :: Vol,Area,dVol

  real*8 :: R_i,R_e,alpha,beta,pi

  pi = 4.*datan(1.d0)

  R_i = dsqrt(R**2-r_**2)
  R_e = dsqrt(R**2+3.*r_**2)

  alpha = dasin(r_/R)
  beta  = pi/n+alpha-datan(2.*r_/R_i)

  Area = 2.*h*((R+r_)*alpha+R_e*(beta-pi/n))

  Vol  = h*(2.*r_*(R_i-R*dsin(pi/2.-alpha))+R_e**2/2.*datan(2.*r_/R_i)- &
       2.*r_**2*alpha+.5*(R_e**2*dacos((R**2+r_**2)/(R*R_e))- &
       3.*alpha*R_i**2))

  dVol = 0.

end subroutine region5_geo

subroutine heat_transfer_mrcvc(myData, globalData, Ucyl, Area, cp, dQ_ht)
  !
  !  Computes the heat losses through the chamber walls
  !
  !  model = 1 -> Annand model
  !  model = 2 -> Woschni model 1
  !  model = 3 -> Woschni model 2
  !  model = 4 -> Wilmer
  !
  !  heat_transfer_mrcvc is called by: heat_transfer
  !  heat_transfer_mrcvc calls the following subroutines and functions:
  !    none
  !  
  implicit none
  
  real*8, intent(in) :: Area,cp
  real*8, dimension(3), intent(in) :: Ucyl
  type(this), intent(in) :: myData
  type(dataSim), intent(in) :: globalData 
  real*8, intent(out) :: dQ_ht
 
  integer :: ispecie,n
  real*8 :: R,r_,h,Twall,rho,p,T
  real*8 :: mu,Pr,kappa,Sp,Re,Nu
  real*8 :: factor,C_h,C_r,dQ_hth,dQ_htr
  
  R     = myData%major_radius  ! Radius of the vanes' center (R)
  r_    = myData%minor_radius  ! Radius of the vanes (r)
  h     = myData%chamber_heigh ! Chamber heigh (h)
  n     = myData%nvanes        ! Number of vanes (n)
  Twall = myData%Twall         ! Wall temperature of the chamber
  
  rho = Ucyl(1)
  p   = Ucyl(2)
  T   = Ucyl(3)
  
  ! computing Cp for a given mixture (isp) at a given temperature
  ! ispecie = 0
  ! cp = compute_cp(T,ispecie,globalData%R_gas)
 
  mu    = compute_visco(T) ! dynamic viscosity
  Pr    = 0.72             ! Prandtl number
  kappa = mu*cp/Pr         ! thermal conductivity
  
  ! average speed
  Sp  = 2.*pi*R*globalData%rpm/60.
  
  Re = rho*dabs(Sp)*2.*r_/mu    ! Reynolds number

  ! Nusselt number
  select case (myData%model_ht)  
  case (1)
     ! According to Annand (J.I. Ramos pp. 25)
     factor = myData%factor_ht * 0.14
     Nu     = factor * Re**0.7
  case (2)
     ! factor is a constant times 0.035
     ! Some authors use a sensibility analysis in order to set this parameter.
     ! For instance, John Abraham uses in his code the value 4.5 times 0.035;
     ! Yacoub and Bata apply between 1 and 3 times this value in their
     ! '98 paper.
     factor = myData%factor_ht * 0.035
     Nu     = factor*Re**0.8*Pr**0.33
  case (3)
     factor = myData%factor_ht * 0.037
     Nu     = factor*Re**0.8*Pr**0.3
  case (4)
     ! Correlacion by Wilmer (J.I. Ramos pp. 25)
     factor = myData%factor_ht*685.2497
     Nu     = factor*(globalData%rpm/1000.*p/1.d5)**0.786*T**(-0.525)*2.*r_/kappa
  end select
  
  ! heat transfer coefficient
  C_h = Nu*kappa/(2.*r_)
  C_r = 0.0d0
  if(myData%type_ig.eq.0) then
     ! SI Engine
     C_r = 4.25d-9*((T**2+Twall**2)*(T+Twall))
  elseif(myData%type_ig.eq.1) then
     ! CI Engine
     C_r = 3.2602d-8*((T**2+Twall**2)*(T+Twall))
  endif
  
  dQ_hth = Area*C_h*(T-Twall)
  dQ_htr = Area*C_r*(T-Twall)
  dQ_ht = dQ_hth + dQ_htr
  if(globalData%save_extras) then
     write(myData%nunit,902) C_h, C_r, dQ_hth, dQ_htr
902  format (F12.4,F12.4,E15.6,E15.6)
  endif

end subroutine heat_transfer_mrcvc

subroutine ideal_cycle_mrcvc(myData, globalData, Q_LHV, T_i, p_i, AF, AVPT)
  !
  !  AVPT = [angle , density , pressure , temperature ] 
  !
  !  ideal_cycle_mrcvc is called by: run_ideal_cycle
  !  ideal_cycle_mrcvc calls the following subroutines and functions: 
  !    region1_geo, region5_geo
  !
  implicit none

  type(this), intent(in) :: myData
  type(dataSim), intent(in) :: globalData
  real*8, intent(in) :: Q_LHV,T_i,p_i,AF
  real*8, dimension(7,4), intent(out) :: AVPT

  real*8 :: cv,R_gas,ga,Vmin,Vmax,Amax,Amin,dVol,r_c
  real*8 :: theta1,theta2,theta3,theta4,theta5,theta6,theta7
  real*8 :: rho1,rho2,rho3,rho4,rho5,rho6,rho7
  real*8 :: V1,V2,V3,V4,V5,V6,V7
  real*8 :: p1,p2,p3,p4,p5,p6,p7
  real*8 :: T1,T2,T3,T4,T5,T6,T7
  real*8 :: m1,m2,m3,m4,m5,m6,m7,mf

  integer :: n
  real*8 :: R,r_,h

  ga    = globalData%ga
  R_gas = globalData%R_gas
  cv    = R_gas/(ga-1.)

  R  = myData%major_radius  ! Radius of the vanes' center (R)
  r_ = myData%minor_radius  ! Radius of the vanes (r)
  h  = myData%chamber_heigh ! Chamber heigh (h)
  n  = myData%nvanes        ! Number of vanes (n)

  call region1_geo(0.d0, n, R, r_, h, Vmax, Amax, dVol)
  call region5_geo((n+2.)/(2.*n)*pi, n, R, r_, h, Vmin, Amin, dVol)

  ! point 1
  V1     = Vmin
  p1     = 0.99999*p_i
  T1     = T_i
  rho1   = p1/R_gas/T1
  m1     = rho1*V1
  theta1 = 0.

  ! point 2
  V2     = Vmax
  p2     = p1
  T2     = T1
  rho2   = rho1
  m2     = 0.7*rho2*V2
  theta2 = (n+2.)/(2.*n)*180.

  mf     = m2/AF
  ! point 3
  V3     = Vmin
  m3     = m2
  rho3   = m3/V3
  p3     = p2*(rho3/rho2)**ga
  T3     = p3/R_gas/rho3
  theta3 = ((n+2.)/n-dasin(r_/R)/pi)*180.

  ! point 4
  V4     = Vmin
  m4     = m3
  rho4   = rho3
  T4     = T3+mf*Q_LHV/m4/cv
  p4     = R_gas*rho4*T4
  theta4 = ((n+2.)/n+dasin(r_/R)/pi)*180.

  ! point 5
  V5     = Vmax
  m5     = m4
  rho5   = m5/V5
  p5     = p4*(rho5/rho4)**ga
  T5     = p5/R_gas/rho5
  theta5 = 3.*(n+2.)/(2.*n)*180.-1.

  ! point 6
  V6     = Vmax
  p6     = 1.00001*p_i
  r_c    = Vmax/Vmin
  T6     = T2*(1.+Q_LHV*mf/m2/cv/T2/r_c**(ga-1.))**(1./ga)
  T6     = 1.*T2
  rho6   = p6/(R_gas*T6)
  m6     = rho6*V6
  theta6 = 3.*(n+2.)/(2.*n)*180.+1.

 ! point 7
  V7     = V1
  p7     = p6
  T7     = T6
  rho7   = rho6
  m7     = rho7*V7
  theta7 = 2.*(n+2.)/n*180.

  AVPT(:,1) = (/theta1,theta2,theta3,theta4,theta5,theta6,theta7/)
  AVPT(:,2) = (/rho1,rho2,rho3,rho4,rho5,rho6,rho7/)
  AVPT(:,3) = (/p1,p2,p3,p4,p5,p6,p7/)
  AVPT(:,4) = (/T1,T2,T3,T4,T5,T6,T7/)

end subroutine ideal_cycle_mrcvc
