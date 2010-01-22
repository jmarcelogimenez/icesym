module def_tube

  use def_simulator
  use utilities, only : pi
  use, intrinsic :: ISO_C_BINDING
  
  type, BIND(C) :: this
     integer(C_INT) :: nnod, ndof,nnod_input
     real(C_DOUBLE) :: longitud, dt_max
  end type this

contains
   subroutine state_initial_tube(myData, atm, globalData, state_ini, xnod, Area, Twall, curvature, dAreax) BIND(C)

    use, intrinsic :: ISO_C_BINDING
    type(this)::myData
    type(dataSim)::globalData
    real(C_DOUBLE) :: state_ini(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE), dimension(myData%nnod) :: xnod, Area, Twall, curvature, dAreax
    real(C_DOUBLE) :: atm(0:2)

    do i=0,((myData%nnod*myData%ndof)-1),3
       state_ini(i)   = atm(0)
       state_ini(i+1) = atm(1)
       state_ini(i+2) = atm(2)
    enddo

   end subroutine state_initial_tube
	
  subroutine solve_tube(myData, globalData, state, new_state, xnod, hnod, &
       Area, Twall, curvature, dAreax, itube) BIND(C)

    use, intrinsic :: ISO_C_BINDING
    type(this)::myData
    type(dataSim)::globalData
    integer(C_INT) :: itube
    real(C_DOUBLE) :: state(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: new_state(0:(myData%nnod*myData%ndof)-1)
    real(C_DOUBLE), dimension(myData%nnod-1) :: hnod
    real(C_DOUBLE), dimension(myData%nnod) :: xnod, Area, Twall, curvature, dAreax
    integer :: i,j
    real*8 :: ga,dt,dx
    real*8, dimension(5) :: prop_g
    real*8, dimension(myData%nnod-1) :: tau,tau1
    real*8, dimension(3,myData%nnod) :: U,Up,U_old,Fa,H,F_TVD
    real*8, dimension(3,myData%nnod_input) :: Uref
    
    ga = globalData%ga
    dt = globalData%dt

    prop_g(2) = ga
    prop_g(3) = globalData%viscous_flow
    prop_g(4) = globalData%heat_flow
    prop_g(5) = globalData%R_gas

    tau  =  dt / hnod
    tau1 = tau / hnod

    do i=0,myData%nnod-1
       do j=0,myData%ndof-1
          Up(j+1,i+1) = state(i*myData%ndof+j)
       enddo
    enddo

    do i=myData%nnod,myData%nnod+1
       do j=0,myData%ndof-1
          Uref(j+1,i-myData%nnod+1) = state(i*myData%ndof+j)
       enddo
    enddo

    ! From primitive basis to conservative basis
    call pv2cv(Up, U, ga, myData%nnod)
    U_old = U

    forall(i = 1:myData%nnod, j = 1:3)
       U(j,i) = U(j,i)*Area(i)
    end forall

    call fluxa(U, Area, ga, myData%nnod, 0, Fa)
    call source(U, dAreax, Area, prop_g, Twall, myData%nnod, H, 0)
    call fluxTVD(U, Fa, H, dt, tau, ga, myData%nnod, F_TVD)

    forall(i = 2:myData%nnod-1, j = 1:3)
       U(j,i)=U(j,i)-tau(i)*(F_TVD(j,i)-F_TVD(j,i-1))+dt*H(j,i)
    end forall

    U(:,1)           = U(:,2)
    U(:,myData%nnod) = U(:,myData%nnod-1)

    forall(i = 1:myData%nnod, j = 1:3)
       U(j,i) = U(j,i)/Area(i)
    end forall

    ! apply boundary conditions
    dx = hnod(1)
    call bc_tubes(U(:,1), Uref(:,1), U_old(:,2), ga, -1, dt, dx)
    dx = hnod(myData%nnod-1)
    call bc_tubes(U(:,myData%nnod), Uref(:,2), &
         U_old(:,myData%nnod-1), ga, 1, dt, dx)
    
    ! From conservative basis to primitive basis
    call cv2pv(U, Up, ga, myData%nnod)

    do i=0,myData%nnod-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = Up(j+1,i+1)
       enddo
    enddo

    ! myData%dt_max = dt
    call critical_dt(globalData, Up, hnod, myData%nnod, myData%dt_max)

    call actualize_valves(itube, Area, Twall, dAreax, myData%nnod)
    
  end subroutine solve_tube
  
  subroutine cv2pv(Uc, Up, ga, nnod)
    !
    !  Conversion from primitive variables to conservative variables
    !
    !  cv2pv is called by: solve_tube
    !  cv2pv calls the following subroutines and functions: none
    !
    implicit none
    
    integer, intent(in) :: nnod
    real*8, intent(in) :: ga
    real*8, dimension(3,nnod), intent(in) :: Uc
    real*8, dimension(3,nnod), intent(out) :: Up

    real*8 :: g1
    
    g1 = ga-1.
    Up(1,:) = Uc(1,:)
    Up(2,:) = Uc(2,:)/Uc(1,:)
    Up(3,:) = (Uc(3,:)-0.5*Up(1,:)*Up(2,:)**2)*g1
    
  end subroutine cv2pv
  
  subroutine pv2cv(Up, Uc, ga, nnod)
    !
    !  Conversion from conservativa variables to primitive variables
    !
    !  pv2cv is called by: solve_tubes, bc_tubes
    !  pv2cv calls the following subroutines and functions: none
    !  
    implicit none
    
    integer, intent(in) :: nnod
    real*8, intent(in) :: ga
    real*8, dimension(3,nnod), intent(in) :: Up
    real*8, dimension(3,nnod), intent(inout) :: Uc

    real*8 :: g1

    g1 = ga-1
    Uc(1,:) = Up(1,:)
    Uc(2,:) = Up(2,:)*Up(1,:)
    Uc(3,:) = Up(3,:)/g1+0.5*Up(1,:)*Up(2,:)**2

  end subroutine pv2cv

  subroutine fluxa(U, area, ga, nnod, scheme, F)
    !
    !  Computes advective flux
    !
    !  fluxa is called by: solve_tube
    !  fluxa calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: nnod,scheme
    real*8, intent(in) :: ga
    real*8, dimension(3,nnod), intent(in) :: U
    real*8, dimension(nnod), intent(in) :: area
    real*8, dimension(3,nnod), intent(out) :: F

    real*8 :: g1
    real*8, dimension(nnod) :: pres

    g1 = ga-1.
    
    if(scheme.eq.0)then
       pres = (U(3,:)-0.5*U(2,:)**2/U(1,:))*g1/area
       F(1,:) = U(2,:)
       F(2,:) = U(2,:)**2/U(1,:)+pres*area
       F(3,:) = (U(3,:) + pres*area)*(U(2,:)/U(1,:))
    else
       pres = (U(3,:)-0.5*U(2,:)**2/U(1,:))*g1
       F(1,:) = U(2,:)
       F(2,:) = U(2,:)**2/U(1,:)+pres
       F(3,:) = (U(3,:) + pres)*(U(2,:)/U(1,:))
    endif

  end subroutine fluxa

  subroutine source(U, dareax, area, prop_g, Twall, nnod, h, scheme)
    !
    !   Computes the source term for 
    !      variable section,
    !      friction flow, and
    !      heat transfer
    !
    !  source is called by: solve_tube
    !  source calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: nnod,scheme
    real*8, dimension(3,nnod), intent(in) :: U
    real*8, dimension(nnod), intent(in) :: dareax,area,Twall
    real*8, dimension(5), intent(in) :: prop_g
    real*8, dimension(3,nnod), intent(out) :: h

    integer :: viscous_flow,heat_flow
    real*8 :: co_mu(2),ga,g1,R_gas,cp,Pr
    real*8 :: th_equiv,cf_codo,alfa_ext,h_cond
    real*8, dimension(nnod) :: pres,temp,rho,rho_vel
    real*8, dimension(nnod) :: dia,mu,Re,ff,cf_rect,cf
    real*8, dimension(nnod) :: ggg,alfa,qptot,alfa_int,K,kappa

    viscous_flow = prop_g(3)
    heat_flow    = prop_g(4)
    ga           = prop_g(2)
    R_gas        = prop_g(5)

    g1 = ga-1.
    cp = R_gas*ga/g1

    h = 0.0d0
    if(scheme.eq.0) then
       pres    = (U(3,:)-0.5*U(2,:)**2/U(1,:))*g1/area
       rho     = U(1,:)/area
       rho_vel = U(2,:)/area
       h(2,:)  = pres*dareax
    else
       pres    = (U(3,:)-0.5*U(2,:)**2/U(1,:))*g1
       rho     = U(1,:)
       rho_vel = U(2,:)
       h(1,:)  = -dareax/area*U(2,:)
       h(2,:)  = -dareax/area*U(2,:)**2/U(1,:)
       h(3,:)  = -dareax/area*(U(3,:)+pres)*U(2,:)/U(1,:)
    endif

    Temp = pres/rho/R_gas

    dia = dsqrt(4.0d0*area/pi)
    ! viscosity coefficient
    co_mu(1) = 1.457e-6
    co_mu(2) = 110.0
    ! dynamic viscosity  in [pascal sec]
    mu = co_mu(1)*Temp**1.5/(co_mu(2)+Temp)
    Re = dabs(rho_vel) * dia / mu
    Re = dmax1(Re,1.0d-10)

    if(viscous_flow.eq.1) then
       cf_codo = 0.0d0
       ff      = 0.2373 / Re**0.25
       cf_rect = 4.*ff /dia
       cf      = cf_codo+cf_rect
       
       ! friction source
       h(2,:) = h(2,:) - cf/2.0*U(2,:)*dabs(U(2,:)/U(1,:))
    endif

    if(heat_flow.eq.1) then
       if(.false.) then
          ! Ecuaciones originales
          ggg      = dabs(rho_vel)*3600.0 ! en kg/m^2/hora
          alfa     = 4.0d-3*ggg**0.8 / dia**0.2
          th_equiv = 4187.0  ! 4187 joules = 1 kcal
          
          if(.true.)then
             alfa_int = alfa
             alfa_ext = 1.0
             h_cond   = 1.0
             K        = 1.0/(1.0/alfa_int + 1.0/h_cond + 1.0/alfa_ext)
             qptot    = 4*K*(Twall-Temp)/dia*area*Th_equiv
          else
             qptot    = 4*alfa*(Twall-Temp)/dia*area*th_equiv
          endif
          
          !  heat transfer source       
          h(3,:) = h(3,:)+qptot / 3600.0
       else
          ! Pruebo la correlacion que puse en la tesis
          Pr    = 0.71
          kappa = mu*cp/Pr
          qptot = pi*0.0297*kappa*Re**0.75*Pr**(1./3.)*(Twall-Temp)
          !  heat transfer source
          h(3,:) = h(3,:)+qptot
       endif
    endif

  end subroutine source

  subroutine fluxTVD(U, fa, c, dt, r, ga, nnod, f_tvd)
    !
    !
    !
    !  fluxTVD is called by: solve_tube
    !  fluxTVD calls the following subroutines and functions: diagsys, qtvd2, ajac
    !
    implicit none

    integer, intent(in) :: nnod
    real*8,intent(in) :: dt,ga
    real*8, dimension(3,nnod), intent(in) :: U,fa,c
    real*8, dimension(nnod), intent(in) :: r
    real*8, dimension(3,nnod-1), intent(out) :: f_tvd

    integer :: i,j,k
    real*8 :: tol
    real*8, dimension(3,3,nnod-1) :: p_half,p1_half,j_half
    real*8, dimension(3,nnod) :: deltau,g_tilde,g_2tilde,theta,g_i
    real*8, dimension(3,nnod-1) :: u_half,c_half,d_half,nu_half,alfa_half
    real*8, dimension(3,nnod-1) :: gamma_half,s_alfa_half,ter2,ter3,ter23,ter4
    real*8, dimension(3,nnod-1) :: q_half,q2_half,g_tilde_half,sigma_half
    real*8, dimension(3,nnod-1) :: s_half,vaux

    tol = 1d-12

    forall(i = 1:nnod-1, j = 1:3)
       deltau(j,i) = U(j,i+1)-U(j,i)
    end forall

    forall(i=1:nnod,j=1:3)
       deltau(j,i) = dsign(max(tol,dabs(deltau(j,i))),deltau(j,i))
    end forall

    forall(i = 1:nnod-1,j=1:3)
       u_half(j,i) = 0.5*(U(j,i+1)+U(j,i))
       c_half(j,i) = 0.5*(c(j,i+1)+c(j,i))
    end forall

    call diagsys(u_half, ga, nnod-1, p_half, p1_half, d_half)

    do j=1,3
       forall(i = 1:nnod-1)
          nu_half(j,i) = r(i)*d_half(j,i)
       end forall
    end do

    alfa_half = 0.0d0
    do k=1,3
       do j=1,3
          forall(i = 1:nnod-1)
             alfa_half(k,i) = alfa_half(k,i)+p1_half(k,j,i)*deltau(j,i)
          end forall
       end do
    end do

    call qtvd2(nu_half, 0*nu_half, nnod-1, q_half)

    g_tilde_half = 0.5d0*(q_half-nu_half**2)*alfa_half

    sigma_half  = 0.5d0*(1.0d0-q_half)
    s_alfa_half = dsign(1.0d0,alfa_half)
    s_half      = dsign(1.0d0,g_tilde_half)

    vaux   = 0.0d0
    do k=1,3
       forall (i = 2:nnod-1)
          vaux(k,i)    = s_half(k,i)*g_tilde_half(k,i-1)
          vaux(k,i)    = min(dabs(g_tilde_half(k,i)),vaux(k,i))
          vaux(k,i)    = max(vaux(k,i),0.0d0)
          g_tilde(k,i) = s_half(k,i)*vaux(k,i)
       end forall
    end do

    vaux   = 0.0d0
    do k=1,3
       forall (i = 2:nnod-1)
          vaux(k,i) = (s_alfa_half(k,i-1)*sigma_half(k,i-1)*alfa_half(k,i-1))
          vaux(k,i) = min(vaux(k,i),(sigma_half(k,i)*dabs(alfa_half(k,i))))
          vaux(k,i) = max(vaux(k,i),0.0d0)
          g_2tilde(k,i) = s_alfa_half(k,i)*vaux(k,i)
          
          theta(k,i) = (dabs(alfa_half(k,i)-alfa_half(k,i-1)) /   &
               (dabs(alfa_half(k,i))+dabs(alfa_half(k,i-1))))
       end forall

       g_tilde(k,1)     = g_tilde(k,2)
       g_tilde(k,nnod)  = g_tilde(k,nnod-1)
       g_2tilde(k,1)    = g_2tilde(k,2)
       g_2tilde(k,nnod) = g_2tilde(k,nnod-1)
       theta(k,1)       = theta(k,2)
       theta(k,nnod)    = theta(k,nnod-1)
    end do

    g_i = g_tilde + theta*g_2tilde

    do j=1,3
       do i=1,nnod-1
          if (dabs(alfa_half(j,i)) > tol) then
             alfa_half(j,i)  = alfa_half(j,i)
             gamma_half(j,i) = (g_tilde(j,i+1)-g_tilde(j,i))/alfa_half(j,i)
          else
             alfa_half(j,i) = tol
          endif
       end do
    end do

    call ajac(u_half, ga, nnod-1, j_half)

    call qtvd2(nu_half, gamma_half, nnod-1, q2_half)

    do j=1,3
       forall(i=1:nnod-1)
          ter2(j,i) = 0.5*q_half(j,i)*(g_i(j,i)+g_i(j,i+1))
          ter3(j,i) = q2_half(j,i)*alfa_half(j,i)
       end forall
    end do

    ter23 = 0.0d0
    do k=1,3
       do j=1,3
          forall(i = 1:nnod-1)
             ter23(k,i) = ter23(k,i)+p_half(k,j,i)*(ter2(j,i)-ter3(j,i))
          end forall
       end do
    end do

    ter4 = 0.0d0
    do k=1,3
       do j=1,3
          forall(i = 1:nnod-1)
             ter4(k,i) = ter4(k,i)+j_half(k,j,i)*c_half(j,i)
          end forall
       end do
    end do

    do j=1,3
       forall(i = 1:nnod-1)
          f_tvd(j,i) = 0.5*( fa(j,i)+fa(j,i+1) +  &
               1.0d0/r(i)*ter23(j,i) + dt*ter4(j,i) )
       end forall
    end do

  end subroutine fluxTVD

  subroutine diagsys(U, ga, n, P, P1, D)
    !
    !  Returns the diagonal system for the flux equations
    !
    !    p = [p11 p21 ... pn1 
    !         p12 p22 ... pn2  
    !		   ... 
    !         p1n p2n ... pnn]
    !
    !  diagsys is called by: fluxTVD, bc_abso, bc_abso_unsteady
    !  diagsys calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: n
    real*8, intent(in) :: ga
    real*8, dimension(3,n), intent(in) :: U
    real*8, dimension(3,n), intent(out) :: D
    real*8, dimension(9,n), intent(out) :: P,P1

    real*8 :: tol,g1
    real*8, dimension(n) :: ucell,pres,cc,vaux

    tol   = 1.0d-16
    g1    =  ga-1.

    ucell =  U(2,:)/U(1,:)
    ! pres  =  (U(3,:)-0.5*U(1,:)*U(2,:))*g1 ! CHEQUEAR
    pres  =  (U(3,:)-0.5*U(2,:)**2/U(1,:))*g1
    cc    =  dsqrt(ga*pres/U(1,:))

    vaux = g1/2./cc**3

    P1(1,:) = vaux*(cc**2/g1+0.5*ucell*cc)*ucell
    P1(2,:) = vaux*(2*cc**2/g1-ucell**2)*cc
    P1(3,:) = vaux*ucell*cc*(0.5*ucell-cc/g1)
    P1(4,:) = vaux*(-(cc**2/g1+ucell*cc))
    P1(5,:) = vaux*2*ucell*cc
    P1(6,:) = vaux*(cc**2/g1-ucell*cc)
    P1(7,:) = vaux*cc
    P1(8,:) = vaux*(-2.*cc)
    P1(9,:) = vaux*cc

    P(1,:) = 1.
    P(2,:) = ucell-cc
    P(3,:) = cc**2/g1+0.5*ucell**2-ucell*cc
    P(4,:) = 1.
    P(5,:) = ucell
    P(6,:) = 0.5*ucell**2
    P(7,:) = 1.
    P(8,:) = ucell+cc
    P(9,:) = cc**2/g1+0.5*ucell**2+ucell*cc

    D(1,:) = ucell-cc
    D(2,:) = ucell
    D(3,:) = ucell+cc

  end subroutine diagsys

  subroutine qtvd2(d, b, n, q)
    !
    !
    !
    !  qtvd2 is called by: fluxTVD
    !  qtvd2 calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: n
    real*8, dimension(3,n), intent(in) :: d,b
    real*8, dimension(3,n), intent(out) :: q

    integer itype
    real*8 :: epsil
    
    itype = 1
    epsil = 0.1d0

    if(itype.eq.0) then
       where (dabs(d) < 2.0d0*epsil)
          q = d**2+epsil+dabs(b)
       elsewhere
          q = dabs(d)+dabs(b)
       end where
    else
       where (dabs(d+b) < epsil)
          q = ((d+b)**2+epsil**2)/(2.0d0*epsil)
       elsewhere
          q = dabs(d+b)
       end where
    endif

  end subroutine qtvd2

  subroutine ajac(U, ga, n, A)
    !
    !   Computes the advective jacobian matrix in conservative variables
    !
    !  ajac is called by: fluxTVD
    !  ajac calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: n
    real*8, intent(in) :: ga
    real*8, dimension(3,n), intent(in) :: U
    real*8, dimension(3,3,n), intent(out) :: A

    real*8, dimension(n) :: ucell,ecell
    real*8 :: g1

    ucell = U(2,:)/U(1,:)
    ecell = U(3,:)/U(1,:)
    g1    = ga-1.

    A = 0.0
    A(1,2,:) = 1.0
    A(2,1,:) = 0.5*(ga-3.0)*ucell**2
    A(2,2,:) = -(ga-3.0)*ucell
    A(2,3,:) = g1
    A(3,1,:) = (g1*ucell**2-ga*ecell)*ucell
    A(3,2,:) = ga*ecell-1.5*g1*ucell**2
    A(3,3,:) = ga*ucell

  end subroutine ajac

  subroutine bc_tubes(U, Uref, Uin, ga, normal, dt, dx)
    !
    !  WARNING! Saco las bc Dirichlet, todas son abso
    !
    !  bc_tubes is called by: solve_tube
    !  bc_tubes calls the following subroutines and functions: bc_abso_unsteady, bc_abso
    !
    implicit none

    integer, intent(in) :: normal
    real*8, intent(in) :: ga,dt,dx
    real*8, dimension(3,1), intent(in) :: Uref,Uin
    real*8, dimension(3,1), intent(inout) :: U

    real*8 :: alfa
    real*8, dimension(3,1) :: Uref_C,Ubc
    logical :: unsteady_absorbent

    unsteady_absorbent = .false.

    call pv2cv(Uref, Uref_C, ga, 1)
    
    if(unsteady_absorbent) then
       call bc_abso_unsteady(U, normal, Uin, Uref_C, ga, Ubc, dt, dx)
    else
       call bc_abso(U, normal, Uin, Uref_C, ga, Ubc)
    endif

    alfa = 0.
    U = alfa*U+(1.-alfa)*Ubc

  end subroutine bc_tubes

  subroutine bc_abso(U, normal, Uin, Uref, ga, U_bc)
    !
    !
    !
    !  bc_abso is called by: bc_tubes
    !  bc_abso calls the following subroutines and functions: diagsys
    !
    implicit none

    integer, intent(in) :: normal
    real*8, intent(in) :: ga
    real*8, dimension(3,1), intent(in) :: U,Uin,Uref
    real*8, dimension(3,1), intent(out) :: U_bc

    integer :: j
    real*8, dimension(3,1) :: D,vin,vref,v
    real*8, dimension(9,1) :: P,P1

    call diagsys(U, ga, 1, P, P1, D)

    forall(j = 1:3)
       D(j,:) = D(j,:)*normal
    end forall

    vin(1,:) = P1(1,:)*Uin(1,:)+P1(4,:)*Uin(2,:)+P1(7,:)*Uin(3,:)
    vin(2,:) = P1(2,:)*Uin(1,:)+P1(5,:)*Uin(2,:)+P1(8,:)*Uin(3,:)
    vin(3,:) = P1(3,:)*Uin(1,:)+P1(6,:)*Uin(2,:)+P1(9,:)*Uin(3,:)

    vref(1,:) = P1(1,:)*Uref(1,:)+P1(4,:)*Uref(2,:)+P1(7,:)*Uref(3,:)
    vref(2,:) = P1(2,:)*Uref(1,:)+P1(5,:)*Uref(2,:)+P1(8,:)*Uref(3,:)
    vref(3,:) = P1(3,:)*Uref(1,:)+P1(6,:)*Uref(2,:)+P1(9,:)*Uref(3,:)

    where(D.le.0.0) v = vref
    where(D.gt.0.0) v = vin

    U_bc(1,:) = P(1,:)*v(1,:)+P(4,:)*v(2,:)+P(7,:)*v(3,:)
    U_bc(2,:) = P(2,:)*v(1,:)+P(5,:)*v(2,:)+P(8,:)*v(3,:)
    U_bc(3,:) = P(3,:)*v(1,:)+P(6,:)*v(2,:)+P(9,:)*v(3,:)

  end subroutine bc_abso

  subroutine bc_abso_unsteady(U, normal, Uin, Uref, ga, U_bc, dt, dx)
    !
    !
    !
    !  bc_abso_unsteady is called by: bc_tubes
    !  bc_abso_unsteady calls the following subroutines and functions: diagsys
    !
    implicit none

    integer, intent(in) :: normal
    real*8, intent(in) :: ga,dt,dx
    real*8, dimension(3,1), intent(in) :: U,Uin,Uref
    real*8, dimension(3,1), intent(out) :: U_bc

    integer :: j
    real*8, dimension(3,1) :: D,vin,vref,v,dvdx
    real*8, dimension(9,1) :: P,P1

    call diagsys(U, ga, 1, P, P1, D)

    forall(j = 1:3)
       D(j,:) = D(j,:)*normal
    end forall

    v(1,:) = P1(1,:)*U(1,:)+P1(4,:)*U(2,:)+P1(7,:)*U(3,:)
    v(2,:) = P1(2,:)*U(1,:)+P1(5,:)*U(2,:)+P1(8,:)*U(3,:)
    v(3,:) = P1(3,:)*U(1,:)+P1(6,:)*U(2,:)+P1(9,:)*U(3,:)

    vin(1,:) = P1(1,:)*Uin(1,:)+P1(4,:)*Uin(2,:)+P1(7,:)*Uin(3,:)
    vin(2,:) = P1(2,:)*Uin(1,:)+P1(5,:)*Uin(2,:)+P1(8,:)*Uin(3,:)
    vin(3,:) = P1(3,:)*Uin(1,:)+P1(6,:)*Uin(2,:)+P1(9,:)*Uin(3,:)

    vref(1,:) = P1(1,:)*Uref(1,:)+P1(4,:)*Uref(2,:)+P1(7,:)*Uref(3,:)
    vref(2,:) = P1(2,:)*Uref(1,:)+P1(5,:)*Uref(2,:)+P1(8,:)*Uref(3,:)
    vref(3,:) = P1(3,:)*Uref(1,:)+P1(6,:)*Uref(2,:)+P1(9,:)*Uref(3,:)

    !where(D.le.0.0) dvdx = (vref-vin)/dx
    where(D.le.0.0) dvdx = (vref-v)/dx
    where(D.gt.0.0) dvdx = (v-vin)/dx

    v = v - dt*D*dvdx

    U_bc(1,:) = P(1,:)*v(1,:)+P(4,:)*v(2,:)+P(7,:)*v(3,:)
    U_bc(2,:) = P(2,:)*v(1,:)+P(5,:)*v(2,:)+P(8,:)*v(3,:)
    U_bc(3,:) = P(3,:)*v(1,:)+P(6,:)*v(2,:)+P(9,:)*v(3,:)

  end subroutine bc_abso_unsteady

  subroutine critical_dt(globalData, Up, hnod, nnod, dt)

    implicit none

    integer, intent(in) :: nnod
    real*8, dimension(nnod-1), intent(in) :: hnod
    real*8, dimension(3,nnod), intent(in) :: Up
    type(dataSim), intent(in) :: globalData
    real*8, intent(out) :: dt

    real*8 :: dtmin,dt_rpm,hmin,velmax,sonic_speed_min,ratio_rpm
    real*8, dimension(nnod) :: sonic_speed,vel_star

    if(globalData%rpm.gt.0.0d0)then
       dt_rpm    = globalData%dtheta_rpm/(6.*globalData%rpm)
       ratio_rpm = globalData%rpm_ini/globalData%rpm
    else
       dt_rpm    = 1e+10
       ratio_rpm = 1.
    endif

    sonic_speed = globalData%ga*Up(3,:)/Up(1,:)

    sonic_speed_min = minval(sonic_speed)

    if(sonic_speed_min.le.0.0d0) then
       write(6,*)' Negative sonic speed @ critical_dt , min(c)=', &
            sonic_speed_min
    endif

    vel_star = dsqrt(abs(sonic_speed))+dabs(Up(2,:))

    hmin   = minval(hnod)
    velmax = maxval(vel_star)
    dtmin  = globalData%Courant*ratio_rpm*hmin/velmax
    dt     = dmin1(dtmin,dt_rpm)

 end subroutine critical_dt

 subroutine actualize_valves(itube, Area, Twall, dAreax, nnod)
   !
   !
   !
    use def_cylinder

    implicit none

    integer, intent(in) :: itube, nnod
    real*8, dimension(nnod), intent(in) :: dAreax,Area,Twall

    integer :: i,k

    if(size(cyl).gt.0) then
       do k=1,size(cyl),1
          do i=1,cyl(k)%nvi,1
             if(itube == cyl(k)%intake_valves(i)%tube) then
                cyl(k)%intake_valves(i)%Area_tube   = Area(nnod)
                cyl(k)%intake_valves(i)%dAreax_tube = dAreax(nnod)
                cyl(k)%intake_valves(i)%twall_tube  = Twall(nnod)
             endif
          enddo
          do i=1,cyl(k)%nve,1
             if(itube == cyl(k)%exhaust_valves(i)%tube) then
                cyl(k)%exhaust_valves(i)%Area_tube   = Area(1)
                cyl(k)%exhaust_valves(i)%dAreax_tube = dAreax(1)
                cyl(k)%exhaust_valves(i)%twall_tube  = Twall(1)
             endif
          enddo
       enddo
    endif

  end subroutine actualize_valves

end module def_tube
