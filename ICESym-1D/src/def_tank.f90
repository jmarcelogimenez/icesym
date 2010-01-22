module def_tank
     
  use def_simulator
  use def_valve
  use gasdyn_utils
  use, intrinsic :: ISO_C_BINDING
  
  type, BIND(C) :: this
    integer(C_INT) :: nnod, ndof,nnod_input,nunit
    real(C_DOUBLE) :: Volume, mass, h_film, Area_wall, T_wall
 end type this

contains
  
  subroutine state_initial_tank(myData, atm, globalData, state_ini, type_end, Cd_ports)BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(this) :: myData
    real(C_DOUBLE) :: state_ini(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    integer(C_INT), dimension(0:myData%nnod_input-1) :: type_end
    real(C_DOUBLE), dimension(0:myData%nnod_input-1) :: Cd_ports
    real(C_DOUBLE), dimension(0:2) :: atm
    type(dataSim) :: globalData

    integer :: i,j
    real*8 :: rho0,p0,T0

    rho0 = atm(0)
    p0   = atm(2)
    T0   = p0/(globalData%R_gas*rho0)

    state_ini(0) = rho0
    state_ini(1) = p0
    state_ini(2) = T0
    do i=1,myData%nnod-1
       do j=0,myData%ndof-1
          state_ini(i*myData%ndof+j) = atm(j)
       end do
    end do

    myData%mass = rho0*myData%Volume

  end subroutine state_initial_tank
  
  subroutine solve_tank(myData, globalData, state, new_state, &
       type_end, Area_tube,  twall_tube, dAreax_tube, Cd_ports) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(this) :: myData
    type(dataSim) :: globalData
    real(C_DOUBLE) :: state(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: new_state(0:(myData%nnod*myData%ndof)-1)  
    integer(C_INT), dimension(0:myData%nnod_input-1) :: type_end
    real(C_DOUBLE), dimension(0:myData%nnod_input-1) :: Area_tube, twall_tube, &
         dAreax_tube, Cd_ports
    
    integer :: ntubes,i,j,itube
    real*8 :: Area_P,Area_T,alpha,dAreax_P,Twall_P
    real*8, dimension(myData%ndof) :: Utank,Utube,Uthroat,RHS
    real*8, dimension(myData%ndof,myData%nnod-1) :: Upipe,Uref,Utpipe
    
    ntubes = myData%nnod-1

    ! Get the tank state from state vector
    do j=0,myData%ndof-1
       Utank(j+1) = state(j)
    enddo
    ! Get the tubes state(s) from state vector
    do i=1,myData%nnod-1
       do j=0,myData%ndof-1
          Upipe(j+1,i) = state(i*myData%ndof+j)
       enddo
    enddo
    ! Get the reference state(s) from state vector
    do i=myData%nnod,myData%nnod+myData%nnod_input-1
       do j=0,myData%ndof-1
          Uref(j+1,i-myData%nnod+1) = state(i*myData%ndof+j)
       enddo
    enddo
    
    do itube=1,ntubes
       Area_P   = Area_tube(itube-1)
       Area_T   = Cd_ports(itube-1)*Area_P
       dAreax_P = dAreax_tube(itube-1)
       Twall_P  = twall_tube(itube-1)
       alpha = 1.
       if(.true.) then
          call solve_valve(globalData, Utank, Uref(:,itube), &
               type_end(itube-1), Area_T, Area_P, Utube, Uthroat)
          alpha = 0.15
       else
          call rhschar(Uref(:,itube), Area_P, dAreax_P, Twall_P, globalData%ga, &
               globalData%R_gas, globalData%dt, globalData%viscous_flow, &
               globalData%heat_flow, RHS)
          Utube = Upipe(:,itube)
          call solve_valve_implicit(Utank, Uref(:,itube), type_end(itube-1), &
               Area_T, Area_P, RHS, globalData%ga, Utube, Uthroat)
       end if
       Upipe(:,itube)  = alpha*Utube + &
            (1.-alpha)*Upipe(:,itube)
       Utpipe(:,itube) = Uthroat
    end do

    call tank_solver(myData, globalData, Utank, Upipe, Utpipe, &
       Area_tube, Cd_ports, type_end)

    ! Actualizes the tank state
    do j=0,myData%ndof-1
       new_state(j) = Utank(j+1)
    enddo
    ! Actualizes the state(s) at the pipes
    do i=1,myData%nnod-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = Upipe(j+1,i)
       enddo
    enddo

  end subroutine solve_tank

  subroutine tank_solver(myData, globalData, Utank, Upipes, Utpipes, &
       Area_P, Cd_ports, type_end)
    !
    !
    !
    !  tank_solver is called by:
    !  tank_solver calls the following subroutines and functions:
    !
    implicit none

    type(this), intent(inout) :: myData
    integer, dimension(myData%nnod-1), intent(in) :: type_end
    real*8, dimension(myData%nnod-1), intent(in) :: Area_P,Cd_ports
    real*8, dimension(3,myData%nnod-1), intent(in) :: Upipes,Utpipes
    type(dataSim), intent(in) :: globalData 
    real*8, dimension(3), intent(inout) :: Utank

    integer :: ispecie,ntubes
    real*8 :: dt,cp,cv,rho_tank,p_tank,T_tank
    real*8 :: mass_old,mass_new,tol
    real*8 :: dQ_ht,edot,ene_old,ene_new
    real*8, dimension(myData%nnod-1) :: mdots,hdots

    tol = 1.0d-6

    dt = globalData%dt

    ntubes = myData%nnod-1

    rho_tank = Utank(1)
    p_tank   = Utank(2)
    T_tank   = Utank(3)

    ispecie = 0
    cp = compute_cp(T_tank, ispecie, globalData%R_gas)
    cv = cp - globalData%R_gas
    
    call flow_rates(ntubes, type_end, globalData%ga, Area_P, Cd_ports, &
         Upipes, Utpipes, Utank, mdots, hdots)
    
    mass_old = myData%mass
    mass_new = dt*sum(mdots) + mass_old
    mass_new = dmax1(mass_new, tol)

    dQ_ht = myData%h_film*myData%Area_wall*(T_tank-myData%T_wall)

    edot    = sum(hdots)-dQ_ht
    ene_old = mass_old*cv*T_tank
    ene_new = ene_old + dt*edot
    ene_new = dmax1(ene_new, tol)
    
    ! Continuity equation
    rho_tank = mass_new/myData%Volume
    ! Energy equation
    T_tank  = ene_new/(mass_new*cv)
    ! State equation
    p_tank  = rho_tank*globalData%R_gas*T_tank
    
    myData%mass = mass_new

    Utank(1) = rho_tank
    Utank(2) = p_tank
    Utank(3) = T_tank
    !if(globalData%save_extras) then
    !   write(myData%nunit,*) globalData%icycle, globalData%theta, globalData%time, mdots,hdots,mass_new,dQ_ht
    !endif
  end subroutine tank_solver

  subroutine flow_rates(ntubes, type_end, ga, Area_P, Cd_ports, Upipe, Utpipe, &
       Utank, mdots, hdots)
    !
    !
    !  flow_rates is called by: tank_solver
    !  flow_rates calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: ntubes
    integer, dimension(ntubes), intent(in) :: type_end
    real*8, intent(in) :: ga
    real*8, dimension(3), intent(in) :: Utank
    real*8, dimension(ntubes), intent(in) :: Area_P,Cd_ports
    real*8, dimension(3,ntubes), intent(in) :: Upipe,Utpipe
    real*8, dimension(ntubes), intent(out) :: mdots,hdots

    integer :: i
    real*8 :: g1,rho,vel,pre

    g1 = ga-1.

    do i=1,ntubes
       rho = Utpipe(1,i)
       vel = Utpipe(2,i)
       pre = Utpipe(3,i)

       mdots(i) = type_end(i)*rho*vel*Area_P(i)*Cd_ports(i)
       
       rho = Upipe(1,i)
       vel = Upipe(2,i)
       pre = Upipe(3,i)

       if(vel*type_end(i).gt.0.) then
          ! Outflow from end pipe
          hdots(i) = mdots(i)*(ga*pre/rho/g1+0.5*vel**2)
       else
          ! Inflow to end pipe
          hdots(i) = mdots(i)*ga*Utank(2)/Utank(1)/g1
       endif
    end do

  end subroutine flow_rates

end module def_tank
