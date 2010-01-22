module def_cylinder
  
  use def_simulator
  use def_valve
  use gasdyn_utils
  use utilities
  use, intrinsic :: ISO_C_BINDING
  
  type, BIND(C) :: fuel
     real(C_DOUBLE) ::Q_fuel, y, hvap_fuel
  end type fuel

  type :: injection
     integer :: pulse, ignition_delay_model
     real*8 :: m_inj, dtheta_inj, T_fuel, theta_inj_ini, theta_id, integral
     real*8, dimension(:,:), pointer :: mfdot_array
  end type injection

  type :: combustion
     integer :: combustion_model
     real*8, dimension(:,:),pointer  :: xbdot_array
     real*8 :: theta_ig_0, dtheta_comb, phi, phi_ig, a_wiebe, m_wiebe
     real*8 :: mass_fuel_ini, mass_air_ini, mass_res_ini
     logical :: start_comb
  end type combustion

  type, BIND(C) :: scavenge
     real(C_DOUBLE) :: val_1, val_2, SRv
     logical(C_BOOL) :: close_cyl
  end type scavenge

  type, BIND(C) :: this
     integer(C_INT) :: nnod_input,nvi,nve,nnod,ndof,model_ht, type_ig, nunit,species_model,ntemp
     real(C_DOUBLE) :: Bore,crank_radius,Vol_clearance,rod_length,head_chamber_area, &
          piston_area,theta_0,delta_ca,Twall,factor_ht
     logical(C_BOOL) :: scavenge, full_implicit,nh_temp
  end type this

  type :: cylinder
     type(fuel)::fuel_data
     type(injection)::injection_data
     type(combustion)::combustion_data
     type(scavenge)::scavenge_data
     type(valve),dimension(:),pointer::intake_valves
     type(valve),dimension(:),pointer::exhaust_valves
     real*8,dimension(:),pointer :: prop, U_crevice, data_crevice
     integer :: nvi, nve
  end type cylinder
  
  type (cylinder),dimension(:),allocatable :: cyl
  integer :: numcyl

contains

  subroutine initialize_cylinders(ncyl) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: ncyl
    allocate(cyl(ncyl))
    numcyl = ncyl
    return
  end subroutine initialize_cylinders

  subroutine initialize_valves(icyl,nvi,nve) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl,nvi,nve
    allocate(cyl(icyl)%intake_valves(nvi))
    allocate(cyl(icyl)%exhaust_valves(nve))
    cyl(icyl)%nvi = nvi
    cyl(icyl)%nve = nve
    return
  end subroutine initialize_valves

  subroutine initialize_arrays(icyl, prop, U_crevice, data_crevice,l1,l2,l3) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl,l1,l2,l3
    real(C_DOUBLE) :: prop(0:l1-1), U_crevice(0:l2-1), data_crevice(0:l3-1)
    allocate(cyl(icyl)%prop(l1))
    allocate(cyl(icyl)%U_crevice(l2))
    allocate(cyl(icyl)%data_crevice(l3))
    
    do i=1,(l1)
       cyl(icyl)%prop(i) = prop(i-1)
    enddo
    do i=1,(l2)
       cyl(icyl)%U_crevice(i) = U_crevice(i-1)
    enddo
    do i=1,(l3)
       cyl(icyl)%data_crevice(i) = data_crevice(i-1)
    enddo
    return
  end subroutine initialize_arrays
  
  subroutine initialize_fuel(icyl,fuelData) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl
    type(fuel)::fuelData
    cyl(icyl)%fuel_data%Q_fuel    = fuelData%Q_fuel
    cyl(icyl)%fuel_data%y         = fuelData%y
    cyl(icyl)%fuel_data%hvap_fuel = fuelData%hvap_fuel
  end subroutine initialize_fuel
  
  subroutine initialize_scavenge(icyl,scavengeData) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl
    type(scavenge)::scavengeData
    cyl(icyl)%scavenge_data%val_1     = scavengeData%val_1
    cyl(icyl)%scavenge_data%val_2     = scavengeData%val_2
    cyl(icyl)%scavenge_data%SRv       = scavengeData%SRv
    cyl(icyl)%scavenge_data%close_cyl = scavengeData%close_cyl
    return
  end subroutine initialize_scavenge
  
  subroutine initialize_injection(icyl,pulse,m_inj,dtheta_inj,T_fuel,theta_inj_ini, &
       theta_id, integral, mfdot_array, ignition_delay_model, l1) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl, pulse, l1,ignition_delay_model
    real(C_DOUBLE) :: m_inj,dtheta_inj,T_fuel,theta_inj_ini, theta_id, integral
    real(C_DOUBLE) :: mfdot_array(0:l1-1)
    cyl(icyl)%injection_data%pulse                = pulse
    cyl(icyl)%injection_data%m_inj                = m_inj
    cyl(icyl)%injection_data%dtheta_inj           = dtheta_inj
    cyl(icyl)%injection_data%T_fuel               = T_fuel
    cyl(icyl)%injection_data%theta_inj_ini        = theta_inj_ini
    cyl(icyl)%injection_data%theta_id             = theta_id
    cyl(icyl)%injection_data%integral             = integral
    cyl(icyl)%injection_data%ignition_delay_model = ignition_delay_model

    allocate(cyl(icyl)%injection_data%mfdot_array(2,l1/2))
    do i=1,(l1/2)
       cyl(icyl)%injection_data%mfdot_array(i,1) = mfdot_array((i-1)*2)
       cyl(icyl)%injection_data%mfdot_array(i,2) = mfdot_array((i-1)*2 +1)
    enddo
    
    return
  end subroutine initialize_injection
  
  subroutine initialize_combustion(icyl,theta_ig_0, dtheta_comb, phi, phi_ig, a_wiebe, m_wiebe, &
       xbdot_array, combustion_model, start_comb, l1) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl,l1,combustion_model
    real(C_DOUBLE) :: xbdot_array(0:(l1-1))
    real(C_DOUBLE) :: theta_ig_0, dtheta_comb, phi, phi_ig, a_wiebe, m_wiebe
    logical(C_BOOL) :: start_comb

    cyl(icyl)%combustion_data%theta_ig_0       = theta_ig_0
    cyl(icyl)%combustion_data%dtheta_comb      = dtheta_comb
    cyl(icyl)%combustion_data%phi              = phi
    cyl(icyl)%combustion_data%phi_ig           = phi_ig
    cyl(icyl)%combustion_data%a_wiebe          = a_wiebe
    cyl(icyl)%combustion_data%m_wiebe          = m_wiebe
    cyl(icyl)%combustion_data%combustion_model = combustion_model
    cyl(icyl)%combustion_data%start_comb       = start_comb

    allocate(cyl(icyl)%combustion_data%xbdot_array(2,l1/2))
    do i=1,(l1/2)
       cyl(icyl)%combustion_data%xbdot_array(i,1) = xbdot_array((i-1)*2)
       cyl(icyl)%combustion_data%xbdot_array(i,2) = xbdot_array((i-1)*2 +1)
    enddo
    
  end subroutine initialize_combustion
  
  subroutine initialize_intake_valves(icyl, ival, Nval, type_dat, angle_VO, angle_VC, &
       Dv, Lvmax, Cd, Lv, valve_model, l1, l2, dx_tube, Area_tube, twall_tube, &
       dAreax_tube, tube) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl,ival,l1,l2,tube
    integer(C_INT) :: Nval, type_dat,valve_model
    real(C_DOUBLE) :: Cd(0:l1-1), Lv(0:l2-1)
    real(C_DOUBLE) :: angle_VO, angle_VC, Dv, Lvmax
    real(C_DOUBLE) :: Area_tube, twall_tube, dAreax_tube, dx_tube
    
    cyl(icyl)%intake_valves(ival)%tube        = tube
    cyl(icyl)%intake_valves(ival)%Nval        = Nval
    cyl(icyl)%intake_valves(ival)%type_dat    = type_dat
    cyl(icyl)%intake_valves(ival)%angle_VO    = angle_VO
    cyl(icyl)%intake_valves(ival)%angle_VO    = angle_VO
    cyl(icyl)%intake_valves(ival)%angle_VC    = angle_VC
    cyl(icyl)%intake_valves(ival)%dx_tube     = dx_tube
    cyl(icyl)%intake_valves(ival)%Area_tube   = Area_tube
    cyl(icyl)%intake_valves(ival)%dAreax_tube = dAreax_tube
    cyl(icyl)%intake_valves(ival)%twall_tube  = twall_tube
    cyl(icyl)%intake_valves(ival)%Dv          = Dv
    cyl(icyl)%intake_valves(ival)%Lvmax       = Lvmax
    cyl(icyl)%intake_valves(ival)%valve_model = valve_model

    allocate(cyl(icyl)%intake_valves(ival)%Cd(l1/2,2))
    allocate(cyl(icyl)%intake_valves(ival)%Lv(l2/2,2))
    do i=1,(l1/2)
       cyl(icyl)%intake_valves(ival)%Cd(i,1) = Cd((i-1)*2)
       cyl(icyl)%intake_valves(ival)%Cd(i,2) = Cd((i-1)*2 +1)
    enddo
    do i=1,(l2/2)
       cyl(icyl)%intake_valves(ival)%Lv(i,1) = Lv((i-1)*2)
       cyl(icyl)%intake_valves(ival)%Lv(i,2) = Lv((i-1)*2 +1)
    enddo

  end subroutine initialize_intake_valves
  
  subroutine initialize_exhaust_valves(icyl, ival, Nval, type_dat, angle_VO, angle_VC, &
       Dv, Lvmax, Cd, Lv, valve_model, l1, l2, dx_tube, Area_tube, twall_tube, dAreax_tube, &
       tube) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT) :: icyl,ival,l1,l2,tube
    integer(C_INT) :: Nval, type_dat,valve_model
    real(C_DOUBLE) :: Cd(0:l1-1), Lv(0:l2-1)
    real(C_DOUBLE) :: angle_VO, angle_VC, Dv, Lvmax
    real(C_DOUBLE) :: Area_tube, twall_tube, dAreax_tube, dx_tube

    cyl(icyl)%exhaust_valves(ival)%tube        = tube
    cyl(icyl)%exhaust_valves(ival)%Nval        = Nval
    cyl(icyl)%exhaust_valves(ival)%type_dat    = type_dat
    cyl(icyl)%exhaust_valves(ival)%angle_VO    = angle_VO
    cyl(icyl)%exhaust_valves(ival)%angle_VC    = angle_VC
    cyl(icyl)%exhaust_valves(ival)%dx_tube     = dx_tube
    cyl(icyl)%exhaust_valves(ival)%Area_tube   = Area_tube
    cyl(icyl)%exhaust_valves(ival)%dAreax_tube = dAreax_tube
    cyl(icyl)%exhaust_valves(ival)%twall_tube  = twall_tube
    cyl(icyl)%exhaust_valves(ival)%Dv          = Dv
    cyl(icyl)%exhaust_valves(ival)%Lvmax       = Lvmax
    cyl(icyl)%exhaust_valves(ival)%valve_model = valve_model
       
    allocate(cyl(icyl)%exhaust_valves(ival)%Cd(l1/2,2))
    allocate(cyl(icyl)%exhaust_valves(ival)%Lv(l2/2,2))
    do i=1,(l1/2)
       cyl(icyl)%exhaust_valves(ival)%Cd(i,1) = Cd((i-1)*2)
       cyl(icyl)%exhaust_valves(ival)%Cd(i,2) = Cd((i-1)*2 +1)
    enddo
    do i=1,(l2/2)
       cyl(icyl)%exhaust_valves(ival)%Lv(i,1) = Lv((i-1)*2)
       cyl(icyl)%exhaust_valves(ival)%Lv(i,2) = Lv((i-1)*2 +1)
    enddo
    
  end subroutine initialize_exhaust_valves

  subroutine state_initial_cylinder(icyl, myData, atm, globalData, state_ini, mass_C, twall) BIND(C)
    !
    !
    !  state_initial_cylinder is called by:
    !  state_initial_cylinder calls the following subroutines and functions:
    !    geometry, run_ideal_cycle, interpolant, fa_ratio
    !
    use, intrinsic :: ISO_C_BINDING
    type(this) :: myData
    type(dataSim) :: globalData
    real(C_DOUBLE) :: atm(0:2)
    integer(C_INT) :: icyl,ntemp
    real(C_DOUBLE) :: state_ini(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: mass_C(6*(myData%nnod - myData%nvi - myData%nve ) )
    real(C_DOUBLE) :: twall(myData%ntemp)


    integer :: i
    real*8 :: theta,Vol,Area,Vdot,phi,y,factor,theta_ig_0,theta_ig_f,EVC
    real*8 :: mass_fuel_ini,mass_air_ini,mass_res_ini,alpha,rho,pres,temp
    real*8, dimension(3) :: atm_state
    real*8, dimension(7,4) :: avpt
    do i=0,2
       atm_state(i+1) = atm(i)
    end do

    theta = myData%theta_0*180/pi

    call geometry(myData, globalData, Vol, Area, Vdot)
    call run_ideal_cycle(icyl, myData, globalData, atm_state, avpt)
    ! The initial state is intepolated from an ideal cycle
    call interpolant(avpt(:,1), avpt(:,2), theta, rho)
    call interpolant(avpt(:,1), avpt(:,3), theta, pres)
    call interpolant(avpt(:,1), avpt(:,4), theta, temp)

    mass_air_ini = Vol*rho
    if(myData%type_ig.eq.0) then
       phi           = cyl(icyl)%combustion_data%phi
       y             = cyl(icyl)%fuel_data%y
       factor        = fa_ratio(phi,y)
       mass_fuel_ini = mass_air_ini*factor

       theta_ig_0 = cyl(icyl)%combustion_data%theta_ig_0*180./pi
       theta_ig_f = theta_ig_0+cyl(icyl)%combustion_data%dtheta_comb*180./pi
    else
       theta_ig_0 = (cyl(icyl)%injection_data%theta_inj_ini + &
            cyl(icyl)%injection_data%theta_id)*180./pi
       theta_ig_f = theta_ig_0 + cyl(icyl)%combustion_data%dtheta_comb*180./pi

       mass_fuel_ini = cyl(icyl)%injection_data%m_inj
       y             = cyl(icyl)%fuel_data%y
       factor        = fa_ratio(1.0d0,y)
       cyl(icyl)%combustion_data%phi_ig = mass_fuel_ini/(mass_air_ini*factor)
    end if
    
    EVC = cyl(icyl)%exhaust_valves(1)%angle_VC*180./pi
    
    alpha = 0.0
    if(globalData%nstroke.eq.4) then
       if((theta.gt.theta_ig_0).and.(theta.le.theta_ig_f)) then
          ! Combustion
          alpha = (theta-theta_ig_0)/(theta_ig_f-theta_ig_0)
       elseif((theta.gt.theta_ig_f).and.(theta.le.720.)) then
          ! Expansion, gas exhaust or partial valve overlap
          alpha = 1.0
       elseif((theta.ge.0).and.(theta.le.EVC))then
          ! Valve overlap with partial load
          alpha = 1.0-(theta-0.)/(EVC-0.)
       endif
    elseif(globalData%nstroke.eq.2) then
       if((theta.gt.theta_ig_0).or.(theta.le.theta_ig_f)) then
          ! Combustion
          alpha = (theta-theta_ig_0-floor((theta-theta_ig_0)/360.)*360.)/ &
               (theta_ig_f-theta_ig_0-floor((theta_ig_f-theta_ig_0)/360.)*360.)
       elseif((theta.gt.theta_ig_f).and.(theta.le.180.)) then
          ! Expansion, gas exhaust or partial valve overlap
          alpha = 1.
       elseif((theta.ge.180.).and.(theta.le.EVC)) then
          ! Valve overlap with partial load
          alpha = 1.0-(theta-180.)/(EVC-180.)
       endif
    endif
    mass_fuel_ini = (1.-alpha)*mass_fuel_ini
    mass_res_ini  = alpha*mass_air_ini
    mass_air_ini  = (1.-alpha)*mass_air_ini
    ! Cylinder state
    state_ini(0) = rho
    state_ini(1) = pres
    state_ini(2) = temp
    ! Valves states
    do i=1,myData%nnod-1
       state_ini(i*myData%ndof+0) = atm_state(1)
       state_ini(i*myData%ndof+1) = atm_state(2)
       state_ini(i*myData%ndof+2) = atm_state(3)
    enddo

    do i=1,myData%nnod-myData%nvi-myData%nve
       mass_C(6*(i-1) + 1) = mass_fuel_ini
       mass_C(6*(i-1) + 2) = mass_air_ini
       mass_C(6*(i-1) + 3) = mass_res_ini
       mass_C(6*(i-1) + 4) = mass_fuel_ini
       mass_C(6*(i-1) + 5) = mass_air_ini
       mass_C(6*(i-1) + 6) = mass_res_ini
    end do

  end subroutine state_initial_cylinder

  subroutine solve_cylinder(icyl,myData, globalData, state, new_state, mass_C, twall) BIND(C)
    !
    !
    !  solve_cylinder is called by: 
    !  solve_cylinder calls the following subroutines and functions: 
    !    area_valve, solve_valve, cylinder_solver
    !
    use, intrinsic :: ISO_C_BINDING
    type(this) :: myData
    type(dataSim) :: globalData
    real(C_DOUBLE) :: state(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: new_state(0:(myData%nnod*myData%ndof)-1)
    real(C_DOUBLE) :: mass_C(6*(myData%nnod - myData%nvi - myData%nve ) )
    real(C_DOUBLE) :: twall(myData%ntemp)
    integer(C_INT) :: icyl

    integer :: nnod_cyl,Nval,type_dat,nstroke,i,j,ival
    real*8 :: theta_g,theta,R_gas
    real*8 :: ga,dx,Area_P,dAreax_P,Twall_P
    real*8 :: Dv,Lvmax,IVO,IVC,EVO,EVC,Area_T
    real*8 :: EGR,alpha
    real*8, dimension(myData%ndof) :: Upipe,Uthroat,RHS
    real*8, dimension(cyl(icyl)%nvi) :: Fiv
    real*8, dimension(cyl(icyl)%nve) :: Fev
    real*8, dimension(3,myData%nnod-myData%nvi-myData%nve) :: mass_cyl,mass_cyl_old
    real*8, dimension(myData%ndof,myData%nnod-myData%nvi-myData%nve) :: Ucyl
    real*8, dimension(myData%ndof,cyl(icyl)%nvi) :: Uriv,Uiv,Utiv,RHSiv
    real*8, dimension(myData%ndof,cyl(icyl)%nve) :: Urev,Uev,Utev,RHSev

! Actualmente Twall se recibe en myData (la primer posicion del arreglo, que en homogeneo ES la twall tradicional)
! Ademas se recibe el arreglo completo por parámetro, y una bandera myData%nh_temp (true=NO homogenea, false=homogenea)
! En Resumen estan los datos para implementar el cylinder con varias temperaturas, y además funciona la versión homogenea tal cual antes

    theta_g = globalData%theta
    nstroke = globalData%nstroke
    ga      = globalData%ga
    R_gas   = globalData%R_gas

    theta = modulo(theta_g+myData%theta_0, &
         nstroke*pi)

    nnod_cyl = myData%nnod-cyl(icyl)%nvi-cyl(icyl)%nve

    EGR = 0.0d0

    if(globalData%save_extras) then
       write(myData%nunit,901) globalData%icycle, theta*180/pi, globalData%time
   901 format (I12,F12.4,E12.4)
    endif

    ! Get the cylinder state(s) from state vector
    do i=0,nnod_cyl-1
       do j=0,myData%ndof-1
          Ucyl(j+1,i+1) = state(i*myData%ndof+j)
       enddo
    enddo
    ! Get the intake valves state(s) from state vector
    do i=nnod_cyl,nnod_cyl+myData%nvi-1
       do j=0,myData%ndof-1
          Uiv(j+1,i-nnod_cyl+1) = state(i*myData%ndof+j)
       enddo
    enddo
    ! Get the exhaust valves state(s) from state vector
    do i=nnod_cyl+myData%nvi,nnod_cyl+myData%nvi+myData%nve-1
       do j=0,myData%ndof-1
          Uev(j+1,i-nnod_cyl-myData%nvi+1) = state(i*myData%ndof+j)
       enddo
    enddo
    ! Get the intake valves reference state(s) from state vector
    do i=myData%nnod,myData%nnod+myData%nvi-1
       do j=0,myData%ndof-1
          Uriv(j+1,i-myData%nnod+1) = state(i*myData%ndof+j)
       enddo
    enddo
    ! Get the exhaust valves state(s) from state vector
    do i=myData%nnod+myData%nvi,myData%nnod+myData%nvi+myData%nve-1
       do j=0,myData%ndof-1
          Urev(j+1,i-myData%nnod-myData%nvi+1) = state(i*myData%ndof+j)
       enddo
    enddo

    do i=1,nnod_cyl
       do j=1,3
          mass_cyl(j,i)     = mass_C(6*(i-1)+j)
          mass_cyl_old(j,i) = mass_C(6*(i-1)+j+3)
       end do
    end do

    ! Intake valves
    do ival=1,cyl(icyl)%nvi
       dx       = cyl(icyl)%intake_valves(ival)%dx_tube
       Area_P   = cyl(icyl)%intake_valves(ival)%Area_tube
       dAreax_P = cyl(icyl)%intake_valves(ival)%dAreax_tube
       Twall_P  = cyl(icyl)%intake_valves(ival)%twall_tube

       Dv       = cyl(icyl)%intake_valves(ival)%Dv
       Lvmax    = cyl(icyl)%intake_valves(ival)%Lvmax
       IVO      = cyl(icyl)%intake_valves(ival)%angle_VO
       IVC      = cyl(icyl)%intake_valves(ival)%angle_VC
       Nval     = cyl(icyl)%intake_valves(ival)%Nval
       type_dat = cyl(icyl)%intake_valves(ival)%type_dat
       call area_valve(theta, Dv, IVO, IVC, Lvmax, Nval, type_dat, &
            cyl(icyl)%intake_valves(ival)%Lv, cyl(icyl)%intake_valves(ival)%Cd, &
            nstroke, Area_T)
       Fiv(ival) = Area_T

        if(cyl(icyl)%intake_valves(ival)%valve_model.eq.1 .or. myData%full_implicit) then
          call rhschar(Uriv(:,ival), Area_P, dAreax_P, Twall_P, ga, R_gas, &
               globalData%dt, globalData%viscous_flow, globalData%heat_flow, RHS)
          RHSiv(:,ival) = RHS
       end if

       if(myData%full_implicit) then
          ! call rhschar(Uriv(:,ival), Area_P, dAreax_P, Twall_P, ga, R_gas, &
          !      globalData%dt, globalData%viscous_flow, globalData%heat_flow, RHS)
          ! RHSiv(:,ival) = RHS
       else
          alpha = 1.          
          if(cyl(icyl)%intake_valves(ival)%valve_model.eq.0) then
             call solve_valve(globalData, Ucyl, Uriv(:,ival), 1, Area_T, Area_P, &
                  Upipe, Uthroat)
             if(Area_T.gt.0.) alpha = 0.15
          else if(cyl(icyl)%intake_valves(ival)%valve_model.eq.1) then
             ! write(*,*) 'INTAKE VALVE - CYL ', icyl
             Upipe = Uiv(:,ival)
             call solve_valve_implicit(Ucyl, Uriv(:,ival), 1, Area_T, Area_P, &
                  RHSiv(:,ival), ga, Upipe, Uthroat)
          else
             stop ' Wrong valve model option '
          end if
          Uiv(:,ival)  = alpha*Upipe + &
               (1.-alpha)*Uiv(:,ival)
          Utiv(:,ival) = Uthroat
       end if
       ! write(10,*) theta, Area_T, Uiv, Utiv
    end do

    ! Exhaust valves
    do ival=1,cyl(icyl)%nve
       dx       = cyl(icyl)%exhaust_valves(ival)%dx_tube
       Area_P   = cyl(icyl)%exhaust_valves(ival)%Area_tube
       dAreax_P = cyl(icyl)%exhaust_valves(ival)%dAreax_tube
       Twall_P  = cyl(icyl)%exhaust_valves(ival)%twall_tube

       Dv       = cyl(icyl)%exhaust_valves(ival)%Dv
       Lvmax    = cyl(icyl)%exhaust_valves(ival)%Lvmax
       EVO      = cyl(icyl)%exhaust_valves(ival)%angle_VO
       EVC      = cyl(icyl)%exhaust_valves(ival)%angle_VC
       Nval     = cyl(icyl)%exhaust_valves(ival)%Nval
       type_dat = cyl(icyl)%exhaust_valves(ival)%type_dat
       call area_valve(theta, Dv, EVO, EVC, Lvmax, Nval, type_dat, &
            cyl(icyl)%exhaust_valves(ival)%Lv, cyl(icyl)%exhaust_valves(ival)%Cd, &
            nstroke, Area_T)
       Fev(ival) = Area_T

       if(cyl(icyl)%exhaust_valves(ival)%valve_model.eq.1 .or. myData%full_implicit) then
          call rhschar(Urev(:,ival), Area_P, dAreax_P, Twall_P, ga, R_gas, &
               globalData%dt, globalData%viscous_flow, globalData%heat_flow, RHS)
          RHSev(:,ival) = RHS
       end if

       if(myData%full_implicit) then
          ! call rhschar(Urev(:,ival), Area_P, dAreax_P, Twall_P, ga, R_gas, &
          !      globalData%dt, globalData%viscous_flow, globalData%heat_flow, RHS)
          ! RHSev(:,ival) = RHS
       else
          alpha = 1.          
          if(cyl(icyl)%exhaust_valves(ival)%valve_model.eq.0) then
             call solve_valve(globalData, Ucyl, Urev(:,ival), -1, Area_T, Area_P, &
                  Upipe, Uthroat)
             if(Area_T.gt.0.) alpha = 0.15
          else if(cyl(icyl)%exhaust_valves(ival)%valve_model.eq.1) then
             ! write(*,*) 'EXHAUST VALVE - CYL ', icyl
             Upipe = Uev(:,ival)
             call solve_valve_implicit(Ucyl, Urev(:,ival), -1, Area_T, Area_P, &
                  RHSev(:,ival), ga, Upipe, Uthroat)
          else
             stop ' Wrong valve model option '
          end if
          Uev(:,ival)  = alpha*Upipe + &
               (1.-alpha)*Uev(:,ival)
          Utev(:,ival) = Uthroat
       end if
       ! write(11,*) theta, Area_T, Uev, Utev
    end do

    if(myData%full_implicit) then
       ! Solves the cylinder and the valves in a monolythic system
       stop ' NOT IMPLEMENTED '
    else
       ! Solves the cylinder
       call cylinder_solver(icyl, myData, globalData, Ucyl, Utiv, Utev, &
            Uiv, Uev, Fiv, Fev, EGR, mass_cyl_old, mass_cyl)
    end if
    
    ! Actualizes the cylinder state(s)
    do i=0,nnod_cyl-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = Ucyl(j+1,i+1)
       enddo
    enddo

    ! Actualizes the intake valve state(s)
    do i=nnod_cyl,nnod_cyl+myData%nvi-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = Uiv(j+1,i-nnod_cyl+1)
       enddo
    enddo

    ! Actualizes the exhaust valve state(s)
    do i=nnod_cyl+myData%nvi,nnod_cyl+myData%nvi+myData%nve-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = Uev(j+1,i-nnod_cyl-myData%nvi+1)
       enddo
    enddo

    do i=1,nnod_cyl
       do j=1,3
          mass_C(6*(i-1)+j)   = mass_cyl(j,i)
          mass_C(6*(i-1)+j+3) = mass_cyl(j,i)
       end do
    end do
    
  end subroutine solve_cylinder

  subroutine geometry(myData, globalData, Vol, Area, Vdot)
    !
    !  Computes the surface chamber area, the volume and 
    !    the volume time derivative
    !
    !  geometry is called by: state_initial_cylinder, ignition_delay,
    !    cylinder_solver
    !  geometry calls the following subroutines and functions: none
    !
    implicit none

    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData
    real*8, intent(out) :: Vol,Area,Vdot

    integer :: nstroke
    real*8 :: beta,rpm,theta_g,theta
    real*8 :: Vc,Bo,l,a,Ach,Ap,s1,sdot1,s2,sdot2

    nstroke = globalData%nstroke
    rpm     = globalData%rpm
    theta_g = globalData%theta

    theta = modulo(theta_g+myData%theta_0, &
         nstroke*pi)

    ! Alternative engine geometry
    Vc  = myData%Vol_clearance     ! clearance volume
    Bo  = myData%Bore              ! Bore
    l   = myData%rod_length        ! connecting road length
    a   = myData%crank_radius      ! cranck radius
    Ach = myData%head_chamber_area ! cylinder head surface area
    Ap  = myData%piston_area       ! piston crown surface Area
    
    ! Computing the stroke (distance between the cranck axis and piston pin axis)
    s1 = a*dcos(theta)+dsqrt(l**2-a**2*dsin(theta)**2)
       
    sdot1 = -(a*dsin(theta)+a**2*dsin(2.0*theta)/2./ &
         dsqrt(l**2-a**2*dsin(theta)**2))
    if(globalData%engine_type.eq.1) then
       beta = myData%delta_ca
       s2 = a*dcos(theta-beta)+dsqrt(l**2-a**2*dsin(theta-beta)**2)
       sdot2 = -(a*dsin(theta-beta)+a**2*dsin(2.0*(theta-beta))/2./ &
            dsqrt(l**2-a**2*dsin(theta-beta)**2))
       Vol  = Vc + pi*Bo**2/4.0*(2.0*(l+a)-(s1+s2))
       Area = Ach + 2.0*Ap + pi*Bo*(2.0*(l+a)-(s1+s2))
       Vdot = -rpm*2.0*pi/60.0*pi*Bo**2/4.0*(sdot1+sdot2)
    else
       !  Computing the engine volume
       Vol  = Vc + pi*Bo**2/4*(l+a-s1)
       !  Computing the engine area
       Area = Ach + Ap + pi*Bo*(l+a-s1)
       !  Computing the volume rate (dV/dt)
       Vdot = -rpm*2*pi/60*pi*Bo**2/4*sdot1
    end if
    
  end subroutine geometry

  subroutine heat_transfer(myData, globalData, Ucyl, Area, dQ_ht)
    !
    !  Computes the heat losses in the cylinder
    !
    !  model = 1 -> Annand model
    !  model = 2 -> Woschni model 1
    !  model = 3 -> Woschni model 2
    !  model = 4 -> Taylor
    !
    !  heat_transfer is called by: cylinder_solver
    !  heat_transfer calls the following subroutines and functions: 
    !    compute_cp, compute_visco, average_piston_speed
    !
    implicit none

    real*8, intent(in) :: Area
    real*8, dimension(3), intent(in) :: Ucyl
    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData 
    real*8, intent(out) :: dQ_ht

    integer :: ispecie
    real*8 :: Bo,l,a,Twall,rho,p,T
    real*8 :: cp,mu,Pr,kappa,Sp,Re,Nu
    real*8 :: factor,C_h,C_r,dQ_hth,dQ_htr

    Bo    = myData%Bore         ! Bore
    l     = myData%rod_length   ! connecting road length
    a     = myData%crank_radius ! cranck radius
    Twall = myData%Twall        ! Wall temperature of the cylinder

    rho = Ucyl(1)
    p   = Ucyl(2)
    T   = Ucyl(3)

    ! computing Cp for a given mixture (isp) at a given temperature
    ispecie = 0
    cp = compute_cp(T,ispecie,globalData%R_gas)

    mu    = compute_visco(T) ! dynamic viscosity
    Pr    = 0.72             ! Prandtl number
    kappa = mu*Cp/Pr         ! thermal conductivity

    Sp = average_piston_speed(globalData%rpm, myData%crank_radius, &
         globalData%engine_type)
    Re = rho*dabs(Sp)*Bo/mu        ! Reynolds number

    ! Nusselt number
    ! According to Woschni but with a velocity different with respect
    ! to the piston speed in order to compute the Reynolds number,
    ! see Heywood pp. 678
    if(myData%model_ht.eq.1) then
       ! According to Annand (Heywood pp. 678)
       factor = myData%factor_ht * 0.49
       Nu     =  factor * Re**0.7
    elseif(myData%model_ht.eq.2) then
       ! factor is a constant times 0.035
       ! Some authors use a sensibility analysis in order to set this parameter.
       ! For instance, John Abraham uses in his code the value 4.5 times 0.035;
       ! Yacoub and Bata apply between 1 and 3 times this value in their
       ! '98 paper.
       factor = myData%factor_ht * 0.035
       Nu    = factor*Re**0.8*Pr**0.33
    elseif(myData%model_ht.eq.3) then
       factor = myData%factor_ht * 0.037
       Nu     = factor*Re**0.8*Pr**0.3
    else
       ! Taylor model
       factor = myData%factor_ht * 10.4
       Nu     = factor*Re**0.75
    endif

    ! Heat transfer coefficient
    C_h = Nu*kappa/Bo
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
   902 format (F12.4,F12.4,F12.4,F12.4)
    endif
    
  end subroutine heat_transfer

  subroutine SI_combustion(icyl, theta, mass_fuel, omega, nstroke, &
    x_burned, dQ_chem, dQ_ht_fuel, save_extras, nunit)
    !
    !  Computes the heat released by combustion for spark-ignition engines
    !    by a Wiebe function pp. 390 Heywood
    !
    !  SI_combustion is called by: heat_released
    !  SI_combustion calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: icyl,nstroke, nunit
    real*8, intent(in) :: theta,mass_fuel,omega
    real*8, intent(out) :: x_burned,dQ_chem,dQ_ht_fuel
    logical, intent(in) :: save_extras

    real*8 :: theta_ig,dtheta_comb,a_wiebe,m_wiebe,Q_fuel,hvap
    real*8 :: pot,xbdot
    a_wiebe     = cyl(icyl)%combustion_data%a_wiebe
    m_wiebe     = cyl(icyl)%combustion_data%m_wiebe
    theta_ig    = cyl(icyl)%combustion_data%theta_ig_0
    dtheta_comb = cyl(icyl)%combustion_data%dtheta_comb

    Q_fuel = cyl(icyl)%fuel_data%Q_fuel
    hvap   = cyl(icyl)%fuel_data%hvap_fuel

    pot      = -a_wiebe*(modulo(theta-theta_ig,nstroke*pi)/dtheta_comb)**(m_wiebe+1.)
    x_burned = 1.0d0 - dexp(pot)

    xbdot = dexp(pot)*a_wiebe*(m_wiebe+1.)* &
         (modulo(theta-theta_ig,nstroke*pi)/dtheta_comb)**m_wiebe/dtheta_comb
    xbdot = omega * xbdot

    dQ_chem    = mass_fuel * Q_fuel * xbdot

    dQ_ht_fuel = mass_fuel * xbdot * hvap

    if(save_extras) then
        write(nunit,904) x_burned, xbdot
    904 format (F12.4,F12.4)
    endif

  end subroutine SI_combustion

  subroutine CI_combustion(icyl, theta, mass_fuel, omega, nstroke, &
       x_burned, dQ_chem, dQ_ht_fuel, save_extras, nunit)
    !
    !  Computes the heat released by combustion for compression-ignition engines
    !
    !    Model 0 -> 'user-defined'
    !    Model 1 -> 'Wiebe-2'
    !    Model 2 -> 'Wiebe-3'
    !    Model 3 -> 'Watson'
    !
    !  CI_combustion is called by: heat_released
    !  CI_combustion calls the following subroutines and functions: interpolant
    !
    implicit none

    integer, intent(in) :: icyl,nstroke, nunit
    real*8, intent(in) :: theta,mass_fuel,omega
    real*8, intent(out) :: x_burned,dQ_chem,dQ_ht_fuel
    logical :: save_extras
    integer :: k,indx
    real*8 :: a,mp,mdi,pot,pote,theta_ig,dtheta_p,dtheta_di
    real*8 :: x_p,x_di,factor,dtheta_comb,xbdot,Q_fuel,hvap
    real*8 :: Fp,Fm,Ft,Ep,Em,Et,Dp,Dm,Dt,WCp,WCm,WCt,pot_p,pot_m,pot_t
    real*8 :: theta_id,phi_ig,tau,beta,theta_deg,rpm
    real*8, dimension(2) :: coef,coefp

    theta_ig    = cyl(icyl)%combustion_data%theta_ig_0
    dtheta_comb = cyl(icyl)%combustion_data%dtheta_comb

    Q_fuel = cyl(icyl)%fuel_data%Q_fuel
    hvap   = cyl(icyl)%fuel_data%hvap_fuel

    if(cyl(icyl)%combustion_data%combustion_model.eq.1) then
       !  Two Wiebe's functions to model the premixed and diffussive
       !    stages in the CI combustion
       a       = 5.
       mp      = 1.5
       mdi     = 1.05
       x_p     = 0.3
       x_di    = 1.-x_p
       factor  = 0.17857

       dtheta_p  = factor*dtheta_comb
       dtheta_di = dtheta_comb-dtheta_p
       pot       = -a*(modulo(theta-theta_ig,nstroke*pi)/dtheta_p)**(mp+1.)
       pote      = -a*(modulo(theta-theta_ig,nstroke*pi)/dtheta_di)**(mdi+1.)
       
       x_burned = 1.0d0 - x_p*dexp(pot) - x_di*dexp(pote)
       
       xbdot = x_p*dexp(pot)*a*(mp+1.)* &
            (modulo(theta-theta_ig,nstroke*pi)/dtheta_p)**mp/dtheta_p
       xbdot = xbdot + x_di*dexp(pote)*a*(mdi+1.)* &
            (modulo(theta-theta_ig,nstroke*pi)/dtheta_di)**mdi/dtheta_di
       xbdot = omega * xbdot
    elseif(cyl(icyl)%combustion_data%combustion_model.eq.2) then
       ! Three Wiebe's functions to model the premixed, main and tail
       ! stages in the CI combustion
       Fp = .02
       Ft = .05
       Fm = 1.-Fp-Ft

       Ep =  .7
       Et = 1.5
       Em =  .9

       Dp = 2.
       Dm = 35.
       Dt = 40.

       WCp = (Dp/2.302**(1./Ep) - .105**(1./Ep))**(-Ep)
       WCm = (Dm/2.302**(1./Em) - .105**(1./Em))**(-Em)
       WCt = (Dt/2.302**(1./Et) - .105**(1./Et))**(-Et)

       pot_p = -WCp*(modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**Ep
       pot_m = -WCm*(modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**Em
       pot_t = -WCt*(modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**Et
       
       x_burned = Fp*(1.0d0 - dexp(pot_p)) + Fm*(1.0d0 - dexp(pot_m)) + &
            Ft*(1.0d0 - dexp(pot_t))
       
       xbdot =         Fp*dexp(pot_p) * WCp * Ep * &
            (modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**(Ep-1.0)
       xbdot = xbdot + Fm*dexp(pot_m) * WCm * Em * &
            (modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**(Em-1.0)
       xbdot = xbdot + Ft*dexp(pot_t) * WCt * Et * &
            (modulo(theta-theta_ig,nstroke*pi)*180.0d0/pi)**(Et-1.0)
       xbdot = omega * xbdot * 180.0d0/pi
    elseif(cyl(icyl)%combustion_data%combustion_model.eq.3) then
       theta_id = cyl(icyl)%injection_data%theta_id
       phi_ig   = cyl(icyl)%combustion_data%phi_ig
       rpm      = 30.*omega/pi

       coef     = (/ 2.0d0 + 1.25d-8*(500.0/3.0*theta_id)**2.4, 5000d0/)
       coefp    = (/ 14.2d0/phi_ig**0.644, 0.79d0*(14.2d0/phi_ig**0.644)**0.25 /)
       beta     = 1d0 - 0.926d0*phi_ig**0.37/(500.0/3.0*theta_id/rpm)**0.26

       tau = modulo(theta-theta_ig*pi/180.,nstroke*pi)/dtheta_comb

       x_burned = beta*(1.-(1.-tau**coef(1))**coef(2)) + &
            (1.-beta)*(1.-dexp(-coefp(1)*tau**coefp(2)))

       xbdot = omega/dtheta_comb*(beta*coef(1)*coef(2)*tau**(coef(1)-1.)* &
            (1.-tau**coef(1))**(coef(2)-1) + &
            (1.-beta)*coefp(1)*coefp(2)*tau**(coefp(2)-1.)*dexp(-coefp(1)*tau**coefp(2)))
    elseif(cyl(icyl)%combustion_data%combustion_model.eq.0) then
       theta_deg = theta*180./pi

       call interpolant(cyl(icyl)%combustion_data%xbdot_array(:,1), &
            cyl(icyl)%combustion_data%xbdot_array(:,2), theta_deg, xbdot)
       k = iminloc(dabs(cyl(icyl)%combustion_data%xbdot_array(:,2)-xbdot))
       if(cyl(icyl)%combustion_data%xbdot_array(k,2).gt.xbdot) &
            k = k-1
       x_burned = 0.
       indx     = k-1
       do k=1,indx
          x_burned = x_burned + .5* &
               (cyl(icyl)%combustion_data%xbdot_array(k+1,1) - &
               cyl(icyl)%combustion_data%xbdot_array(k,1))* &
               (cyl(icyl)%combustion_data%xbdot_array(k+1,2) + &
               cyl(icyl)%combustion_data%xbdot_array(k,2))
       end do
       x_burned = x_burned + .5* &
            (theta_deg-cyl(icyl)%combustion_data%xbdot_array(k,1))* &
            (xbdot+cyl(icyl)%combustion_data%xbdot_array(k,2))

       xbdot = omega * xbdot * 180./pi
    end if

    dQ_chem = mass_fuel * Q_fuel * xbdot

    dQ_ht_fuel = mass_fuel * xbdot * hvap

    if(save_extras) then
        write(nunit,903) x_burned, xbdot
    903 format (F12.4,F12.4)
    endif
    
  end subroutine CI_combustion

  subroutine comp_scavenge(icyl, mass_cyl, SE)
    !
    !
    !  comp_scavenge is called by: cylinder_solver
    !  comp_scavenge calls the following subroutines and functions: fa_ratio
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, dimension(3), intent(inout) :: mass_cyl
    real*8, intent(out) :: SE

    real*8 :: M,C,mass_tr,SRv,phi,y,factor

    M   = cyl(icyl)%scavenge_data%val_1
    C   = cyl(icyl)%scavenge_data%val_2
    SRv = cyl(icyl)%scavenge_data%SRv

    phi = cyl(icyl)%combustion_data%phi
    y   = cyl(icyl)%fuel_data%y

    factor = fa_ratio(phi,y)
    
    SE = dmax1(1.0d0 - dexp(M*SRv+C),0.0d0)
    if(cyl(icyl)%scavenge_data%close_cyl.and.(SRv.gt.0.)) then
       mass_tr     = sum(mass_cyl(1:3))
       mass_cyl(2) = mass_tr*SE
       mass_cyl(1) = mass_cyl(2)*factor
       mass_cyl(3) = mass_tr-(mass_cyl(1)+mass_cyl(2))
       SRv = 0.0d0
       mass_cyl(3) = dmax1(mass_cyl(3),0.0d0)
    end if
    cyl(icyl)%scavenge_data%SRv = SRv

  end subroutine comp_scavenge

  subroutine scavenge_ratio(UC, UP, F_P, V_C, dt, SRv)
    !
    !  Computes the scavenge ratio by volume
    !
    !  scavenge_ratio is called by: 
    !  scavenge_ratio calls the following subroutines and functions: none
    !
    implicit none

    real*8, intent(in) :: F_P,V_C,dt
    real*8, dimension(3), intent(in) :: UC,UP
    real*8, intent(inout) :: SRv

    real*8 :: u_P,p_P,p_C,expand

    u_P = UP(2)
    p_P = UP(3)
    p_C = UC(2)

    expand = p_P/p_C

    SRv = SRv + expand*dt*u_P*F_P/V_C

  end subroutine scavenge_ratio

  subroutine ignition_delay(icyl, myData, globalData, Ucyl)
    !
    !
    !
    !  ignition_delay is called by: cylinder_solver
    !  ignition_delay calls the following subroutines and functions:
    !    geometry, average_piston_speed
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, dimension(3), intent(in) :: Ucyl
    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData

    integer :: nstroke,ignition_delay_model
    real*8 :: theta_g,theta,rpm,dt
    real*8 :: theta_id,angle_aux,theta_inj_ini,dtheta_comb
    real*8 :: A_d,n_d,Ta,t_d,Vch,Vol,Area,Vdot
    real*8 :: nc,T_TC,p_TC,Sp,CN,Ea

    dt      = globalData%dt
    rpm     = globalData%rpm
    nstroke = globalData%nstroke
    theta_g = globalData%theta

    theta_id             = cyl(icyl)%injection_data%theta_id
    theta_inj_ini        = cyl(icyl)%injection_data%theta_inj_ini
    dtheta_comb          = cyl(icyl)%combustion_data%dtheta_comb
    ignition_delay_model = cyl(icyl)%injection_data%ignition_delay_model

    theta = modulo(theta_g+myData%theta_0, &
         nstroke*pi)

    angle_aux = 0d0
    if(ignition_delay_model.eq.0) then
       ! 'Integral' ignition delay model by Livengood and Wu
       if(theta_id.eq.0. .and. modulo(theta-theta_inj_ini,nstroke*pi) .ge. 0. .and. &
            modulo(theta-theta_inj_ini,nstroke*pi) .le. dtheta_comb) then
          A_d = 3.45d0
          n_d = 1.02d0
          Ta  = 2100.0d0
          t_d = A_d/((Ucyl(2)*1.0197162d-5)**n_d)*dexp(Ta/Ucyl(3))*1d-3
          if(cyl(icyl)%injection_data%integral + dt/t_d .lt. 1.) then
             cyl(icyl)%injection_data%integral = cyl(icyl)%injection_data%integral + dt/t_d
             angle_aux = dtheta_comb
          else
             angle_aux = 0d0
             theta_id = theta-theta_inj_ini
          endif
       endif
    elseif(ignition_delay_model.eq.1) then
       ! Correlation for the ignition delay proposed by Hardenber and Hase
       Vch = myData%Vol_clearance
       if(theta_id.eq.0. .and. modulo(theta,nstroke*pi)-theta_inj_ini .ge. 0. .and. &
            modulo(theta,nstroke*pi)-theta_inj_ini .le. dtheta_comb) then
          ! Average piston speed
          Sp = average_piston_speed(rpm, myData%crank_radius, globalData%engine_type)
          CN = 40.
          Ea = 618840./(CN+25.)
          ! We estimates the temperature and pressure at TDC from the start of injection
          call geometry(myData, globalData, Vol, Area, Vdot)
          nc = 1.4
          T_TC = Ucyl(3)*(Vol/Vch)**(nc-1.)
          p_TC = Ucyl(2)*(Vol/Vch)**nc
          theta_id = (0.36d0+0.22d0*Sp)*dexp(Ea*(1d0/(8.3143d0*T_TC)-1d0/17190d0)* &
               (21.2d0/(1d-5*p_TC-12.4d0))**0.63d0)
          theta_id = theta_id*pi/180.
          write(*,*) theta_id, T_TC, p_TC
       endif
    elseif(ignition_delay_model.eq.2) then
       ! A constant user-defined ignition delay angle
    else
       stop ' Unknown ignition delay model '
    endif

    cyl(icyl)%combustion_data%theta_ig_0 = &
         modulo(theta_inj_ini+theta_id+angle_aux,nstroke*pi)
    cyl(icyl)%injection_data%theta_id    = theta_id

  end subroutine ignition_delay

  subroutine restart_injection_var(icyl, nstroke, theta, theta_id, &
       theta_ig_0, dtheta_comb, mass_cyl)
    !
    !
    !
    !  restart_injection_var is called by: cylinder_solver
    !  restart_injection_var calls the following subroutines and functions: none
    !
    implicit none

    integer, intent(in) :: icyl,nstroke
    real*8, intent(in) :: theta,theta_id,theta_ig_0,dtheta_comb
    real*8, dimension(3), intent(inout) :: mass_cyl

    if(theta_id.ne.0.0 .and. & 
         modulo(theta-theta_ig_0,nstroke*pi) .gt. dtheta_comb .and. &
         .not.(cyl(icyl)%combustion_data%start_comb)) then
       mass_cyl(1)     = 0.0d0

       cyl(icyl)%injection_data%theta_id    = 0.0d0
       cyl(icyl)%injection_data%integral    = 0d0
       cyl(icyl)%combustion_data%phi        = 0d0
       cyl(icyl)%combustion_data%start_comb = .true.
    end if

  end subroutine restart_injection_var

  subroutine heat_released(icyl, type_ig, theta, omega, nstroke, &
       dQ_chem, dQ_ht_fuel, mass_cyl, save_extras, nunit)
    !
    !
    !
    !  heat_released is called by: cylinder_solver
    !  heat_released calls the following subroutines and functions:
    !    SI_combustion, CI_combustion
    !
    implicit none

    integer, intent(in) :: icyl,nstroke,type_ig,nunit
    real*8, intent(in) :: theta,omega
    real*8, dimension(3), intent(inout) :: mass_cyl
    real*8, intent(out) :: dQ_chem,dQ_ht_fuel
    logical, intent(in) :: save_extras
    real*8 :: dtheta_comb,theta_ig,xb
    real*8 :: mass_fuel,mass_air,mass_res

    dQ_chem    = 0.
    dQ_ht_fuel = 0.

    dtheta_comb = cyl(icyl)%combustion_data%dtheta_comb
    theta_ig    = cyl(icyl)%combustion_data%theta_ig_0
    !write(*,*) "cyl: ",icyl
    !write(*,*) "modulo(theta-theta_ig,nstroke*pi) (ge 0): ",modulo(theta-theta_ig,nstroke*pi)
    !write(*,*) "menor igual: ",modulo(theta-theta_ig,nstroke*pi) ,modulo(theta-(theta_ig+dtheta_comb),nstroke*pi)
    if(modulo(theta-theta_ig,nstroke*pi).ge.0. .and. &
         modulo(theta-theta_ig,nstroke*pi).le. &
         modulo(theta-(theta_ig+dtheta_comb),nstroke*pi)) then
       mass_fuel = cyl(icyl)%combustion_data%mass_fuel_ini
       mass_air  = cyl(icyl)%combustion_data%mass_air_ini
       mass_res  = cyl(icyl)%combustion_data%mass_res_ini

       if(type_ig.eq.0) then
          call SI_combustion(icyl, theta, mass_fuel, omega, nstroke, &
               xb, dQ_chem, dQ_ht_fuel,save_extras, nunit)
          cyl(icyl)%combustion_data%start_comb = .false.
       else
          if(cyl(icyl)%combustion_data%start_comb) &
               cyl(icyl)%combustion_data%phi_ig = cyl(icyl)%combustion_data%phi
          call CI_combustion(icyl, theta, mass_fuel, omega, nstroke, &
               xb, dQ_chem, dQ_ht_fuel, save_extras, nunit)
          cyl(icyl)%combustion_data%start_comb = .false.
       endif

       mass_cyl(1) = (1.-xb)*mass_fuel
       mass_cyl(2) = (1.-xb)*mass_air
       mass_cyl(3) = mass_res+xb*(mass_fuel+mass_air)
       
       mass_cyl(1) = dmax1(mass_cyl(1),0.0d0)
       mass_cyl(2) = dmax1(mass_cyl(2),0.0d0)
       mass_cyl(3) = dmax1(mass_cyl(3),0.0d0)
       
    else
       if(.not.cyl(icyl)%combustion_data%start_comb) &
            cyl(icyl)%combustion_data%start_comb = .true.
       if(save_extras) then
          write(nunit,905) 0,0
      905 format(I12,I12)
       endif
    end if

  end subroutine heat_released

  subroutine fuel_injection(icyl, nstroke, theta, rpm, dt, mass_cyl_old, mass_cyl)
    !
    !
    !
    !  fuel_injection is called by: cylinder_solver
    !  fuel_injection calls the following subroutines and functions:
    !    mass_injection, fa_ratio
    !
    implicit none

    integer, intent(in) :: icyl,nstroke
    real*8, intent(in) :: theta,rpm,dt
    real*8, dimension(3), intent(in) :: mass_cyl_old
    real*8, dimension(3), intent(inout) :: mass_cyl

    real*8 :: theta_inj_ini,dtheta_inj,delta_mass_fuel
    real*8 :: f_st,phi
    
    theta_inj_ini = cyl(icyl)%injection_data%theta_inj_ini
    dtheta_inj    = cyl(icyl)%injection_data%dtheta_inj

    delta_mass_fuel = 0.0

    if(modulo(theta-theta_inj_ini,nstroke*pi) .ge. 0.0 .and. &
         modulo(theta-theta_inj_ini,nstroke*pi) .le. dtheta_inj) then
       call mass_injection(icyl, theta, rpm, dt, delta_mass_fuel)
       
       f_st = fa_ratio(1.0d0, cyl(icyl)%fuel_data%y)
       phi  = cyl(icyl)%combustion_data%phi
       phi  = phi + (1.+phi*f_st)/((mass_cyl(1)+mass_cyl(2))*f_st)* &
            (delta_mass_fuel - 0.*phi*f_st*(mass_cyl(2)-mass_cyl_old(2)))
       phi  = dmax1(phi,0.0d0)
       cyl(icyl)%combustion_data%phi = phi
    end if

    mass_cyl(1) = mass_cyl(1)+delta_mass_fuel

  end subroutine fuel_injection

  subroutine mass_injection(icyl, theta, rpm, dt, dmass)
    !
    !
    !
    !  mass_injection is called by: fuel_injection
    !  mass_injection calls the following subroutines and functions:
    !    interpolant
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, intent(in) :: theta,rpm,dt
    real*8, intent(out) :: dmass

    real*8 :: m_inj,theta_inj_ini,dtheta_inj
    real*8 :: mass_dot,omega
 
    omega = rpm*pi/ 30.

    theta_inj_ini = cyl(icyl)%injection_data%theta_inj_ini
    dtheta_inj    = cyl(icyl)%injection_data%dtheta_inj
    m_inj         = cyl(icyl)%injection_data%m_inj

    if(cyl(icyl)%injection_data%pulse.eq.1) then
       mass_dot = 2.0d0*m_inj*omega/dtheta_inj* &
            (dsin(pi*theta/dtheta_inj))**2
    elseif(cyl(icyl)%injection_data%pulse.eq.2) then
       mass_dot = m_inj*omega/dtheta_inj
    elseif(cyl(icyl)%injection_data%pulse.eq.3) then
       call interpolant(cyl(icyl)%injection_data%mfdot_array(:,1), &
            cyl(icyl)%injection_data%mfdot_array(:,2), theta*180./pi, mass_dot)
    else
       stop 'Bad pulse type for fuel injection'
    end if

    dmass = dt * mass_dot

  end subroutine mass_injection

  subroutine comp_cyl_masses(icyl, dt, mass_in, mass_out, EGR, mass_cyl_old, &
       mass_cyl)
    !
    !
    !
    !  comp_cyl_masses is called by: cylinder_solver
    !  comp_cyl_masses calls the following subroutines and functions:
    !    fa_ratio
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, intent(in) :: dt,mass_in,mass_out,EGR
    real*8, dimension(3), intent(in) :: mass_cyl_old
    real*8, dimension(3), intent(out) :: mass_cyl
    
    real*8 :: mass_old,xr,phi,y,factor
    real*8 :: delta_mass_air,delta_mass_fuel,delta_mass_res

    mass_old = sum(mass_cyl_old)
    xr       = mass_cyl_old(3)/mass_old

    phi = cyl(icyl)%combustion_data%phi
    y   = cyl(icyl)%fuel_data%y

    factor = fa_ratio(phi,y)

    if(mass_in.gt.0.) then
       delta_mass_air = dt*(mass_in*(1.-EGR) + &
            (1.-xr)*mass_out)/(1.+factor)
    else
       delta_mass_air = dt*(mass_in*(1.-xr) + &
            (1.-xr)*mass_out)/(1.+factor)
    endif
    
    delta_mass_fuel = factor*delta_mass_air
    if(mass_in.lt.0.) then
       delta_mass_res = xr*dt*(mass_out+mass_in)
    else
       delta_mass_res = dt*(xr*mass_out+EGR*mass_in)
    endif
    
    mass_cyl(1) = dmax1(mass_cyl_old(1)+delta_mass_fuel,0.0d0)
    mass_cyl(2) = dmax1(mass_cyl_old(2)+delta_mass_air,0.0d0)
    mass_cyl(3) = dmax1(mass_cyl_old(3)+delta_mass_res,0.0d0)

  end subroutine comp_cyl_masses

  subroutine cylinder_solver(icyl, myData, globalData, Ucyl, Uiv, Uev, &
       Upiv, Upev, Fiv, Fev, EGR, mass_cyl_old, mass_cyl)
    !
    !
    !
    !  cylinder_solver is called by: solve_cylinder
    !  cylinder_solver calls the following subroutines and functions:
    !    valve_flow, geometry, heat_transfer, comp_cyl_masses,
    !    ignition_delay, restart_injection_var, fuel_injection,
    !    heat_released, compute_cp
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, intent(in) :: EGR
    real*8, dimension(3), intent(in) :: mass_cyl_old
    real*8, dimension(cyl(icyl)%nvi), intent(in) :: Fiv
    real*8, dimension(cyl(icyl)%nve), intent(in) :: Fev
    real*8, dimension(3,cyl(icyl)%nvi), intent(in) :: Uiv,Upiv
    real*8, dimension(3,cyl(icyl)%nve), intent(in) :: Uev,Upev
    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData 
    real*8, dimension(3), intent(inout) :: Ucyl
    real*8, dimension(3), intent(out) :: mass_cyl

    integer :: ispecie,nstroke,type_ig, nunit, engine_type
    real*8 :: dt,rpm,theta_g,theta,crank_radius
    real*8 :: cp,cv,Vol,Area,Vdot,mass_old,mass_new
    real*8 :: omega,rho_cyl,p_cyl,T_cyl
    real*8 :: dQ_ht,dQ_chem,dQ_ht_fuel
    real*8 :: edot,ene_old,ene_new
    real*8 :: mass_in,mass_out,h_in,h_out
    real*8 :: SE, Torque
    logical :: save_extras

	nstroke = globalData%nstroke
    dt      = globalData%dt
    rpm     = globalData%rpm
    theta_g = globalData%theta
    save_extras = globalData%save_extras
    engine_type = globalData%engine_type
    nunit =  myData%nunit
    type_ig = myData%type_ig


    theta = modulo(theta_g+myData%theta_0, &
         nstroke*pi)

    omega = rpm*pi/30.

    rho_cyl = Ucyl(1)
    p_cyl   = Ucyl(2)
    T_cyl   = Ucyl(3)

    ispecie = 0
    cp = compute_cp(T_cyl, ispecie, globalData%R_gas)
    cv = cp - globalData%R_gas

    call valve_flow(cyl(icyl)%nvi, 1, globalData%ga, Fiv, Uiv, Upiv, Ucyl, &
         mass_in, h_in)
    call valve_flow(cyl(icyl)%nve, -1, globalData%ga, Fev, Uev, Upev, Ucyl, &
         mass_out, h_out)

    ! Engine geometrical data
    call geometry(myData, globalData, Vol, Area, Vdot)
    ! Computing heat losses
    call heat_transfer(myData, globalData, Ucyl, Area, dQ_ht)
    ! Computing the fuel mass, air mass and gas residual mass into the cylinder
    call comp_cyl_masses(icyl, dt, mass_in, mass_out, EGR, mass_cyl_old, mass_cyl)
    ! Computing scavenge
    if(myData%scavenge) call comp_scavenge(icyl, mass_cyl, SE)
    if(myData%type_ig.eq.1) then
       ! Computing the ignition delay (for CI engines only)
       call ignition_delay(icyl, myData, globalData, Ucyl)
       call restart_injection_var(icyl, nstroke, theta, cyl(icyl)%injection_data%theta_id, &
            cyl(icyl)%combustion_data%theta_ig_0, cyl(icyl)%combustion_data%dtheta_comb, &
            mass_cyl)
       ! Computing the fuel injection (for CI engines only)
       call fuel_injection(icyl, nstroke, theta, rpm, dt, mass_cyl_old, mass_cyl)
    end if
    ! Computing the heat released by combustion
    if(cyl(icyl)%combustion_data%start_comb) then
       cyl(icyl)%combustion_data%mass_fuel_ini = mass_cyl_old(1)
       cyl(icyl)%combustion_data%mass_air_ini  = mass_cyl_old(2)
       cyl(icyl)%combustion_data%mass_res_ini  = mass_cyl_old(3)
    end if
    call heat_released(icyl, type_ig, theta, omega, nstroke, &
         dQ_chem, dQ_ht_fuel, mass_cyl, save_extras, nunit)

    mass_old = sum(mass_cyl_old)
    mass_new = sum(mass_cyl)

    edot    = (h_in+h_out)-p_cyl*Vdot+dQ_chem-(dQ_ht+dQ_ht_fuel)
    ene_old = mass_old*cv*T_cyl
    ene_new = ene_old + dt*edot
    
    ! Continuity equation
    rho_cyl = mass_new/Vol
    ! Energy equation
    T_cyl  = ene_new/(mass_new*cv)
    ! State equation
    p_cyl  = rho_cyl*globalData%R_gas*T_cyl
    
    Ucyl(1) = rho_cyl
    Ucyl(2) = p_cyl
    Ucyl(3) = T_cyl

    Torque = 0.0
    call fun_cyl_Torque(myData, globalData, p_cyl, theta, Torque)
    
    if(globalData%save_extras) then
       write(myData%nunit,900) mass_in, mass_out, Vol, mass_cyl, dQ_ht, dQ_chem, Torque
    900 format (E12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4)
    endif

  end subroutine cylinder_solver

  subroutine run_ideal_cycle(icyl, myData, globalData, atm_state, AVPT)
    !
    !
    !
    !  run_ideal_cycle is called by: state_initial_cylinder
    !  run_ideal_cycle calls the following subroutines and functions:
    !    fa_ratio, ideal_cycle_4S, ideal_cycle_2S
    !
    implicit none

    integer, intent(in) :: icyl
    real*8, dimension(3), intent(in) :: atm_state
    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData
    real*8, dimension(:,:), intent(out) :: AVPT
    
    real*8 :: Bo,a,Vd,Vc,Vmin,Vmax
    real*8 :: rho_i,T_i,p_i,phi,y,ga,R_gas,Q_fuel,factor,AF

    ga    = globalData%ga
    R_gas = globalData%R_gas

    Vc  = myData%Vol_clearance  ! Clearance volume
    Bo  = myData%Bore           ! Bore
    a   = myData%crank_radius   ! Cranck radius

    Q_fuel = cyl(icyl)%fuel_data%Q_fuel
    phi    = cyl(icyl)%combustion_data%phi
    y      = cyl(icyl)%fuel_data%y

    Vd   = a*pi*Bo**2/2.
    if(globalData%engine_type.eq.1) &
         Vd = 2.*Vd
    Vmin = Vc
    Vmax = Vc+Vd

    factor = fa_ratio(phi,y)
    AF     = 1./factor

    rho_i = atm_state(1)
    p_i   = atm_state(3)
    T_i   = p_i/(R_gas*rho_i)

    if(globalData%nstroke.eq.4) then
       call ideal_cycle_4S(Vmin, Vmax, R_gas, ga, 0.3*Q_fuel, T_i, p_i, &
            AF, AVPT)
    elseif(globalData%nstroke.eq.2) then
       call ideal_cycle_2S(Vmin, Vmax, R_gas, ga, 0.*Q_fuel, T_i, p_i, &
            AF, AVPT)
    endif

  end subroutine run_ideal_cycle

  subroutine ideal_cycle_4S(Vmin, Vmax, R_gas, ga, Q_LHV, T_i, p_i, &
       AF, AVPT)
    !
    !  AVPT = [angle , density , pressure , temperature ] 
    !
    !  ideal_cycle_4S is called by: run_ideal_cycle
    !  ideal_cycle_4S calls the following subroutines and functions: none
    !
    implicit none

    real*8, intent(in) :: R_gas,ga,Vmin,Vmax,Q_LHV,T_i,p_i,AF
    real*8, dimension(7,4), intent(out) :: AVPT

    real*8 :: cv
    real*8 :: theta1,theta2,theta3,theta4,theta5,theta6,theta7
    real*8 :: rho1,rho2,rho3,rho4,rho5,rho6,rho7
    real*8 :: V1,V2,V3,V4,V5,V6,V7
    real*8 :: p1,p2,p3,p4,p5,p6,p7
    real*8 :: T1,T2,T3,T4,T5,T6,T7
    real*8 :: m1,m2,m3,m4,m5,m6,m7,mf

    cv = R_gas/(ga-1.)

    ! point 1
    V1     = Vmin
    p1     = p_i
    T1     = T_i
    rho1   = p1/R_gas/T1
    m1     = rho1*V1
    theta1 = 0.

    ! point 2
    V2     = Vmax
    p2     = p1
    T2     = T1
    rho2   = rho1
    m2     = rho2*V2
    theta2 = 180.

    mf = m2/AF
    ! point 3
    V3     = Vmin
    m3     = m2
    rho3   = m3/V3
    p3     = p2*(rho3/rho2)**ga
    T3     = p3/R_gas/rho3
    theta3 = 360.-1.

    ! point 4
    V4     = Vmin
    m4     = m3
    rho4   = rho3
    T4     = T3+mf*Q_LHV/m4/cv
    p4     = R_gas*rho4*T4
    theta4 = 360.+1.

    ! point 5
    V5     = Vmax
    m5     = m4
    rho5   = m5/V5
    p5     = p4*(rho5/rho4)**ga
    T5     = p5/R_gas/rho5
    theta5 = 540.-1.

    ! point 6
    V6     = Vmax
    p6     = p_i
    T6     = T_i
    rho6   = p6/R_gas/T6
    m6     = rho6*V6
    theta6 = 540.+1.

    ! point 7
    V7     = V1
    p7     = p1
    T7     = T1
    rho7   = rho1
    m7     = m1
    theta7 = 720.

    AVPT(:,1) = (/theta1,theta2,theta3,theta4,theta5,theta6,theta7/)
    AVPT(:,2) = (/rho1,rho2,rho3,rho4,rho5,rho6,rho7/)
    AVPT(:,3) = (/p1,p2,p3,p4,p5,p6,p7/)
    AVPT(:,4) = (/T1,T2,T3,T4,T5,T6,T7/)
    
  end subroutine ideal_cycle_4S

  subroutine ideal_cycle_2S(Vmin, Vmax, R_gas, ga, Q_LHV, T_i, p_i, &
       AF, AVPT)
    !
    !      AVPT = [angle , density , pressure , temperature ] 
    !
    !  ideal_cycle_2S is called by: run_ideal_cycle
    !  ideal_cycle_2S calls the following subroutines and functions: none
    !
    implicit none

    real*8, intent(in) :: R_gas,ga,Vmin,Vmax,Q_LHV,T_i,p_i,AF
    real*8, dimension(7,4), intent(out) :: AVPT

    real*8 :: cv,mass,mf,rho_i,Vivc,m1
    real*8 :: theta1,theta2,theta3,theta4,theta5,theta6,theta7
    real*8 :: rho1,rho2,rho3,rho4,rho5,rho6,rho7
    real*8 :: V1,V2,V3,V4,V5,V6,V7
    real*8 :: p1,p2,p3,p4,p5,p6,p7
    real*8 :: T1,T2,T3,T4,T5,T6,T7

    cv = R_gas/(ga-1)

    Vivc  = 2./3.*Vmax ! Hay que calcularlo bien!!!!
    rho_i = p_i/(R_gas*T_i)
    mass  = rho_i*Vivc
    mf    = mass/AF

    ! point 1
    V1     = Vmin
    m1     = mass
    rho1   = m1/V1
    T1     = p_i*(rho1/rho_i)**ga/(R_gas*rho1)+mf*Q_LHV/m1/cv
    p1     = R_gas*rho1*T1
    theta1 = 0.

    ! point 2
    V2     = Vivc
    rho2   = m1/V2
    p2     = p1*(rho2/rho1)**ga
    T2     = p2/(R_gas*rho2)
    theta2 = 135.-1.

    ! point 3
    V3     = Vivc
    p3     = p_i
    T3     = T_i
    rho3   = p3/R_gas/T3
    theta3 = 135.+1.

    ! point 4
    V4     = Vmax
    p4     = p_i
    T4     = T_i
    rho4   = rho_i
    theta4 = 180.-1.

    ! point 5
    V5     = Vmax
    p5     = p_i
    T5     = T_i
    rho5   = rho_i
    theta5 = 180.+1.

    ! point 6
    V6     = Vivc
    rho6   = rho_i
    p6     = p_i
    T6     = T_i
    theta6 = 225.

    ! point 7
    V7     = Vmin
    rho7   = mass/Vmin
    p7     = p_i*(rho7/rho_i)**ga
    T7     = p7/(R_gas*rho7)
    theta7 = 360.

    AVPT(:,1) = (/theta1,theta2,theta3,theta4,theta5,theta6,theta7/)
    AVPT(:,2) = (/rho1,rho2,rho3,rho4,rho5,rho6,rho7/)
    AVPT(:,3) = (/p1,p2,p3,p4,p5,p6,p7/)
    AVPT(:,4) = (/T1,T2,T3,T4,T5,T6,T7/)

  end subroutine ideal_cycle_2S

  function average_piston_speed(rpm, a, engine_type)
    !
    !  Computes the average piston speed
    !
    implicit none

    integer, intent(in) :: engine_type
    real*8, intent(in) :: rpm,a
    real*8 :: average_piston_speed

    average_piston_speed = rpm/30.*(2.*a)
    if(engine_type.eq.1) &
         average_piston_speed = 2.*average_piston_speed

  end function average_piston_speed

  function fa_ratio(phi, y)
    !
    !
    !
    implicit none

    real*8, intent(in) :: phi,y
    real*8 :: fa_ratio

    fa_ratio = (12.011+1.008*y)/(34.56*(4+y))*phi

  end function fa_ratio


  subroutine fun_cyl_Torque(myData, globalData, pcyl, theta, Torque)

    implicit none

    real*8, intent(in) :: pcyl, theta
    type(this), intent(in) :: myData
    type(dataSim), intent(in) :: globalData 
    real*8, intent(inout) :: Torque
    
    real*8 :: lambda,Area_piston
    real*8 :: crank_radius,rod_length
    real*8 :: exce,piston_mass,rod_mass,rod_a,rod_b,mass_A,mass_B,mass_t
    real*8 :: sin_phi,cos_phi,tan_phi
    real*8 :: kphi,ks,kphiprime,ksprime
    real*8 :: inercia_J_AB,crankcase_pressure, inercia_crank	

    Area_piston  = pi/4.0*myData%Bore**2
    crank_radius = myData%crank_radius
    rod_length   = myData%rod_length
    lambda       = crank_radius/rod_length
    exce         = 0.0
    piston_mass = 0.0
    rod_mass    = 0.0
    rod_a       = 0.0*rod_length;  
    rod_b       = rod_length-rod_a

    mass_B      = rod_a/rod_length*rod_mass 
    mass_A      = rod_b/rod_length*rod_mass
    mass_t      = piston_mass+mass_B  

    inercia_J_AB  = 0.0
    inercia_crank = 0.0

    !theta     = globalData%theta_sim1d
    !theta_cyl = globalData%theta+myData%theta_0;
    ! Modificado Black (04/12/07)
    ! theta_cyl = theta_cyl-floor(theta_cyl/(4*pi))*(4*pi);
    !theta_cyl = modulo(theta_cyl,myData%nstroke*pi)
    ! Fin modificado Black (04/12/07)
    !theta     = theta_cyl

    crankcase_pressure = 0.0

    sin_phi   = (lambda*dsin(theta)-exce/rod_length)
    cos_phi   = dsqrt(1.0d0-sin_phi**2)
    tan_phi   = sin_phi/cos_phi
    kphi      = lambda*dcos(theta)/dsqrt(1.0d0-(lambda*dsin(theta)-exce/rod_length)**2);	                
    ks        = -crank_radius*(dsin(theta)+kphi/lambda*sin_phi);
    kphiprime = -lambda*dsin(theta)/cos_phi+kphi**2*tan_phi;
    ksprime   = -crank_radius*(dcos(theta) + & 
         1.0d0/lambda*kphiprime*sin_phi+1.0d0/lambda*kphi**2*cos_phi);
    ! Agregado Black (04/12/07)
    ! Ojo, solo para calcular el torque del gas
    ! No chequeado para calcular la dinamica del mecanismo
    if(globalData%engine_type .eq. 1) then
       sin_phi   = (lambda*dsin(theta-myData%delta_ca)-exce/rod_length)
       cos_phi   = dsqrt(1.0d0-sin_phi**2)
       kphi      = lambda*dcos(theta-myData%delta_ca)/ &
            dsqrt(1.0d0-(lambda*dsin(theta-myData%delta_ca)-exce/rod_length)**2)
       ks        = ks - crank_radius*(dsin(theta-myData%delta_ca)+kphi/lambda*sin_phi)
    end if
    ! Fin agregado Black (04/12/07)

    Torque = (pcyl-crankcase_pressure)*Area_piston*ks;

  end subroutine fun_cyl_Torque

end module def_cylinder
