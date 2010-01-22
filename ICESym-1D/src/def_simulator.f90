module def_simulator
  use, intrinsic :: ISO_C_BINDING
     
    type, BIND(C) :: dataSim
    real(C_DOUBLE) :: time, dt, dtheta_rpm, Courant, ga, R_gas, &
         theta, omega, rpm, rpm_ini
    integer(C_INT) :: iter_sim1d, viscous_flow, heat_flow, icycle, nstroke, engine_type
	logical(C_BOOL) :: save_extras
    
 end type dataSim

end module def_simulator
