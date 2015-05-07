module def_simulator
  use, intrinsic :: ISO_C_BINDING
  
  type, BIND(C) :: dataSim
     real(C_DOUBLE) :: time, dt, dtheta_rpm, Courant, ga, R_gas, &
          theta, omega, rpm, rpm_ini, theta_cycle
     real(C_DOUBLE) :: ga_intake, ga_exhaust
     integer(C_INT) :: iter_sim1d, viscous_flow, heat_flow, icycle, &
          engine_type, ncyl
     logical(C_BOOL) :: save_extras, use_global_gas_prop
  end type dataSim
  
end module def_simulator
