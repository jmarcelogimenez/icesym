typedef struct
	{
	  double time, dt, dtheta_rpm, Courant, ga, R_gas, theta, omega, rpm, rpm_ini;
	int iter_sim1d, viscous_flow, heat_flow, icycle, nstroke, engine_type;
	bool save_extras;
	}dataSim;

typedef struct
	{   
	int nnod_input, nvi,nve,nnod,ndof, model_ht, type_ig, nunit,species_model,ntemp;
	double Bore,crank_radius,Vol_clearance,rod_length,head_chamber_area,piston_area,theta_0,delta_ca,Twall,factor_ht;
	bool scavenge, full_implicit,nh_temp;
	}dataCylinder;


typedef struct
	{   
	int nnod,ndof,nnod_input;
	double longitud,dt_max;
	}dataTube;

typedef struct
	{   
	int nnod,ndof,nnod_input,nunit;
	double Volume, mass, h_film, Area_wall, T_wall;
	}dataTank;

typedef struct
	{   
	//vector<int> type_end;
	int nnod,ndof, nnod_input,modelo_junc;
	//vector<double> state, new_state, Area;
	}dataJunction;
