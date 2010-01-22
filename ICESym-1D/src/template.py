# Plantilla para inicializar la simulacion

# ----------------------------------------------------
#
# Datos Globales Simulacion
#

Simulator = dict()
Simulator['dt']				=	# double - > 0 - requerido
Simulator['tf']				=	# double - > 0 - opcional: 1.0
Simulator['rpms'] 			=  	# list [rpm1, rpm2,...,rpmN] - requerido
Simulator['iter_sim1d'] 	= 	# int - >0 - opcional: 1
Simulator['nsave'] 			=	# int - >0 - opcional: 0
Simulator['nappend'] 		=	# int - >0 - opcional: 0
Simulator['get_state']		=	# int - requerido, 0: desde aca, 1:desde archivo, 2:desde fortran
Simulator['dtheta_rpm']		=   # double - > 0 - requerido
Simulator['inicia'] 		=	# true = 1, false = 0 - opcional: 0
Simulator['Courant']		=	# double - > 0 - requerido
Simulator['ga']				=	# double - > 0 - requerido
Simulator['viscous_flow'] 	=	# int - > 0 - requerido
Simulator['heat_flow'] 		=	# int - > 0 - requerido
Simulator['R_gas'] 			=	# double - > 0 - opcional: 286.9
Simulator['nstroke'] 		= 	# int - >0 - opcional: 4
Simulator['ncycles'] 		= 	# int - >0 - requerido
Simulator['engine_type']	= 	# opcional int - 0: alternative (default), 1: "Opposed-Piston"
Simulator['filein_state']	= 	# string - opcional: ''
Simulator['filesave_state']	= 	# string - opcional: ''
Simulator['filein_spd']		= 	# string - opcional: ''
Simulator['filesave_spd']	= 	# string - opcional: ''
Simulator['ig_order']		= 	# list - requerido
Simulator['calc_engine_data']= 	# booleano - opcional (defecto True)
Simulator['folder_name']	= 	# string - opcional "testDefault"

# ----------------------------------------------------
#
# Datos Cylinder (completar uno por cada cilindro). Los nodos se ordenaran [nodes_cylinder, valves_intake, valves_exhaust]
#

Cylinders = []


Cyl1 = dict()

Cyl['nnod'] 			=	# deber ser la suma de nodos internos, mas uno por valvula
Cyl['ndof'] 			=	# int - >0 - requerido
Cyl['implicit']			=	# int - 0: false - 1:true
Cyl['state_ini']		=	# ndarray [nnod x ndof] 
Cyl['histo']			= 	# lista de los nodos a guardar historial
Cyl['Bore']				=	# double  - >0 - requerido
Cyl['crank_radius']		=	# double  - >0 - requerido
Cyl['Vol_clearance']	=	# double  - >0 - requerido
Cyl['rod_length']		=	# double  - >0 - requerido
Cyl['head_chamber_area']=	# double  - >0 - requerido
Cyl['piston_area']		=	# double  - >0 - opcional: pi*Bore²
Cyl['delta_ca']			=	# double  - opcional: 0
Cyl['type_ig']			=	# int - opcional: 0
Cyl['twall']			=	# double - >0 - requerido
Cyl['prop']				=	# ndarray [n]
Cyl['U_crevice']		=	# ndarray [n]
Cyl['data_crevice']		=	# ndarray [n]
Cyl['mass_C']			=	# ndarray [n]
Cyl['model_ht']			=	# int - opcional: 0
Cyl['factor_ht']		=	# double - opcional: 1.0	
Cyl['scavenge']			=	# int - 0: false - 1:true
Cyl['scavenge_type']	=	# string - opcional: 'uniflow'
Cyl['SRv']				=	# double  - >0 - requerido
Cyl['full_implicit']	=	# int - 0: false - 1:true, default: 1
Cyl['extras']			=	# int - 0: false - 1:true, default: 0
Cyl['species_model']	=	# int - 0: "<none>" - 1:"single component", default: 0
#datos para compatibilizar con GUI
Cyl['nvi']			= #nro valvulas de admision
Cyl['nve']			= #nro valvulas de escape
#fuel
Cyl['fuel'] = dict()
Cyl['fuel']['Q_fuel']	=	
Cyl['fuel']['y']		=
Cyl['fuel']['hvap_fuel']=

#array de valves, para ser compatible con interfaz grafica
Valves = []

#array de valves de admision
Cyl['intake_valves'] = []

Val = dict()
Val['Nval']			=
Val['type_dat']			=
Val['angle_V0']			=
Val['angle_VC']			=
Val['Dv']			=
Val['Lv']			=		#lista de pares [angulo,valor] ejp: [[0,0],[0.5,1],[1,2]]
Val['Lvmax']			=
Val['valve_model']		=		#int - opcional - 0: "Toyota" (default),  1: "Alessandri"
Val['Cd']			=		#lista de pares [angulo,valor] ejp: [[0,0],[0.5,1],[1,2]]
Val['tube']			=		#numero del tube conectado
Cyl['intake_valves'].append(Val) 
Valves.append(Val)
#array de valves de escape
Cyl['exhaust_valves'] = []

Val = dict()
Val['Nval']			=
Val['type_dat']			=
Val['angle_V0']			=
Val['angle_VC']			=
Val['Dv']			=
Val['Lv']			=		#lista de pares [angulo,valor] ejp: [[0,0],[0.5,1],[1,2]]
Val['Lvmax']			=
Val['valve_model']		=		#int - opcional - 0: "Toyota" (default),  1: "Alessandri"
Val['Cd']			=		#lista de pares [angulo,valor] ejp: [[0,0],[0.5,1],[1,2]]
Val['tube']			=		#numero del tube conectado
Valves.append(Val)
Cyl['exhaust_valves'].append(Val) 		

#injection (solo si diesel, si naftero dejar solo la declaracion del dict())
Cyl['injection'] = dict()
Cyl['injection']['pulse']			=
Cyl['injection']['m_inj']			=
Cyl['injection']['dtheta_inj']			=
Cyl['injection']['T_fuel']			=
Cyl['injection']['theta_inj_ini']		=
Cyl['injection']['theta_id']			=
Cyl['injection']['ignition_delay_model']	=	#0: "L-W" (default), 1: "H-H", 2: "user-defined"
Cyl['injection']['integral']			=
Cyl['injection']['mfdot_array']			=

#combustion
Cyl['combustion']['theta_ig_0']			=
Cyl['combustion']['dtheta_comb']		=
Cyl['combustion']['phi']			=
Cyl['combustion']['phi_ig']			=
Cyl['combustion']['a_wiebe']			=
Cyl['combustion']['m_wiebe']			=
Cyl['combustion']['combustion_model']		=	# 0: "user-defined", 1: "Wiebe-2" (default), 2: "Wiebe-3", 3: "Watson", 4: "<none>"
Cyl['combustion']['start_comb']			=
Cyl['combustion']['xbdot_array']		=

Cylinders.append(Cyl1)


# ----------------------------------------------------
#
# Datos Junction (completar uno por cada unión), no hay un orden predefinido para los nodos, en node2tube se debe dar esa relacion
#

Junctions = []

Junc = dict()
Junc['nnod']				=
Junc['ndof']				=
Junc['implicit']			=
Junc['state_ini']			=
Junc['histo']				=
Junc['modelo_junc']			=
Junc['type_end']			=			# ej [1, -1, 1] 1:conectado a admision, -1: conectado a escape
Junc['node2tube']			= 			# ej [0 2 3] conectado a los tubes 0, 2 y 3, en ese orden!
Junctions.append(Junc)

# ----------------------------------------------------
#
# Datos Tube (completar uno por cada tubo), los nodos se consideran de izquierda a derecha
#

Tubes = []

Tube = dict()
Tube['nnod']				=
Tube['ndof']				=
Tube['implicit']			=
Tube['state_ini']			=
Tube['histo']				=
Tube['longitud']			=
Tube['xnod']				=
Tube['Area']				=
Tube['diameter']			=		# solo sirve si no esta definido "Area"
Tube['twall']				=
Tube['dAreax']				=		# opcional, si no se pasa, se calcula por dif finitas
Tube['curvature']			=
Tube['tleft']				=		# 'junction', 'cylinder', 'tube', 'atmosphere', 'tank'
Tube['nleft']				=		#  posicion en el arreglo 
Tube['tright']				=		# 'junction', 'cylinder', 'tube', 'atmosphere', 'tank'
Tube['nright']				=		#  posicion en el arreglo

Tubes.append(Tube)

# ----------------------------------------------------
#
# Datos Tank (completar uno por cada tanque), los nodos se ordenan asi: [intake_nodes, exhaust_nodes]
#

Tanks = []

Tank = dict()
Tank['nnod']				=
Tank['ndof']				=
Tank['implicit']			=
Tank['state_ini']			=
Tank['histo']				=
Tank['Volume']				=
Tank['mass']				=
Tank['h_film']				=
Tank['Area_wall']			=
Tank['T_wall']				=
Tank['type_end']			=
Tank['Cd_ports']			=
Tank['int2tube']			= 		# lista de indices de los tubos de admision
Tank['exh2tube']			=		# lista de indices de los tubos de escape
Tank['extras']			=	# int - 0: false - 1:true, default: 0
			
Tanks.append(Tank)

# ----------------------------------------------------
#
# Datos Atmosphere (completar uno por cada atmosfera)
#

Atmospheres = []

Atmosphere = dict()
Atmosphere['nnod']			=		#opcional: 1
Atmosphere['ndof']			=		#[rho_atm, vel_atm, p_atm]
Atmosphere['implicit']			=		
Atmosphere['state_ini']			=
Atmosphere['histo']			=

Atmospheres.append(Atmosphere)

# ----------------------------------------------------
#
# Empaquetamiento de los datos
#

kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}
