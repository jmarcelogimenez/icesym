#### ---- ####
# Archivo generado por SimulatorGUI
# CIMEC - Santa Fe - Argentina 
# Adecuado para levantar desde Interfaz Grafica 
# O para correr desde consola mediante $python main.py 
#### ---- ####

from numpy import add, alltrue, arange, array, ones, zeros
from math import pi, sqrt

#
# Bicilindrico con un tanque en la admision y escape 2 a 1:
#   Tube-Tank-2*Tube-2*Intake_valve-2*Cylinder-2*Exhaust_valve-2*Tube-Junction-Tube
#

# Plantilla para inicializar la simulacion

# ----------------------------------------------------
#
# Datos Globales Simulacion
#

rho_ref = 1.
p_ref   = 1e5
u_ref   = 1.0

T_ref = 1e5/287

Simulator = dict()
Simulator['dt']			= 1e-7
Simulator['tf']			= 10
Simulator['rpms'] 		= [1000,1200,1400,1600,1800,2000]
Simulator['dtheta_rpm']		= 1.0
Simulator['inicia'] 		= 1
Simulator['Courant']		= 0.8
Simulator['ga']			= 1.4
Simulator['viscous_flow'] 	= 1
Simulator['heat_flow'] 		= 0
Simulator['R_gas'] 		= 287.0
Simulator['nstroke'] 		= 4
Simulator['ncycles'] 		= 4
Simulator['engine_type']	= 0
Simulator['ig_order']		= [0,1]
Simulator['nsave']			= 20
Simulator['nappend']		= 50
Simulator['get_state']		= 2
Simulator['folder_name']	= 'test6'
Simulator['filein_state']	= 'test6_state.dat'
Simulator['filesave_state']	= 'test6_state.dat'
Simulator['calc_engine_data'] = 1

# ----------------------------------------------------
#
# Datos Cylinder (completar uno por cada cilindro). Los nodos se ordenaran [nodes_cylinder, valves_intake, valves_exhaust]
#

Cylinders = []

#array de valves, para ser compatible con interfaz grafica
Valves = []

Cyl0 = dict()
Cyl0['nnod'] 			= 3
Cyl0['ndof'] 			= 3
Cyl0['implicit']		= 1
Cyl0['state_ini']		= zeros((3,3))
Cyl0['histo']			= [0,1,2]
Cyl0['extras']			= 1
Cyl0['label']			= "Cyl0"
Cyl0['Bore']			= 0.086
Cyl0['crank_radius']		= 0.034
Cyl0['Vol_clearance']		= 4.4e-5
Cyl0['rod_length']		= 0.118
Cyl0['head_chamber_area']	= 0.008
Cyl0['piston_area']		= 0.006
Cyl0['theta_0']			= 0.0
Cyl0['delta_ca']		= 0.0
Cyl0['type_ig']			= 0
Cyl0['twall']			= [600.0]
Cyl0['prop']			= [1]
Cyl0['U_crevice']		= [1]
Cyl0['data_crevice']		= [1]
Cyl0['mass_C']			= zeros(6)
Cyl0['model_ht']		= 3
Cyl0['factor_ht']		= 0.6
Cyl0['scavenge']		= 0
Cyl0['scavenge_type']		= 'uniflow'
Cyl0['SRv']			= 0.01
Cyl0['full_implicit']		= 0
Cyl0['nvi']						= 1
Cyl0['nve']						= 1
#fuel
Cyl0['fuel'] = dict()
Cyl0['fuel']['Q_fuel']		= 46.0e6
Cyl0['fuel']['y']		= 2.25
Cyl0['fuel']['hvap_fuel']	= 0.01

#array de valves de admision
Cyl0['intake_valves'] = []

Val0 = dict()
Val0['Nval']			= 1
Val0['type_dat']		= 0
Val0['angle_V0']		= 12.1126
Val0['angle_VC']		= 4.0317
Val0['Dv']			= 0.0303
Val0['Lvmax']			= 0.009
Val0['valve_model']		= 0
Val0['Cd']			= [[0.0,0.0], [0.00151,0.13904794], [0.00303,0.26121148], [0.00454,0.36053144],[0.00606,0.45389219 ],[0.00757,0.48666778],[0.00909,0.51050457],[0.0106,0.52308509]]
Val0['tube']			= 1		#numero del tube conectado
Cyl0['intake_valves'].append(Val0)
Valves.append(Val0)
#array de valves de escape
Cyl0['exhaust_valves'] = []

Val1 = dict()
Val1['Nval']			= 1
Val1['type_dat']		= 0
Val1['angle_V0']		= 8.587
Val1['angle_VC']		= 0.4363
Val1['Dv']			= 0.03
Val1['Lvmax']			= 0.0085
Val1['valve_model']		= 0
Val1['Cd']			= [[0.0, 0.0], [0.0015, 0.14082964], [0.003, 0.2691636], [0.0045, 0.36507635], [0.006, 0.412695], [0.0075, 0.43262054], [0.009, 0.44140128], [0.0105, 0.44545394]]
Val1['tube']			= 3		#numero del tube conectado
Cyl0['exhaust_valves'].append(Val1)
Valves.append(Val1)
#injection
Cyl0['injection'] = dict()

#combustion
Cyl0['combustion'] = dict()
Cyl0['combustion']['theta_ig_0']	       	= 5.94
Cyl0['combustion']['dtheta_comb']		= 1.4835
Cyl0['combustion']['phi']			= 1.3904
Cyl0['combustion']['combustion_model']		= 1

Cylinders.append(Cyl0)

Cyl1 = dict()
Cyl1['nnod']                    = 3
Cyl1['ndof']                    = 3
Cyl1['implicit']                = 1
Cyl1['state_ini']               = zeros((3,3))
Cyl1['histo']                   = [0,1,2]
Cyl1['extras']					= 1
Cyl1['label']					= "Cyl1"
Cyl1['Bore']                    = 0.086
Cyl1['crank_radius']            = 0.034
Cyl1['Vol_clearance']           = 4.4e-5
Cyl1['rod_length']              = 0.118
Cyl1['head_chamber_area']       = 0.008
Cyl1['piston_area']             = 0.006
Cyl1['theta_0']                 = 2.*pi
Cyl1['delta_ca']                = 0.0
Cyl1['type_ig']                 = 0
Cyl1['twall']                   = [600.0]
Cyl1['prop']                    = [1]
Cyl1['U_crevice']               = [1]
Cyl1['data_crevice']            = [1]
Cyl1['mass_C']                  = zeros(6)
Cyl1['model_ht']                = 3
Cyl1['factor_ht']               = 0.6
Cyl1['scavenge']                = 0
Cyl1['scavenge_type']           = 'uniflow'
Cyl1['SRv']                     = 0.01
Cyl1['full_implicit']           = 0
Cyl1['nvi']						= 1
Cyl1['nve']						= 1
#fuel
Cyl1['fuel'] = dict()
Cyl1['fuel']['Q_fuel']          = 46.0e6
Cyl1['fuel']['y']               = 2.25
Cyl1['fuel']['hvap_fuel']       = 0.01

#array de valves de admision
Cyl1['intake_valves'] = []

Val2 = dict()
Val2['Nval']                    = 1
Val2['type_dat']                = 0
Val2['angle_V0']                = 12.1126
Val2['angle_VC']                = 4.0317
Val2['Dv']                      = 0.0303
Val2['Lvmax']                   = 0.009
Val2['valve_model']             = 0
Val2['Cd']                      = [[0.0,0.0], [0.00151,0.13904794], [0.00303,0.26121148], [0.00454,0.36053144], [0.00606,0.45389219], [0.00757,0.48666778], [0.00909,0.51050457], [0.0106,0.52308509]] 
Val2['tube']                    = 2             #numero del tube conectado

Cyl1['intake_valves'].append(Val2)
Valves.append(Val2)
#array de valves de escape
Cyl1['exhaust_valves'] = []

Val3 = dict()
Val3['Nval']                    = 1
Val3['type_dat']                = 0
Val3['angle_V0']                = 8.587
Val3['angle_VC']                = 0.4363
Val3['Dv']                      = 0.03
Val3['Lvmax']                   = 0.0085
Val3['valve_model']             = 0
Val3['Cd']                      = [[0.0, 0.0], [0.0015, 0.14082964], [0.003, 0.2691636], [0.0045, 0.36507635], [0.006, 0.412695], [0.0075, 0.43262054], [0.009, 0.44140128], [0.0105, 0.44545394]]
Val3['tube']                    = 4             #numero del tube conectado

Cyl1['exhaust_valves'].append(Val3)
Valves.append(Val3)
#injection
Cyl1['injection'] = dict()

#combustion
Cyl1['combustion'] = dict()
Cyl1['combustion']['theta_ig_0']                = 5.94
Cyl1['combustion']['dtheta_comb']               = 1.4835
Cyl1['combustion']['phi']                       = 1.3904
Cyl1['combustion']['combustion_model']          = 1

Cylinders.append(Cyl1)

# ----------------------------------------------------
#
# Datos Junction (completar uno por cada union), no hay un orden predefinido para los nodos, en node2tube se debe dar esa relacion
#

Junctions = []

Junc0 = dict()
Junc0['nnod']                           = 3
Junc0['ndof']                           = 3
Junc0['implicit']                       = 1
Junc0['state_ini']                      = zeros((3,3))
Junc0['histo']                          = [0,1,2]
Junc0['label']							= "Junc0"
Junc0['modelo_junc']                    = 1
Junc0['type_end']                       = [1,1,-1]
Junc0['node2tube']                      = [3,4,5]

Junctions.append(Junc0)

# ----------------------------------------------------
#
# Datos Tube (completar uno por cada tubo), los nodos se consideran de izquierda a derecha
#

N = 50
state = zeros((N+1,3))
state[:,0] = rho_ref
state[:,1] = u_ref
state[:,2] = p_ref

Tubes = []

Tube0 = dict()
Tube0['nnod']				= N+1
Tube0['ndof']				= 3
Tube0['implicit']			= 0
Tube0['state_ini']			= state.copy()
Tube0['histo']				= [0,N]
Tube0['label']				= "Tube0"
Tube0['longitud']			= .25
Tube0['xnod']				= .25/N*arange(N+1)
Tube0['Area']				= 0.065**2*pi/4*ones(N+1)
Tube0['twall']				= 330.*ones(N+1)
Tube0['curvature']			= zeros(N+1)
Tube0['tleft']				= 'atmosphere'
Tube0['nleft']				= 0		#  posicion en el arreglo
Tube0['tright']				= 'tank'
Tube0['nright']				= 0		#  posicion en el arreglo

Tubes.append(Tube0)

N = 120
state = zeros((N+1,3))
state[:,0] = rho_ref
state[:,1] = u_ref
state[:,2] = p_ref

Tube1 = dict()
Tube1['nnod']				= N+1
Tube1['ndof']				= 3
Tube1['implicit']			= 0
Tube1['state_ini']			= state.copy()
Tube1['histo']				= [0,N]
Tube1['label']				= "tube1"
Tube1['longitud']			= .6
Tube1['xnod']				= .6/N*arange(N+1)
Tube1['Area']				= .0303**2*pi/4*ones(N+1)
Tube1['twall']				= 330.*ones(N+1)
Tube1['curvature']			= zeros(N+1)
Tube1['tleft']				= 'tank'
Tube1['nleft']				= 0
Tube1['tright']				= 'cylinder'
Tube1['nright']				= 0

Tubes.append(Tube1)

Tube2 = dict()
Tube2['nnod']                           = N+1
Tube2['ndof']                           = 3
Tube2['implicit']                       = 0
Tube2['state_ini']                      = state.copy()
Tube2['histo']                          = [0,N]
Tube2['label']                          = "Tube2"
Tube2['longitud']                       = .6
Tube2['xnod']                           = .6/N*arange(N+1)
Tube2['Area']                           = .0303**2*pi/4*ones(N+1)
Tube2['twall']                          = 330.*ones(N+1)
Tube2['curvature']                      = zeros(N+1)
Tube2['tleft']                          = 'tank'
Tube2['nleft']                          = 0
Tube2['tright']                         = 'cylinder'
Tube2['nright']                         = 1

Tubes.append(Tube2)

N = 60
state = zeros((N+1,3))
state[:,0] = rho_ref
state[:,1] = u_ref
state[:,2] = p_ref

Tube3 = dict()
Tube3['nnod']				= N+1
Tube3['ndof']				= 3
Tube3['implicit']			= 0
Tube3['state_ini']			= state.copy()
Tube3['histo']				= [0]
Tube3['label']				= "Tube3"
Tube3['longitud']			= .3
Tube3['xnod']				= .3/N*arange(N+1)
Tube3['Area']				= .03**2*pi/4*ones(N+1)
Tube3['twall']				= 600.*ones(N+1)
Tube3['curvature']			= zeros(N+1)
Tube3['tleft']				= 'cylinder'
Tube3['nleft']				= 0
Tube3['tright']				= 'junction'
Tube3['nright']				= 0

Tubes.append(Tube3)

Tube4 = dict()
Tube4['nnod']                           = N+1
Tube4['ndof']                           = 3
Tube4['implicit']                       = 0
Tube4['state_ini']                      = state.copy()
Tube4['histo']                          = [0,N]
Tube4['label']                          = "Tube4"
Tube4['longitud']                       = .3
Tube4['xnod']                           = .3/N*arange(N+1)
Tube4['Area']                           = .03**2*pi/4*ones(N+1)
Tube4['twall']                          = 600.*ones(N+1)
Tube4['curvature']                      = zeros(N+1)
Tube4['tleft']                          = 'cylinder'
Tube4['nleft']                          = 1
Tube4['tright']                         = 'junction'
Tube4['nright']                         = 0

Tubes.append(Tube4)

N = 100
state = zeros((N+1,3))
state[:,0] = rho_ref
state[:,1] = u_ref
state[:,2] = p_ref

Tube5 = dict()
Tube5['nnod']                           = N+1
Tube5['ndof']                           = 3
Tube5['implicit']                       = 0
Tube5['state_ini']                      = state.copy()
Tube5['histo']                          = [0]
Tube5['label']                          = "Tube5"
Tube5['longitud']                       = .5
Tube5['xnod']                           = .5/N*arange(N+1)
Tube5['Area']                           = .03**2*pi/4*ones(N+1)
Tube5['twall']                          = 400.*ones(N+1)
Tube5['curvature']                      = zeros(N+1)
Tube5['tleft']                          = 'junction'
Tube5['nleft']                          = 0
Tube5['tright']                         = 'atmosphere'
Tube5['nright']                         = 0

Tubes.append(Tube5)

# ----------------------------------------------------
#
# Datos Tank (completar uno por cada tanque), los nodos se ordenan asi: [intake_nodes, exhaust_nodes]
#

Tanks = []

Tank = dict()
Tank['nnod']				= 4
Tank['ndof']				= 3
Tank['implicit']			= 1
Tank['state_ini']			= zeros((4,3))
Tank['Volume']				= 0.0015
Tank['mass']				= 0.0015
Tank['h_film']				= 300
Tank['Area_wall']			= 0.09
Tank['T_wall']				= 330
Tank['type_end']			= [1,-1,-1]
Tank['Cd_ports']			= [.5, .5, .5]
Tank['int2tube']			= [0]		# lista de indices de los tubos de admision
Tank['exh2tube']			= [1,2]		# lista de indices de los tubos de escape
Tank['histo']               = [0,1,2]
Tank['label']               = "Tank0"

Tanks.append(Tank)

# ----------------------------------------------------
#
# Datos Atmosphere (completar uno por cada atmosfera)
#

Atmospheres = []

Atmosphere = dict()
Atmosphere['nnod']			= 1				#opcional: 1
Atmosphere['ndof']			= 3		
Atmosphere['implicit']		= 0		
Atmosphere['state_ini']		= [[1.1647,0.1,101330.0]]		#[rho_atm, vel_atm, p_atm]
Atmosphere['histo']			= []

Atmospheres.append(Atmosphere)

# ----------------------------------------------------
#
# Empaquetamiento de los datos
#

kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}
