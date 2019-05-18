#### ---- ####
# Archivo generado por SimulatorGUI
# CIMEC - Santa Fe - Argentina 
# Adecuado para levantar desde Interfaz Grafica 
# O para correr desde consola mediante $python main.py 
#### ---- ####

Simulator0 = dict()
Simulator0['dtheta_rpm'] = 1.0
Simulator0['filein_state'] = 'None'
Simulator0['Courant'] = 0.8
Simulator0['heat_flow'] = 1.0
Simulator0['R_gas'] = 287.0
Simulator0['rpms'] = [1000, 1400, 1800, 2200]
Simulator0['calc_engine_data'] = 1
Simulator0['filesave_state'] = 'KamAZ7405'
Simulator0['ncycles'] = 5
Simulator0['nsave'] = 5
Simulator0['folder_name'] = 'KamAZ7405'
Simulator0['ga'] = 1.4
Simulator0['viscous_flow'] = 1.0
Simulator0['filesave_spd'] = ''
Simulator0['ig_order'] = [0]
Simulator0['get_state'] = 2
Simulator0['nappend'] = 5.0
Simulator0['engine_type'] = 0
Simulator0['nstroke'] = 4

Simulator = Simulator0


#--------- Inicializacion de Cylinders

Cylinders = []

Cylinders0 = dict()
Cylinders0['crank_radius'] = 0.06
Cylinders0['type_ig'] = 1
Cylinders0['ndof'] = 3
Cylinders0['full_implicit'] = 0.0
Cylinders0['model_ht'] = 1
Cylinders0['factor_ht'] = 0.36
Cylinders0['piston_area'] = 0.011309734
Cylinders0['exhaust_valves'] = []
Cylinders0['ownState'] = 1
Cylinders0['mass_C'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Cylinders0['nnod'] = 3
Cylinders0['label'] = 'cyl'
Cylinders0['twall'] = [459.0]
Cylinders0['state_ini'] = [[1.1769, 101330.0, 300.0], [1.3195, 0.1, 128000.0], [1.1769, 0.1, 101330.0]]
Cylinders0['nve'] = 1
Cylinders0['intake_valves'] = []
Cylinders0['head_chamber_area'] = 0.012817698
Cylinders0['type_temperature'] = 0
Cylinders0['rod_length'] = 0.225
Cylinders0['species_model'] = 0
Cylinders0['nvi'] = 1
Cylinders0['delta_ca'] = 0.0
Cylinders0['Vol_clearance'] = 9.04778684234e-05
Cylinders0['extras'] = 1
Cylinders0['histo'] = [0]
Cylinders0['position'] = (529, 199)
Cylinders0['Bore'] = 0.12
Cylinders0['scavenge'] = 0.0

#--------- Inicializacion de combustion

combustion0 = dict()
combustion0['m_wiebe'] = 1.64
combustion0['combustion_model'] = 1
combustion0['theta_ig_0'] = 6.1261056745
combustion0['a_wiebe'] = 6.0
combustion0['dtheta_comb'] = 1.20427718388
Cylinders0['combustion'] = combustion0


#--------- FIN Inicializacion de combustion


#--------- Inicializacion de fuel

fuel0 = dict()
fuel0['y'] = 2.16667
fuel0['hvap_fuel'] = 350000.0
fuel0['Q_fuel'] = 42500000.0
Cylinders0['fuel'] = fuel0


#--------- FIN Inicializacion de fuel


#--------- Inicializacion de injection

injection0 = dict()
injection0['ignition_delay_model'] = 1
injection0['dtheta_inj'] = 0.130009575981
injection0['m_inj'] = 0.0001218
injection0['theta_inj_ini'] = 6.0388392119
injection0['T_fuel'] = 300.0
injection0['pulse'] = 1
Cylinders0['injection'] = injection0


#--------- FIN Inicializacion de injection


Cylinders.append(Cylinders0)


#--------- FIN Inicializacion de Cylinders


#--------- Inicializacion de Valves

Valves = []

Valves0 = dict()
Valves0['Lvmax'] = 0.0088452
Valves0['tube'] = 0
Valves0['angle_VC'] = 4.01425727959
Valves0['ncyl'] = 0
Valves0['label'] = 'ivalve'
Valves0['histo'] = 0
Valves0['Nval'] = 1
Valves0['Dv'] = 0.04
Valves0['position'] = (397, 199)
Valves0['type_dat'] = 0
Valves0['Cd'] = [[0.0, 1.0], [0.009, 1.0]]
Valves0['type'] = 0
Valves0['typeVal'] = 'int'
Valves0['angle_V0'] = 12.3045712266
Valves0['valve_model'] = 1
Cylinders0['intake_valves'].append(Valves0)


Valves.append(Valves0)

Valves1 = dict()
Valves1['Lvmax'] = 0.00880308
Valves1['tube'] = 1
Valves1['angle_VC'] = 0.261799387799
Valves1['ncyl'] = 0
Valves1['label'] = 'evalve'
Valves1['histo'] = 0
Valves1['Nval'] = 1
Valves1['Dv'] = 0.04
Valves1['position'] = (661, 199)
Valves1['type_dat'] = 0
Valves1['Cd'] = [[0.0, 1.0], [0.009, 1.0]]
Valves1['type'] = 1
Valves1['typeVal'] = 'exh'
Valves1['angle_V0'] = 8.29031394697
Valves1['valve_model'] = 1
Cylinders0['exhaust_valves'].append(Valves1)


Valves.append(Valves1)


#--------- FIN Inicializacion de Valves


#--------- Inicializacion de Tubes

Tubes = []

Tubes0 = dict()
Tubes0['diameter'] = [0.044, 0.0442666666667, 0.0445333333333, 0.0448, 0.0450666666667, 0.0453333333333, 0.0456, 0.0458666666667, 0.0461333333333, 0.0464, 0.0466666666667, 0.0469333333333, 0.0472, 0.0474666666667, 0.0477333333333, 0.048, 0.0482666666667, 0.0485333333333, 0.0488, 0.0490666666667, 0.0493333333333, 0.0496, 0.0498666666667, 0.0501333333333, 0.0504, 0.0506666666667, 0.0509333333333, 0.0512, 0.0514666666667, 0.0517333333333, 0.052]
Tubes0['longitud'] = 0.15
Tubes0['nnod'] = 31
Tubes0['numNorm'] = 3
Tubes0['twall'] = [292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0, 292.0]
Tubes0['ndof'] = 3
Tubes0['state_ini'] = [[1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0], [1.3195, 0.1, 128000.0]]
Tubes0['tright'] = 'cylinder'
Tubes0['label'] = 'tube0'
Tubes0['histo'] = [0, 15, 30]
Tubes0['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15]
Tubes0['nleft'] = 0
Tubes0['position'] = (265, 199)
Tubes0['tleft'] = 'atmosphere'
Tubes0['posNorm'] = [0.0, 0.5, 1.0]
Tubes0['nright'] = 0
Tubes0['typeSave'] = 1

Tubes.append(Tubes0)

Tubes1 = dict()
Tubes1['diameter'] = [0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038]
Tubes1['longitud'] = 0.4
Tubes1['nnod'] = 81
Tubes1['numNorm'] = 3
Tubes1['twall'] = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
Tubes1['ndof'] = 3
Tubes1['state_ini'] = [[1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tubes1['tright'] = 'atmosphere'
Tubes1['label'] = 'tube1'
Tubes1['histo'] = [0, 40, 80]
Tubes1['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 0.36, 0.365, 0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 0.4]
Tubes1['nleft'] = 0
Tubes1['position'] = (793, 199)
Tubes1['tleft'] = 'cylinder'
Tubes1['posNorm'] = [0.0, 0.5, 1.0]
Tubes1['nright'] = 1
Tubes1['typeSave'] = 1

Tubes.append(Tubes1)


#--------- FIN Inicializacion de Tubes


#--------- Inicializacion de Tanks

Tanks = []


#--------- FIN Inicializacion de Tanks


#--------- Inicializacion de Junctions

Junctions = []


#--------- FIN Inicializacion de Junctions


#--------- Inicializacion de Atmospheres

Atmospheres = []

Atmospheres0 = dict()
Atmospheres0['nnod'] = 1
Atmospheres0['ndof'] = 3
Atmospheres0['state_ini'] = [1.3195, 0.1, 128000.0]
Atmospheres0['position'] = (133, 199)

Atmospheres.append(Atmospheres0)

Atmospheres1 = dict()
Atmospheres1['nnod'] = 1
Atmospheres1['ndof'] = 3
Atmospheres1['state_ini'] = [1.4351, 0.1, 120330.0]
Atmospheres1['position'] = (925, 199)

Atmospheres.append(Atmospheres1)


#--------- FIN Inicializacion de Atmospheres

kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}