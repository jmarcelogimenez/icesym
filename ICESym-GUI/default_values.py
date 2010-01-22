"""
	This page set the default values for each forms.
"""


#names_simulation[number] -> contiene el nombre de la variable
#default_simulation[number] -> contiene el valor default de la variable

#names_simulation={}
#defaults_simulation={}
dict_simulation = {}

"""
names_simulation[0]='rpms'
names_simulation[1]='ncycles'
names_simulation[2]='nsave'
names_simulation[3]='dtheta_rpm'
names_simulation[4]='Courant'
names_simulation[5]='ga'
names_simulation[6]='viscous_flow'
names_simulation[7]='heat_flow'
names_simulation[8]='R_gas'
names_simulation[9]='nstroke'
names_simulation[10]='engine_type'
names_simulation[11]='ig_order'
names_simulation[12]='filesave_state'
names_simulation[13]='filesave_spd'
names_simulation[14]='get_state'
names_simulation[15]='filein_state'
names_simulation[16]='calc_engine_data'
names_simulation[17]='nappend'
names_simulation[18]='folder_name'

defaults_simulation[0]=[1000,2000]
defaults_simulation[1]=0
defaults_simulation[2]=0
defaults_simulation[3]=0
defaults_simulation[4]=0
defaults_simulation[5]=0
defaults_simulation[6]=0
defaults_simulation[7]=0
defaults_simulation[8]=286.9
defaults_simulation[9]=1
defaults_simulation[10]=0
defaults_simulation[11]=[]
defaults_simulation[12]=''
defaults_simulation[13]=''
defaults_simulation[14]=0
defaults_simulation[15]=''
defaults_simulation[16]=bool(1)
defaults_simulation[17]=10
defaults_simulation[18]="test6"
"""
dict_simulation['rpms']				=[1000,2000]
dict_simulation['ncycles']			= 0
dict_simulation['nsave']			= 0
dict_simulation['dtheta_rpm']		= 0
dict_simulation['Courant']			= 0
dict_simulation['ga']				= 0
dict_simulation['viscous_flow']		= 0
dict_simulation['heat_flow']		= 0
dict_simulation['R_gas']			= 286.9
dict_simulation['nstroke']			= 4
dict_simulation['engine_type']		= 0
dict_simulation['ig_order']			= []
dict_simulation['filesave_state']	= ''
dict_simulation['filesave_spd']		= ''
dict_simulation['get_state']		= 0
dict_simulation['filein_state']		= ''
dict_simulation['calc_engine_data'] = True
dict_simulation['nappend']			= 10
dict_simulation['folder_name']		= "test_default"


"""
names_cylinder = {}
defaults_cylinder = {}

names_cylinder[0]='nnod'
names_cylinder[1]='nvi'
names_cylinder[2]='nve'
names_cylinder[3]='ndof'
names_cylinder[4]='Bore'
names_cylinder[5]='crank_radius'
names_cylinder[6]='Vol_clearance'
names_cylinder[7]='head_chamber_area'
names_cylinder[8]='piston_area'
names_cylinder[9]='delta_ca'
names_cylinder[10]='type_ig'
names_cylinder[11]='twall'
names_cylinder[12]='model_ht'
names_cylinder[13]='factor_ht'
names_cylinder[14]='scavenge'
names_cylinder[15]='scavenge_type'
names_cylinder[16]='SRv'
names_cylinder[17]='full_implicit'
names_cylinder[18]='mass_C'
names_cylinder[19]='data_crevice'
names_cylinder[20]='U_crevice'
names_cylinder[21]='rod_length'
names_cylinder[22]='state_ini'
names_cylinder[23]='ownState'
names_cylinder[24]='extras'
names_cylinder[25]='histo'
names_cylinder[26]='label'
names_cylinder[27]='species_model'
names_cylinder[28]='type_temperature'

defaults_cylinder[0]=1
defaults_cylinder[1]=1
defaults_cylinder[2]=1
defaults_cylinder[3]=3
defaults_cylinder[4]=0
defaults_cylinder[5]=0
defaults_cylinder[6]=0
defaults_cylinder[7]=0
defaults_cylinder[8]=0
defaults_cylinder[9]=0
defaults_cylinder[10]=0
defaults_cylinder[11]=0
defaults_cylinder[12]=0
defaults_cylinder[13]=1.0
defaults_cylinder[14]=bool(0)
defaults_cylinder[15]=0
defaults_cylinder[16]=0
defaults_cylinder[17]=bool(1)
defaults_cylinder[18]=0
defaults_cylinder[19]=0
defaults_cylinder[20]=0
defaults_cylinder[21]=0
defaults_cylinder[22]=[[0,0,0]]
defaults_cylinder[23]=bool(1)
defaults_cylinder[24]=bool(1)
defaults_cylinder[25]=[0]
defaults_cylinder[26]='default'
defaults_cylinder[27]=0
defaults_cylinder[28]=0
"""

dict_cylinder = {}
dict_cylinder['nnod'] = 1
dict_cylinder['nvi'] = 1
dict_cylinder['nve'] = 1
dict_cylinder['ndof'] = 3
dict_cylinder['Bore'] = 0
dict_cylinder['crank_radius'] = 0
dict_cylinder['Vol_clearance'] = 0
dict_cylinder['head_chamber_area'] = 0
dict_cylinder['piston_area'] = 0
dict_cylinder['delta_ca'] = 0
dict_cylinder['type_ig'] = 0
dict_cylinder['twall'] = 0
dict_cylinder['model_ht'] = 0
dict_cylinder['factor_ht'] = 1.0
dict_cylinder['scavenge'] = False
dict_cylinder['scavenge_type'] = 0
dict_cylinder['SRv'] = 0
dict_cylinder['full_implicit'] = True
dict_cylinder['mass_C'] = 0
dict_cylinder['data_crevice'] = 0
dict_cylinder['U_crevice'] = 0
dict_cylinder['rod_length'] = 0
dict_cylinder['state_ini'] = [[0,0,0]]
dict_cylinder['ownState'] = True
dict_cylinder['extras'] = True
dict_cylinder['histo'] = [0]
dict_cylinder['label'] = 'default'
dict_cylinder['species_model'] = 0
dict_cylinder['type_temperature'] = 0
"""
names_fuel={}
defaults_fuel={}

names_fuel[0] = 'Q_fuel'
names_fuel[1] = 'y'
names_fuel[2] = 'hvap_fuel'

defaults_fuel[0] = 44e6
defaults_fuel[1] = 2.25
defaults_fuel[2] = 350000
"""
dict_fuel={}

dict_fuel['Q_fuel'] = 44e6
dict_fuel['y'] = 2.25
dict_fuel['hvap_fuel'] = 350000
"""
names_injection={}
defaults_injection={}

names_injection[0] = 'pulse'
names_injection[1] = 'm_inj'
names_injection[2] = 'dtheta_inj'
names_injection[3] = 'T_fuel'
names_injection[4] = 'theta_inj_ini'
names_injection[5] = 'theta_id'
names_injection[6] = 'ignition_delay_model'
names_injection[7] = 'mfdot_array'

defaults_injection[0] = 0
defaults_injection[1] = 0
defaults_injection[2] = 0
defaults_injection[3] = 0
defaults_injection[4] = 0
defaults_injection[5] = 0.0
defaults_injection[6] = 0
defaults_injection[7] = 0.0
"""
dict_injection={}
dict_injection['pulse'] = 0
dict_injection['m_inj'] = 0
dict_injection['dtheta_inj'] = 0
dict_injection['T_fuel'] = 0
dict_injection['theta_inj_ini'] = 0 
dict_injection['theta_id'] = 0
dict_injection['ignition_delay_model'] = 0
dict_injection['mfdot_array'] = 0

"""
names_combustion={}
defaults_combustion={}

names_combustion[0] = 'theta_ig_0'
names_combustion[1] = 'dtheta_comb'
names_combustion[2] = 'phi'
names_combustion[3] = 'a_wiebe'
names_combustion[4] = 'm_wiebe'
names_combustion[5] = 'combustion_model'
names_combustion[6] = 'xbdot_array'

defaults_combustion[0] = 0
defaults_combustion[1] = 0
defaults_combustion[2] = 0
defaults_combustion[3] = 6.02
defaults_combustion[4] = 1.64
defaults_combustion[5] = 1
defaults_combustion[6] = 0
"""

dict_combustion={}
dict_combustion['theta_ig_0'] = 0
dict_combustion['dtheta_comb'] = 0
dict_combustion['phi'] = 0
dict_combustion['a_wiebe'] = 6.02 
dict_combustion['m_wiebe'] = 1.64
dict_combustion['combustion_model'] = 1
dict_combustion['xbdot_array'] = 0

"""
names_junction={}
defaults_junction={}

names_junction[0] = 'ndof'
names_junction[1] = 'nnod'
names_junction[2] = 'modelo_junc'
#names_junction[3] = 'state_ini'
names_junction[3] = 'histo'
names_junction[4] = 'extras'
names_junction[5] = 'label'

defaults_junction[0] = 3 
defaults_junction[1] = 0
defaults_junction[2] = 0
#defaults_junction[3] = [[1,2],[1,2],[1,2]]
defaults_junction[3] = []
defaults_junction[4] = bool(0)
defaults_junction[5] = 'default'
"""
dict_junction={}
dict_junction['ndof'] = 3 
dict_junction['nnod'] = 0
dict_junction['modelo_junc'] = 0 
dict_junction['histo'] = []
dict_junction['extras'] = False
dict_junction['label'] = 'default_junc'
"""
names_tank={}
defaults_tank={}

names_tank[0] = 'nnod'
names_tank[1] = 'ndof'
names_tank[2] = 'Volume'
names_tank[3] = 'mass'
names_tank[4] = 'h_film'
names_tank[5] = 'Area_wall'
names_tank[6] = 'T_wall'
names_tank[7] = 'Cd_ports'
names_tank[8] = 'state_ini'
names_tank[9] = 'histo'
names_tank[10] = 'extras'
names_tank[11] = 'label'

defaults_tank[0] = 1 
defaults_tank[1] = 3
defaults_tank[2] = 0
defaults_tank[3] = 0.0
defaults_tank[4] = 0.0
defaults_tank[5] = 0
defaults_tank[6] = 0
defaults_tank[7] = 0
defaults_tank[8] = 0
defaults_tank[9] = [1]
defaults_tank[10] = bool(0)
defaults_tank[11] = 'default'
"""
dict_tank={}
dict_tank['nnod'] = 1
dict_tank['ndof'] = 3
dict_tank['Volume'] = 0 
dict_tank['mass'] = 0.0
dict_tank['h_film'] = 0.0
dict_tank['Area_wall'] = 0
dict_tank['T_wall'] = 0
dict_tank['Cd_ports'] = 0 
dict_tank['state_ini'] = 0
dict_tank['histo'] = [1]
dict_tank['extras'] = False
dict_tank['label'] = 'default_tank'
"""
names_tube={}
defaults_tube={}

names_tube[0] = 'nnod'
names_tube[1] = 'ndof'
names_tube[2] = 'longitud'
names_tube[3] = 'xnod'
names_tube[4] = 'diameter'
names_tube[5] = 'twall'
names_tube[6] = 'curvature'
names_tube[7] = 'state_ini'
names_tube[8] = 'typeSave'
names_tube[9] = 'numNorm'
names_tube[10] = 'posNorm'
names_tube[11] = 'histo'
names_tube[12] = 'label'

defaults_tube[0] = 2 
defaults_tube[1] = 3
defaults_tube[2] = 1
defaults_tube[3] = [0,1]
defaults_tube[4] = 0 
defaults_tube[5] = 0
defaults_tube[6] = 0
defaults_tube[7] = 0
defaults_tube[8] = 0
defaults_tube[9] = 2
defaults_tube[10] = [0,1]
defaults_tube[11] = [0,1]
defaults_tube[12] = 'default'
"""
dict_tube={}
dict_tube['nnod'] = 2
dict_tube['ndof'] = 3
dict_tube['longitud'] = 1 
dict_tube['xnod'] = [0,1]
dict_tube['diameter'] = 0
dict_tube['twall'] = 0
dict_tube['curvature'] = 0
dict_tube['state_ini'] = 0
dict_tube['typeSave'] = 0
dict_tube['numNorm'] = 2
dict_tube['posNorm'] = [0,1]
dict_tube['histo'] = [0,1]
dict_tube['label'] = 'default_tube'
"""
names_atmosphere={}
defaults_atmosphere={}

names_atmosphere[0] = 'rho'
names_atmosphere[1] = 'u'
names_atmosphere[2] = 'p'

defaults_atmosphere[0] = 0
defaults_atmosphere[1] = 0
defaults_atmosphere[2] = 0
"""
dict_atmosphere={}
dict_atmosphere['rho'] = 1
dict_atmosphere['u'] = 1
dict_atmosphere['p'] = 1e5

"""
names_valve={}
defaults_valve={}

names_valve[0] = 'Nval'
names_valve[1] = 'type_dat'
names_valve[2] = 'angle_V0'
names_valve[3] = 'angle_VC'
names_valve[4] = 'Dv'
names_valve[5] = 'valve_model'
names_valve[6] = 'type'
names_valve[7] = 'Lvmax'
names_valve[8] = 'Lv'
names_valve[9] = 'Cd'
#names_valve[10] = 'state_ini'
names_valve[10] = 'histo'
names_valve[11] = 'label'

defaults_valve[0] = 1
defaults_valve[1] = 0
defaults_valve[2] = 0
defaults_valve[3] = 0
defaults_valve[4] = 0
defaults_valve[5] = 0
defaults_valve[6] = 0
defaults_valve[7] = 0
defaults_valve[8] = 0
defaults_valve[9] = [[0.0, 1.0],[10,1.0]]
#defaults_valve[10] = [[0,0,0]]
defaults_valve[10] = 0
defaults_valve[11] = 'default'
"""
dict_valve={}
dict_valve['Nval'] = 1
dict_valve['type_dat'] = 0 
dict_valve['angle_V0'] = 0
dict_valve['angle_VC'] = 0
dict_valve['Dv'] = 0
dict_valve['valve_model'] = 0
dict_valve['type'] = 0
dict_valve['Lvmax'] = 0
dict_valve['Lv'] = 0
dict_valve['Cd'] = [[0.0, 1.0],[10,1.0]]
dict_valve['histo'] = 0
dict_valve['label'] = 'default_tank'


giveMeNames = dict()
giveMeNames['Tubes'] = "tube"
giveMeNames['Atmospheres'] = "atmosphere"
giveMeNames['Tanks'] = "tank"
giveMeNames['Valves'] = "valve"
giveMeNames['Cylinders'] = "cylinder"
giveMeNames['Junctions'] = "junction"
giveMeNames['<none>'] = "<none>"

name2element = dict()
name2element['tube'] = "Tubes"
name2element['tank'] = "Tanks"
name2element['junction'] = "Junctions"
name2element['atmosphere'] = "Atmospheres"

