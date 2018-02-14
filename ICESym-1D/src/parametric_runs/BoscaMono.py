#	MODELO MONOCILINDRICO DE MOTOR - BOSCAROL	#20/06
#	Escuela de Ingenieria Mecanica - FCEIA - UNR 


##Parametros geometricos
L0 = 0.14;	D0I = 0.24;	D0F = 0.186;
L1 = 0.26;
D1 = 0.08;
L2a = 0.15;	
L2b = 0.225;
D2aI = 0.05;	D2aF = 0.043;	D2bF = 0.042;
L3 = 0.92;
D3I = 0.044;	D3F = 0.044;
L4 = 0.5;
D4I = 0.060;	D4F = 0.078;	##Antes eran los 2 de 0.65

##Diametros valvulas: 
DValInt = 0.03; DValExh = 0.023; #diametros minimos de asiento segun planos

##Valores empleados para inicializacion:
rho_ini = 1.1647;	p_ini = 101330.0;	T_ini = 330.15;	Q = 0.001; 	#IC basadas en Q para mantener la condicion de continuidad en los canhos
pi = 3.14159;
PHI = 1.4318;
rpms = [8000];
#inicio combustion
theta0 = 5.65;
#duracion combustion
dtheta = 1.6755;
AWiebe = 5;
MWiebe = 2;

#Cd del plenum, para prueba arandela de restriccion
CdPlen = 1;

CdI0 = 0.65; #0.9;
CdE0 = 0.62; #0.9;

CdI = [[0., 0.], [1e-8, CdI0],[0.004, 0.65],[0.012, 0.4], [0.015, 0.44]]
CdE = [[0., 0.], [1e-8, CdE0], [0.007, 0.6], [0.015, 0.4]]

#--------------------Fin definicion parametros del modelo----------------------##
##Librerias empleadas para manejo de arreglos
from numpy import ones, arange, zeros, concatenate

#--------- Inicializacion de Simulator
Simulator0 = dict()							#Se crea el diccionario 'simulator'
Simulator0['filein_state'] = 'mono.dat'	#Archivo para cargar el estado inicial.
Simulator0['filesave_spd'] = ''				#Nombre de archivo para guardar la concentracion quimica de las especies en cada paso de tiempo.
Simulator0['filesave_state'] = 'mono.dat'	#Nombre de archivo para guardar el estado en cada paso del tiempo.
Simulator0['folder_name'] = 'BoscaMono'		#Nombre de la carpeta para guardar la simulacion. Se debe buscarla para realizar el post-proceso.
Simulator0['calc_engine_data'] = 1
Simulator0['heat_flow'] = 1.0				#Calcular transferencia de calor. 0 = no, 1 = si.
Simulator0['viscous_flow'] = 1.0			#Calcular friccion. 0 = no, 1 = si.
Simulator0['get_state'] = 2					#Select how you initialize the state: \nFrom Here-> You must indicates the State in each component \nFrom File-> You must indicates the file to load the state for all components. (Used when you need to continue a simulation) \nFrom Fortran-> ICESym-1D calculates the initial state numerically with a predictor.
Simulator0['nappend'] = 50	#100.0				#Cantidad de iteraciones para guardar la informacion de los componentes. Para ver correctamente los resultados en el post-proceso conviene disminuir este valor. Para mas velocidad, aumentarlo. 0 = nunca.
Simulator0['nsave'] = 10					#Cantidad de iteraciones para guardar el estado de la simulacion. 0 = Nunca.
Simulator0['Courant'] = 0.7
Simulator0['dtheta_rpm'] = 0.05				#Variacion del angulo de ciguenal en cada paso de tiempo.
Simulator0['dt'] = 0.366666666e-7	#0.3666666666667e-7	#
Simulator0['rpms'] = rpms 		#Lista de RPM a simular.
Simulator0['ncycles'] = 6					#Numero de ciclos (termodinamicos) a simular en cada RPM.
Simulator0['nstroke'] = 4					#Tiempos por ciclo (4 tiempos, o 2 tiempos)
Simulator0['R_gas'] = 286.9					#Constante de los gases.
Simulator0['ga'] = 1.36						#Relacion de calores especificos.
Simulator0['ga_intake'] = 1.37
Simulator0['ig_order'] = [0]			#Orden de encendido.
Simulator0['engine_type'] = 0
Simulator0['use_global_gas_properties'] = 1
Simulator = Simulator0
#--------- FIN Inicializacion de Simulator


#--------- Inicializacion de Cylinders
Cylinders = []
#			CILINDRO 0
Cylinders0 = dict()					#Crea el Diccionario 'Cylinders0'.
Cylinders0['label'] = 'CILINDRO0'	#Etiqueta.
Cylinders0['position'] = (965,1100)	#Ubicacion en la pantalla de la GUI.
Cylinders0['nvi'] = 1				#Numero de valvulas de admision por cilindro.
Cylinders0['nve'] = 1				#Numero de valvulas de escape por cilindro.
Cylinders0['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Cylinders0['type_ig'] = 0			#Tipo de Encendido: 0 = por chispa, 1 = por compresion.
Cylinders0['full_implicit'] = 0.0	
Cylinders0['model_ht'] = 3			#Modelo de calculo de transferencia de calor. 3 = Taylor.
Cylinders0['factor_ht'] = 0.8
Cylinders0['ownState'] = 1
Cylinders0['nnod'] = 3				#Numero de nodos. El cilindro tiene 3.
Cylinders0['species_model'] = 0		
Cylinders0['scavenge'] = 0.0
Cylinders0['extras'] = 1
Cylinders0['histo'] = [0, 1, 2]		#Lista de nodos para los cuales se guardan datos.
Cylinders0['delta_ca'] = 0.0		#Solo aplica a motores con pistones opuestos.
Cylinders0['exhaust_valves'] = []
Cylinders0['intake_valves'] = []
Cylinders0['crank_radius'] = 0.044	#Radio de Ciguenal.
Cylinders0['rod_length'] = 0.14		#Longitud de biela.
Cylinders0['piston_area'] = 0.0057	#Superficie de la cabeza de piston.
Cylinders0['Bore'] = 0.0865			#Diametro de cilindro.
Cylinders0['twall'] = [600.0]		#Temperatura de pared del cilindro.
Cylinders0['SRv'] = 0.01			#Scavenge ratio by volume.
Cylinders0['head_chamber_area'] = 0.005741 	#Superficie de pared de la camara de combustion.
Cylinders0['Vol_clearance'] = 5.944e-05		#Volumen de la camara de combustion.
Cylinders0['mass_C'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #Masa de combustible, aire y gases residuales en el cilindro.
Cylinders0['state_ini'] = [[rho_ini, p_ini, T_ini], [rho_ini, 1.0, p_ini], [rho_ini, 1.0, p_ini]]
#--------- Inicializacion de fuel
fuel0 = dict()					#Crea el diccionario 'Fuel0'.
fuel0['y'] = 2.25				#Relacion Hidrogeno/Carbono del Combustible.
fuel0['hvap_fuel'] = 0.01 		#Calor de vaporizacion del combustible [J/kg].
fuel0['Q_fuel'] = 	46000000.0	#Poder calorifico inferior del combustible [J/kg].
Cylinders0['fuel'] = fuel0
#--------- FIN Inicializacion de fuel
#--------- Inicializacion de combustion
combustion0 = dict()			#Crea el diccionario 'Combustion0'.
combustion0['phi'] = PHI 	#Relacion Combustible/Aire.
combustion0['combustion_model'] = 1 #Modelo de combustion.
combustion0['a_wiebe'] = AWiebe		#Parametro de eficiencia de combustion para la funcion Wiebe.
combustion0['m_wiebe'] = MWiebe		#Parametro de forma para la funcion Wiebe.
combustion0['dtheta_comb'] = dtheta #Duracion de la combustion en angulo de ciguenal, en radianes
combustion0['theta_ig_0'] = theta0		#Angulo de inicio de combustion, en radianes. 2*pi= PMS de ciclo de potencia
Cylinders0['combustion'] = combustion0
#--------- FIN Inicializacion de combustion
#--------- Inicializacion de injection
injection0 = dict()	#Crea el diccionario 'Injection0'.
Cylinders0['injection'] = injection0
#--------- FIN Inicializacion de injection
Cylinders.append(Cylinders0) #Agrega el Cilindro 0 al  conjunto de cilindros (sin esta linea el Cilindro 0 no existe).

#--------- Inicializacion de Valves
Valves = []
from numpy import loadtxt
LvI = loadtxt("Lv_I-mb.dat");
LvE = loadtxt("Lv_E-mb.dat")
IVO = 11.3446;	#11.8159
IVC = 4.93928;	#4.41568
EVO = 7.6969;	#8.22050
EVC = 1.3089;	#0.82030
Lvmaxi = max(LvI[:,1]);		Lvmaxe = max(LvE[:,1])
#Lvmaxi = 0.01419;	Lvmaxe = 0.01419

#			VALVULA 0
Valves0 = dict()				#Crea el diccionario 'Valves0'.
Valves0['label'] = 'VALVULA0'	#Etiqueta
Valves0['position'] = (860,1170)#Ubicacion en la pantalla de la GUI.
Valves0['tube'] = 2				#Tubo conectado a la valvula.
Valves0['ncyl'] = 0				#Cilindro conectado a la valvula.
Valves0['typeVal'] = 'int'		#Tipo de valvula: 'int' = admision, 'exh' = escape.
Valves0['type'] = 0				#Tipo de valvula: 0 = admision, 1 = escape.
Valves0['Nval'] = 2				#Numero de valvulas.
Valves0['histo'] = 1			#Lista de nodos para los cuales se guardan datos.
Valves0['valve_model'] = 0  #Modelo de valvula, 0 o 1
Valves0['type_dat'] = 1			#Tipo de Perfil de Leva: 0 = seno cuadrado, 1 = definido por usuario.
Valves0['Lv'] =  LvI			#Arreglos de pares [Angulo(deg),Altura(m)].
Valves0['Lvmax'] = Lvmaxi		#Maxima altura de apertura.
Valves0['angle_V0'] = IVO		#Angulo Apertura.
Valves0['angle_VC'] = IVC		#Angulo Cierre.
Valves0['Dv'] = DValInt			#Diametro cabeza de valvula.
Valves0['Cd'] = CdI				#Coeficiente de descarga.
Cylinders0['intake_valves'].append(Valves0)
Valves.append(Valves0)
#			VALVULA 1
Valves1 = dict()				#Crea el diccionario 'Valves1'.
Valves1['label'] = 'VALVULA1'	#Etiqueta
Valves1['position'] = (1065,1170)#Ubicacion en la pantalla de la GUI.
Valves1['tube'] = 3				#Tubo conectado a la valvula.
Valves1['ncyl'] = 0				#Cilindro conectado a la valvula.
Valves1['typeVal'] = 'exh'		#Tipo de valvula: 'int' = admision, 'exh' = escape.
Valves1['type'] = 1				#Tipo de valvula: 0 = admision, 1 = escape.
Valves1['Nval'] = 2				#Numero de valvulas.
Valves1['histo'] = 1			#Lista de nodos para los cuales se guardan datos.
Valves1['valve_model'] = 0  #Modelo de valvula, 0 o 1
Valves1['type_dat'] = 1			#Tipo de Perfil de Leva: 0 = seno cuadrado, 1 = definido por usuario.
Valves1['Lv'] = LvE				#Arreglos de pares [Angulo(deg),Altura(m)].
Valves1['Lvmax'] = Lvmaxe		#Maxima altura de apertura.
Valves1['angle_V0'] = EVO		#Angulo Apertura.
Valves1['angle_VC'] = EVC		#Angulo Cierre.
Valves1['Dv'] = DValExh			#Diametro cabeza de valvula.
Valves1['Cd'] = CdE				#Coeficiente de descarga. #Coeficiente de descarga.
Cylinders0['exhaust_valves'].append(Valves1)
Valves.append(Valves1)

#--------- Inicializacion de Tubes
Tubes = []
#			Notas respecto a los tubos: Un nodo cada 5 mm.
#			Se conectan solo a atmosferas, juntas, tanques y cilindros.

#			TUBO 0
N0 = 28;
diameter0 = D0I*ones(N0+1)+(D0F-D0I)/N0*arange(N0+1)
StateIni0 = zeros((N0+1,3))
StateIni0[:,0] = rho_ini
StateIni0[:,1] = Q/(pi*(diameter0/2)*(diameter0/2))
StateIni0[:,2] = p_ini
Tubes0 = dict() #Crea el diccionario 'Tubes0'.
Tubes0['label'] = 'TUBO0'		#Etiqueta.
Tubes0['position'] = (200,625)	#Ubicacion en la pantalla de la GUI.
Tubes0['typeSave'] = 1
Tubes0['tleft'] = 'atmosphere'	#Tipo de elemento a la izquierda.
Tubes0['nleft'] = 0				#Indice del elemento de la izquierda.
Tubes0['tright'] = 'tank'		#Tipo de elemento a la derecha.
Tubes0['nright'] = 0			#Indice del elemento de la derecha.
Tubes0['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tubes0['nnod'] = N0+1			#Cantidad de nodos.
Tubes0['numNorm'] = 2			#Se deja 2.
Tubes0['histo'] = [0, N0] 		#Lista de nodos para los cuales se guardan datos.
Tubes0['posNorm'] = [0.0, 1.0]
Tubes0['longitud'] = L0			#Longitud del tubo.
Tubes0['diameter'] = diameter0.copy()
Tubes0['twall'] = T_ini*ones(N0+1)
Tubes0['xnod'] = L0/N0*arange(N0+1)
Tubes0['state_ini'] = StateIni0.copy() #[Densidad, Velocidad, Presion] x Cantidad de Nodos.
Tubes0['type'] = 'intake'
Tubes.append(Tubes0)

#			TUBO 1
Tubes1 = dict() #Crea el diccionario 'Tubes1'.
N1 = int(L1/0.005);
diameter1 = D1*ones(N1+1)
StateIni1 = zeros((N1+1,3))
StateIni1[:,0] = rho_ini
StateIni1[:,1] = Q/(pi*(diameter1/2)*(diameter1/2))
StateIni1[:,2] = p_ini
Tubes1['label'] = 'TUBO1'		#Etiqueta.
Tubes1['position'] = (435,625)	#Ubicacion en la pantalla de la GUI.
Tubes1['typeSave'] = 1
Tubes1['tleft'] = 'tank'		#Tipo de elemento a la izquierda.
Tubes1['nleft'] = 0				#Indice del elemento de la izquierda.
Tubes1['tright'] = 'tank'		#Tipo de elemento a la derecha.
Tubes1['nright'] = 1			#Indice del elemento de la derecha.
Tubes1['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tubes1['nnod'] = N1+1			#Cantidad de nodos.
Tubes1['numNorm'] = 2			#Se deja 2.
Tubes1['histo'] = [0, N1] 		#Lista de nodos para los cuales se guardan datos.
Tubes1['posNorm'] = [0.0, 1.0]
Tubes1['longitud'] = L1
Tubes1['diameter'] = diameter1.copy()
Tubes1['twall'] = T_ini*ones(N1+1)
Tubes1['xnod'] = L1/N1*arange(N1+1)
Tubes1['state_ini'] = StateIni1.copy() #[Densidad, Velocidad, Presion] x Cantidad de Nodos.
Tubes1['type'] = 'intake'
Tubes.append(Tubes1)

#			TUBO 2

N2a = int(L2a/0.005);
N2b = int(L2b/0.005);
L2 = L2a+L2b;
N2 = N2a+N2b+1
diameter2 = D2aI*ones(N2a+1)+(D2aF-D2aI)/N2a*arange(N2a+1)
diameter2 = concatenate((diameter2,D2aF*ones(N2b+1)+(D2bF-D2aF)/N2b*arange(N2b+1)),axis=0)
StateIni2 = zeros((N2+1,3))
StateIni2[:,0] = rho_ini
StateIni2[:,1] = Q/(pi*(diameter2/2)*(diameter2/2))
StateIni2[:,2] = p_ini
Tubes2 = dict() #Crea el diccionario 'Tubes2'.
Tubes2['label'] = 'TUBO2'		#Etiqueta.
Tubes2['position'] = (730,1170)	#Ubicacion en la pantalla de la GUI.
Tubes2['typeSave'] = 1
Tubes2['tleft'] = 'tank'	#Tipo de elemento a la izquierda.
Tubes2['nleft'] = 1				#Indice del elemento de la izquierda.
Tubes2['tright'] = 'cylinder'	#Tipo de elemento a la derecha.
Tubes2['nright'] = 0			#Indice del elemento de la derecha.
Tubes2['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tubes2['nnod'] = N2+1			#Cantidad de nodos.
Tubes2['numNorm'] = 2			#Se deja 2.
Tubes2['histo'] = [0, N2] 		#Lista de nodos para los cuales se guardan datos.
Tubes2['posNorm'] = [0.0, 1.0]
Tubes2['longitud'] = L2			#Longitud del tubo.
Tubes2['diameter'] = diameter2.copy()
Tubes2['twall'] = T_ini*ones(N2+1)
Tubes2['xnod'] = L2/N2*arange(N2+1)
Tubes2['state_ini'] = StateIni2.copy() #[Densidad, Velocidad, Presion] x Cantidad de Nodos.
Tubes2['type'] = 'intake'
Tubes.append(Tubes2)

#			TUBO 3
N3 = int(L3/0.005)
diameter3 = D3I*ones(N3+1)+(D3F-D3I)/N3*arange(N3+1)
StateIni3 = zeros((N3+1,3))
StateIni3[:,0] = rho_ini
StateIni3[:,1] = Q/(pi*(diameter3/2)*(diameter3/2))
StateIni3[:,2] = p_ini
Tubes3 = dict() #Crea el diccionario 'Tubes3'.
Tubes3['label'] = 'TUBO3'		#Etiqueta.
Tubes3['position'] = (730,880)	#Ubicacion en la pantalla de la GUI.
Tubes3['typeSave'] = 1
Tubes3['tleft'] = 'cylinder'	#Tipo de elemento a la izquierda.
Tubes3['nleft'] = 0				#Indice del elemento de la izquierda.
Tubes3['tright'] = 'tank'	#Tipo de elemento a la derecha.
Tubes3['nright'] = 2			#Indice del elemento de la derecha.
Tubes3['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tubes3['nnod'] = N3+1			#Cantidad de nodos.
Tubes3['numNorm'] = 2			#Se deja 2.
Tubes3['histo'] = [0, N3] 		#Lista de nodos para los cuales se guardan datos.
Tubes3['posNorm'] = [0.0, 1.0]
Tubes3['longitud'] = L3			#Longitud del tubo.
Tubes3['diameter'] = diameter3.copy()
Tubes3['twall'] = T_ini*ones(N3+1)
Tubes3['xnod'] = L3/N3*arange(N3+1)
Tubes3['state_ini'] = StateIni3.copy() #[Densidad, Velocidad, Presion] x Cantidad de Nodos.
Tubes3['type'] = 'exhaust'
Tubes.append(Tubes3)

#			TUBO 4
N4 = int(L4/0.005)
diameter4 = D4I*ones(N4+1)+(D4F-D4I)/N4*arange(N4+1)
StateIni4 = zeros((N4+1,3))
StateIni4[:,0] = rho_ini
StateIni4[:,1] = Q/(pi*(diameter4/2)*(diameter4/2))
StateIni4[:,2] = p_ini
Tubes4 = dict() #Crea el diccionario 'Tubes4'.
Tubes4['label'] = 'TUBO4'		#Etiqueta.
Tubes4['position'] = (730,480)	#Ubicacion en la pantalla de la GUI.
Tubes4['typeSave'] = 1
Tubes4['tleft'] = 'tank'	#Tipo de elemento a la izquierda.
Tubes4['nleft'] = 2				#Indice del elemento de la izquierda.
Tubes4['tright'] = 'atmosphere'	#Tipo de elemento a la derecha.
Tubes4['nright'] = 1			#Indice del elemento de la derecha.
Tubes4['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tubes4['nnod'] = N4+1			#Cantidad de nodos.
Tubes4['numNorm'] = 2			#Se deja 2.
Tubes4['histo'] = [0, N4] 		#Lista de nodos para los cuales se guardan datos.
Tubes4['posNorm'] = [0.0, 1.0]
Tubes4['longitud'] = L4			#Longitud del tubo.
Tubes4['diameter'] = diameter4.copy()
Tubes4['twall'] = T_ini*ones(N4+1)
Tubes4['xnod'] = L4/N4*arange(N4+1)
Tubes4['state_ini'] = StateIni4.copy() #[Densidad, Velocidad, Presion] x Cantidad de Nodos.
Tubes4 ['type'] = 'exhaust'
Tubes.append(Tubes4)
#--------- Fin Inicializacion de Tubos


#--------- Inicializacion de Tubes
Junctions = []	#en este modelo no hay juntas
#
#--------- Fin Inicializacion de Tubos

#--------- Inicializacion de Tanks
Tanks = []
#			TANQUE 0
Tanks0 = dict()				#Crea el diccionario 'Tanks0'.
Tanks0['label'] = 'TANQUE0'	#Etiqueta.
Tanks0['position'] = (315,625)	#Ubicacion en la pantalla de la GUI.
Tanks0['implicit']  = 1
Tanks0['int2tube'] = [0]		#Tubos que Entran.
Tanks0['exh2tube'] = [1]		#Tubos que Salen.
Tanks0['type_end']  = [1, -1]
Tanks0['nnod'] = 3				#Cantidad de nodos. nnod = tubos+1
Tanks0['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tanks0['extras'] = 0			
Tanks0['histo'] = [0, 1, 2]		#Lista de nodos para los cuales se guardan datos.
Tanks0['Area_wall'] = 0.175	#Superficie de pared.
Tanks0['Volume'] = 0.0075		#Volumen del tanque.
Tanks0['mass'] = 0.01		#Masa de gas contenida en el tanque.
Tanks0['Cd_ports'] = [1.0, 1.0]	#Coeficiente de descarga. Para juntas y tanques: 1.0 #0.144 referencia bounous para la toma de aire
Tanks0['T_wall'] = 330.0		#Temperatura de pared.
Tanks0['h_film'] = 300.0		#Coeficiente de transferencia de calor [W/m^2 K].
Tanks0['state_ini'] = [[rho_ini, p_ini, T_ini], [rho_ini, 1.0, p_ini], [rho_ini, 1.0, p_ini]]
Tanks0['type'] = 'intake'
Tanks.append(Tanks0)

#			TANQUE 1
Tanks1 = dict()					#Crea el diccionario 'Tanks1'.
Tanks1['label'] = 'TANQUE1'		#Etiqueta.
Tanks1['position'] = (565,625)	#Ubicacion en la pantalla de la GUI.
Tanks1['implicit']  = 1
Tanks1['int2tube'] = [1]		#Tubos que Entran.
Tanks1['exh2tube'] = [2]		#Tubos que Salen.
Tanks1['type_end']  = [1, -1]
Tanks1['nnod'] = 3				#Cantidad de nodos. nnod = tubos+1
Tanks1['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tanks1['extras'] = 0			
Tanks1['histo'] = [0, 1, 2]	#Lista de nodos para los cuales se guardan datos.
Tanks1['Area_wall'] = 0.025	#Superficie de pared.
Tanks1['Volume'] = 0.00037504		#Volumen del tanque.
Tanks1['mass'] = 0.00045		#Masa de gas contenida en el tanque.
Tanks1['Cd_ports'] = [1.0, CdPlen]	#Coeficiente de descarga. Para juntas y tanques: 1.0. Para el tubo 1 se contempla la restriccion de la brida
Tanks1['T_wall'] = 330.0		#Temperatura de pared.
Tanks1['h_film'] = 300.0		#Coeficiente de transferencia de calor [W/m^2 K].
Tanks1['state_ini'] = [[rho_ini, p_ini, T_ini], [rho_ini, 1.0, p_ini], [rho_ini, 1.0, p_ini]]
Tanks1['type'] = 'intake'
Tanks.append(Tanks1)

#			TANQUE 2
Tanks2 = dict()					#Crea el diccionario 'Tanks2'.
Tanks2['label'] = 'TANQUE2'		#Etiqueta.
Tanks2['position'] = (1900,625)	#Ubicacion en la pantalla de la GUI.
Tanks2['implicit']  = 1
Tanks2['int2tube'] = [3]		#Tubos que Entran.
Tanks2['exh2tube'] = [4]		#Tubos que Salen.
Tanks2['type_end']  = [1, -1]
Tanks2['nnod'] = 3				#Cantidad de nodos. nnod = tubos+1
Tanks2['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Tanks2['extras'] = 0			
Tanks2['histo'] = [0, 1, 2]		#Lista de nodos para los cuales se guardan datos.
Tanks2['Area_wall'] = 0.061261	#Superficie de pared.
Tanks2['Volume'] = 0.00183775		#Volumen del tanque.
Tanks2['mass'] = 0.0021425		#Masa de gas contenida en el tanque.
Tanks2['Cd_ports'] = [1.0, 1.0]	#Coeficiente de descarga. Para juntas y tanques: 1.0
Tanks2['T_wall'] = 330.0		#Temperatura de pared.
Tanks2['h_film'] = 300.0		#Coeficiente de transferencia de calor [W/m^2 K].
Tanks2['state_ini'] = [[rho_ini, p_ini, T_ini], [rho_ini, 1.0, p_ini], [rho_ini, 1.0, p_ini]]
Tanks2['type'] = 'exhaust'
Tanks.append(Tanks2)
#--------- FIN Inicializacion de Tanks

#--------- Inicializacion de Atmospheres
Atmospheres = []

#			ATMOSFERA 0
Atmospheres0 = dict()					#Crea el diccionario 'Atmosphere0'.
Atmospheres0['position'] = (170,625)	#Ubicacion en la pantalla de la GUI.
Atmospheres0['implicit']  = 0
Atmospheres0['nnod'] = 1				#Cantidad de nodos. nnod = 1.
Atmospheres0['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Atmospheres0['state_ini'] = [rho_ini, 2.2222, p_ini]	#[Densidad, Velocidad, Presion Abs].
Atmospheres.append(Atmospheres0)

#			ATMOSFERA 1
Atmospheres1 = dict()					#Crea el diccionario 'Atmosphere1'.
Atmospheres1['position'] = (2150,625)	#Ubicacion en la pantalla de la GUI.
Atmospheres1['implicit']  = 0
Atmospheres1['nnod'] = 1				#Cantidad de nodos. nnod = 1.
Atmospheres1['ndof'] = 3				#Numero de grados de libertad: Siempre ndof=3.
Atmospheres1['state_ini'] = [rho_ini, 4.666, p_ini]	#[Densidad, Velocidad, Presion Abs].
Atmospheres.append(Atmospheres1)
#--------- FIN Inicializacion de Atmospheres


#UBICA LOS DATOS
kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}
