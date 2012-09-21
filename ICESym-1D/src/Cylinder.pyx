from numpy import array, zeros, sign, arange, equal
from validator import *
cdef extern from "cylinder.h":

    # templates de la stl
	ctypedef struct intvec "std::vector<int>":
		void (* push_back)(int elem)
		int (* at)(unsigned int index)
		void acc "operator[]"(unsigned int index)
	intvec intvec_factory "std::vector<int>"(int len)
    
	ctypedef struct doublevec "std::vector<double>":
		void (* push_back)(double elem)
		double (* at)(unsigned int index)
		void acc "operator[]"(unsigned int index)
	doublevec doublevec_factory "std::vector<double>"(int len)
    
    # estructuras internas de cylinder
	ctypedef struct fuel:
		double Q_fuel
		double y
		double hvap_fuel
    
	ctypedef struct valve:
		int Nval
		int type_dat
		doublevec Cd
		doublevec Lv
		double angle_V0
		double angle_VC
		double Dv
		double Lvmax
		double* dx_tube
		double* area_tube
		double* twall_tube
		double* dAreax_tube
		int valve_model
		int tube
	
	ctypedef struct Scavenge:
		double val_1
		double val_2
		double SRv
		int close_cyl
	
	ctypedef struct injection:
		int pulse
		double m_inj
		double dtheta_inj
		double T_fuel
		double theta_inj_ini
		double theta_id
		int ignition_delay_model
		double integral
		doublevec mfdot_array

	ctypedef struct combustion:
		double theta_ig_0
		double dtheta_comb
		double phi
		double phi_ig
		double a_wiebe
		double m_wiebe
		int combustion_model
		int start_comb
		doublevec xbdot_array
	            
	ctypedef struct valvevec "std::vector<valve>":
		void (* push_back)(valve elem)
		valve (* at)(unsigned int index)
		valve acc "operator[]"(unsigned int index)
	valvevec valvevec_factory "std::vector<valve>"(int len)
    	            
    #tipo de dato clase cylinder
	ctypedef struct c_Cylinder "Cylinder":    
		unsigned int nnod
		unsigned int ndof
		unsigned int nnod_input
		int implicit
		doublevec state_ini
		intvec histo
		char* label
		double Bore
		double crank_radius
		double Vol_clearance
		double rod_length
		double head_chamber_area
		double piston_area
		double theta_0
		double delta_ca
		doublevec Twall
		doublevec prop
		doublevec U_crevice
		doublevec data_crevice
		doublevec mass_C
		int model_ht
		double factor_ht
		int scavenge
		char* scavenge_type
		int type_ig
		int full_implicit
		int extras
		int species_model
		
		Scavenge scavenge_data
		fuel fuel_data
		injection injection_data
		combustion combustion_data
		valvevec intake_valves
		valvevec exhaust_valves
		
	c_Cylinder *new_Cylinder "new Cylinder" (unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, 
						 doublevec state_ini, intvec histo, char* label, double Bore, double crank_radius, 
						 double Vol_clearance, double rod_length, double head_chamber_area, 
						 double piston_area, double theta_0, double delta_ca, doublevec Twall, doublevec prop, 
						 doublevec U_crevice, doublevec data_crevice, doublevec mass_C, int model_ht, 
						 double factor_ht, int scavenge,char* scavenge_type, int type_ig, int full_implicit, 
						 fuel fuel_data, combustion combustion_data, injection injection_data, 
						 valvevec intake_valves, valvevec exhaust_valves, Scavenge scavenge_data, int extras, int species_model)
	void del_Cylinder "delete" (c_Cylinder *cyl)
	c_Cylinder copyCylinder "new Cylinder" (c_Cylinder* c)
	
#defino la clase
cdef class Cylinder:
	cdef c_Cylinder *thisptr
    
	def __cinit__(self, **kargs):
    	
		cdef unsigned int nnod 			= validatePositive(kargs,'nnod','Cylinder')
		cdef unsigned int ndof 			= validatePositive(kargs,'ndof','Cylinder')
		kargs['implicit']				= assignOptional(kargs,'implicit',1)
		cdef int implicit 				= boolean(kargs,'implicit','Cylinder')
    	
		onlyAssert(kargs,'state_ini','Cylinder')
		cdef doublevec state_ini = doublevec_factory(0)
		for i in range(nnod):
			for j in range(ndof):
				state_ini.push_back(kargs['state_ini'][i][j])
		
		onlyAssert(kargs,'histo','Cylinder')
		cdef intvec histo = intvec_factory(0)
		#print kargs['histo'], range(nnod)
		for i in range(len(kargs['histo'])):
			if(kargs['histo'][i] in range(nnod)):
				histo.push_back(kargs['histo'][i])
			else:
				print 'Fail inicialitation in [Cylinder], some node to histo not exists'
				sys.exit()	
		s =  assignOptional(kargs,'label','cyl_default')
		cdef char* label =  s

		vargs = kargs['intake_valves']
		nintake = len(vargs)
		vargs = vargs + kargs['exhaust_valves']
		 
		cdef double Bore         		= validatePositive(kargs,'Bore','Cylinder')
		cdef double crank_radius 		= validatePositive(kargs,'crank_radius','Cylinder')
		cdef double Vol_clearance		= validatePositive(kargs,'Vol_clearance','Cylinder')
		cdef double rod_length      	= validatePositive(kargs,'rod_length','Cylinder')
		cdef double head_chamber_area 	= validatePositive(kargs,'head_chamber_area','Cylinder')
		cdef double piston_area 		= validatePositive(kargs,'piston_area','Cylinder',(3.1416*Bore**2)/4)
		# ver calculo de theta_0 !!!!!!!!!!!! 
		cdef double theta_0 			= kargs['theta_0']
		cdef double delta_ca	 		= assignOptional(kargs,'delta_ca',0)

		onlyAssert(kargs,'twall','Cylinder')
		cdef doublevec twall	 		= doublevec_factory(0)
		for i in range(len(kargs['twall'])):
			twall.push_back(kargs['twall'][i])
		#cdef double twall	 			= validatePositive(kargs,'twall','Cylinder')

		cdef int type_ig				= assignOptional(kargs,'type_ig',0)
		cdef int full_implicit 			= boolean(kargs,'full_implicit','Cylinder',1)
		cdef doublevec prop 			= doublevec_factory(0)
		#onlyAssert(kargs,'prop','Cylinder')
		kargs['prop'] = assignOptional(kargs,'prop',[1])
		for i in range(len(kargs['prop'])):
			prop.push_back(kargs['prop'][i])
    	
		cdef doublevec U_crevice 		= doublevec_factory(0)
		#onlyAssert(kargs,'U_crevice','Cylinder')
		kargs['U_crevice'] = assignOptional(kargs,'U_crevice',[1])
		for i in range(len(kargs['U_crevice'])):
			U_crevice.push_back(kargs['U_crevice'][i])
	
		cdef doublevec data_crevice 	= doublevec_factory(0)
		#onlyAssert(kargs,'data_crevice','Cylinder')
		kargs['data_crevice'] = assignOptional(kargs,'data_crevice',[1])
		for i in range(len(kargs['data_crevice'])):
			data_crevice.push_back(kargs['data_crevice'][i])
    	
		cdef doublevec mass_C			 = doublevec_factory(0)
		validateSize(kargs,'mass_C','Cylinder',6*(nnod-len(vargs)))
		onlyAssert(kargs,'mass_C','Cylinder')
		for i in range(len(kargs['mass_C'])):
			mass_C.push_back(kargs['mass_C'][i])
    	
		cdef int model_ht				=  assignOptional(kargs,'model_ht',0)
		cdef double factor_ht			=  assignOptional(kargs,'factor_ht',1.0)
		cdef int scavenge				=  assignOptional(kargs,'scavenge',0)
		s								=  assignOptional(kargs,'scavenge_type','uniflow')
		cdef char* scavenge_type		=  s
		cdef int extras					=  boolean(kargs,'extras','Cylinder',0)
		cdef int species_model			=  assignOptional(kargs,'species_model',0)
    	
    	#condiciones para scavenge
		cdef Scavenge scavenge_data
		if(not(kargs['scavenge']==0)):
			#print "tengo scavenge"
			if(kargs['scavenge_type']=='scre'):
				scavenge_data.val_1 = -1.6709
				scavenge_data.val_2 = 0.1899
			else: 
				if(kargs['scavenge_type']=='yam1'):
					scavenge_data.val_1 = -1.6993
					scavenge_data.val_2 = 0.3053
				else:
					if(kargs['scavenge_type']=='yam6'):
						scavenge_data.val_1 = -1.3516
						scavenge_data.val_2 = 0.1435
					else:
						if(kargs['scavenge_type']=='cd'):
							scavenge_data.val_1 = -1.0104
							scavenge_data.val_2 = -0.117
						else: 
							if(kargs['scavenge_type']=='qubcr'):
								scavenge_data.val_1 = -1.6325
								scavenge_data.val_2 = 0.1397
							else :									#'uniflow'
								scavenge_data.val_1 = -1.7827
								scavenge_data.val_2 = 0.2094
		else:
			scavenge_data.val_1 = 0
			scavenge_data.val_2 = 0
		scavenge_data.SRv	=  validatePositive(kargs,'SRv','Cylinder',0)
		scavenge_data.close_cyl =  0
		
		#condiciones para fuel
		fargs						= kargs['fuel']
		cdef fuel fuel_data
		fuel_data.Q_fuel    = validatePositive(fargs,'Q_fuel','Fuel-Cylinder',44e6)
		fuel_data.y	    = validatePositive(fargs,'y','Fuel-Cylinder',2.25)
		fuel_data.hvap_fuel = validatePositive(fargs,'hvap_fuel','Fuel-Cylinder',350000)
    	
		#condiciones para valves
		vargs = kargs['intake_valves']
		nintake = len(vargs)
		vargs = vargs + kargs['exhaust_valves']
		cdef valvevec intake_valves 	= valvevec_factory(0)
		cdef valvevec exhaust_valves 	= valvevec_factory(0)
		cdef valve auxValve
		for i in range(len(vargs)):
			auxValve.tube	     = onlyAssert(vargs[i],'tube','Valve-Cylinder')
			auxValve.Nval	     = validatePositive(vargs[i],'Nval','Valve-Cylinder',1)
			auxValve.type_dat    = onlyAssert(vargs[i],'type_dat','Valve-Cylinder')
			auxValve.angle_V0    = onlyAssert(vargs[i],'angle_V0','Valve-Cylinder')
			auxValve.angle_VC    = onlyAssert(vargs[i],'angle_VC','Valve-Cylinder')
			auxValve.Dv	     = validatePositive(vargs[i],'Dv','Valve-Cylinder')
			auxValve.Lv 	     = doublevec_factory(0)
			auxValve.Lvmax 	     = validatePositive(vargs[i],'Lvmax','Valve-Cylinder',0)
			auxValve.valve_model = validateInList(vargs[i],'valve_model','Valve-Cylinder',[0,1],0)
			if (auxValve.type_dat == 0): 				# ley senoidal cuadrada, a calcular en Fortran
				auxValve.Lv.push_back(-1)
				auxValve.Lv.push_back(-1)
			else:
				if (auxValve.type_dat == -1):			# ley exponencial, a calcular en Fortran
					auxValve.Lv.push_back(-1)
					auxValve.Lv.push_back(-1)
			#	else:   # se debe ingresar un array(2x721) y lo mapea a un array unidimensional [ang,val,ang,val...]
			#		# if (validateSize(vargs[i],'Lv','Valve-Cylinder',721)):
			#		auxValve.Lv = []
			#		for j in range(len(vargs[i]['Lv'][0])):
			#			auxValve.Lv.push_back(vargs[i]['Lv'][0][j])
			#			auxValve.Lv.push_back(vargs[i]['Lv'][1][j])
			# ahora se recibe pares [angulo,valor]
				else:
					auxValve.Lv = []
					for j in range(len(vargs[i]['Lv'])):
						auxValve.Lv.push_back(vargs[i]['Lv'][j][0])
						auxValve.Lv.push_back(vargs[i]['Lv'][j][1])			
				
			auxValve.Cd = []
			#for j in range(len(vargs[i]['Cd'][0])):
			#	auxValve.Cd.push_back(vargs[i]['Cd'][0][j])
			#	auxValve.Cd.push_back(vargs[i]['Cd'][1][j])
			for j in range(len(vargs[i]['Cd'])):
				auxValve.Cd.push_back(vargs[i]['Cd'][j][0])
				auxValve.Cd.push_back(vargs[i]['Cd'][j][1])		
			#print "Cd: ", vargs[i]['Cd']
			
			if(i<nintake):				
				intake_valves.push_back(auxValve)
			else:
				exhaust_valves.push_back(auxValve)
		
		cdef unsigned int nnod_input = validatePositive(kargs,'nnod_input','Cylinder',len(kargs['intake_valves'])+len(kargs['exhaust_valves']))
					
		#condiciones para injection
		cdef injection injection_data
		injection_data.mfdot_array = doublevec_factory(0)
		if(type_ig==1):
			iargs 				    = kargs['injection']
			injection_data.pulse		    = assignOptional(iargs,'pulse',2)
			injection_data.m_inj		    = validatePositive(iargs,'m_inj','Injection-Cylinder')
			injection_data.dtheta_inj	    = validatePositive(iargs,'dtheta_inj','Injection-Cylinder')
			injection_data.T_fuel		    = validatePositive(iargs,'T_fuel','Injection-Cylinder')
			injection_data.theta_inj_ini	    = validatePositive(iargs,'theta_inj_ini','Injection-Cylinder')
			injection_data.theta_id		    = validatePositive(iargs,'theta_id','Injection-Cylinder',0.0)
			injection_data.ignition_delay_model = validateInList(iargs,'ignition_delay_model','Injection-Cylinder',[0,1,2],0)
			injection_data.integral		    = validatePositive(iargs,'integral','Injection-Cylinder',0.0)
			if(injection_data.pulse == 3):
				mfdot = onlyAssert(iargs,'mfdot_array','Injection-Cylinder')
				#for i in range(len(mfdot[0])): #[key1,value1,key2,value2,....,keyN,valueN]
				#	injection_data.mfdot_array.push_back(mfdot[0][i])
				#	injection_data.mfdot_array.push_back(mfdot[1][i])	
				for i in range(len(mfdot)): #[key1,value1,key2,value2,....,keyN,valueN]
					injection_data.mfdot_array.push_back(mfdot[i][0])
					injection_data.mfdot_array.push_back(mfdot[i][1])			

		#condiciones para combustion
		cdef combustion combustion_data
		cargs 				= kargs['combustion']
		combustion_data.xbdot_array	= doublevec_factory(0)
		combustion_data.theta_ig_0	= validatePositive(cargs,'theta_ig_0','Combustion-Cylinder')
		combustion_data.dtheta_comb	= validatePositive(cargs,'dtheta_comb','Combustion-Cylinder')
		
		if(type_ig==0):
			combustion_data.phi	= onlyAssert(cargs,'phi','Combustion-Cylinder')
		else:
			combustion_data.phi	= 0
			
		combustion_data.phi_ig	= assignOptional(cargs,'phi_ig',0.0) #como asignar internamente??
		combustion_data.a_wiebe	= assignOptional(cargs,'a_wiebe',6.02)
		combustion_data.m_wiebe	= assignOptional(cargs,'m_wiebe',1.64)
		
		if(type_ig==1):
			combustion_data.combustion_model = validateInList(cargs,'combustion_model','Combustion-Cylinder',[0,1,2,3,4],1)
		else:
			combustion_data.combustion_model = validateInList(cargs,'combustion_model','Combustion-Cylinder',[0,1,2,3,4],4)
		
		combustion_data.start_comb	= assignOptional(cargs,'start_comb',1)
		if(combustion_data.combustion_model==0):
			xbdot = onlyAssert(cargs,'xbdot_array','Combustion-Cylinder')
			#for i in range(len(xbdot[0])): #[key1,value1,key2,value2,....,keyN,valueN]
			#	combustion_data.xbdot_array.push_back(xbdot[0][i])
			#	combustion_data.xbdot_array.push_back(xbdot[1][i])		
			for i in range(len(xbdot)): #[key1,value1,key2,value2,....,keyN,valueN]
				combustion_data.xbdot_array.push_back(xbdot[i][0])
				combustion_data.xbdot_array.push_back(xbdot[i][1])		
					
		#instancio la clase
		self.thisptr = new_Cylinder(nnod, ndof, nnod_input, implicit, state_ini, histo, label, Bore, crank_radius,
					    Vol_clearance, rod_length, head_chamber_area, piston_area, theta_0, 
					    delta_ca, twall, prop, U_crevice, data_crevice, mass_C, model_ht, 
					    factor_ht, scavenge, scavenge_type, type_ig, full_implicit,fuel_data, 
					    combustion_data, injection_data, intake_valves, exhaust_valves, scavenge_data, extras, species_model)

	def __dealloc__(self):
		del_Cylinder(self.thisptr)
