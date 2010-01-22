from numpy import array, zeros, sign, arange, equal
cdef extern from "tank.h":
    
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
            
    #tipo de dato clase tank
	ctypedef struct c_Tank "Tank":
    
		unsigned int nnod
		unsigned int ndof
		unsigned int nnod_input
		int implicit
		doublevec state_ini
		intvec histo
		char* label
		double Volume
		double mass
		double h_film
		double Area_wall
		double T_wall
		intvec type_end
		doublevec Cd_ports
		intvec int2tube
		intvec exh2tube
		int extras
	
	c_Tank *new_Tank "new Tank" (unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, doublevec state_ini, intvec histo, char* label,
		double Volume, double mass, double h_film, double Area_wall, double T_wall, intvec type_end, doublevec Cd_ports, intvec int2tube, intvec exh2tube, int extras)
	c_Tank copyTank "new Tank" (c_Tank* t)
	void del_Tank "delete" (c_Tank *tank)

#defino la clase
cdef class Tank:
	cdef c_Tank *thisptr
    
	def __cinit__(self, **kargs):
    	
		cdef unsigned int nnod 	     = validatePositive(kargs,'nnod','Tank')
		cdef unsigned int ndof 	     = validatePositive(kargs,'ndof','Tank')
		cdef unsigned int nnod_input = validatePositive(kargs,'nnod_input','Tank',nnod-1)
		kargs['implicit']			= assignOptional(kargs,'implicit',1)
		cdef int implicit 	     	= boolean(kargs,'implicit','Tank')
    	
		onlyAssert(kargs,'state_ini','Tank')
		cdef doublevec state_ini = doublevec_factory(0)
		for i in range(nnod):
			for j in range(ndof):
				state_ini.push_back(kargs['state_ini'][i][j])
				
		onlyAssert(kargs,'histo','Tank')
		cdef intvec histo = intvec_factory(0)
		for i in range(len(kargs['histo'])):
			if(kargs['histo'][i] in range(nnod)):
				histo.push_back(kargs['histo'][i])
			else:
				print 'Fail inicialitation in [Tank], some node to histo not exists'
				sys.exit()
		
		s =  assignOptional(kargs,'label','tank_default')
		cdef char* label =  s

		cdef double Volume = validatePositive(kargs,'Volume','Tank')
		cdef double mass   = validatePositive(kargs,'mass','Tank',0.0) #inicio masa en cero??
		cdef double h_film = assignOptional(kargs,'h_film',0.0)
    	
		cdef double Area_wall = 0
		cdef double T_wall    = 0
		if(h_film>0.0):
			Area_wall = validatePositive(kargs,'Area_wall','Tank')
			T_wall    = validatePositive(kargs,'T_wall','Tank')
    		
		# asignacion que puede ser mejorada, la mejoré??
		cdef intvec type_end = intvec_factory(0)
		#for i in range(len(kargs['type_end'])):
		#	type_end.push_back(kargs['type_end'][i])
		for i in range(len(kargs['int2tube'])):
			type_end.push_back(1)
		for i in range(len(kargs['exh2tube'])):
			type_end.push_back(-1)

		cdef doublevec Cd_ports = doublevec_factory(0)
		for i in range(len(kargs['Cd_ports'])):
			Cd_ports.push_back(kargs['Cd_ports'][i])
		
		onlyAssert(kargs,'int2tube','Tank')
		onlyAssert(kargs,'exh2tube','Tank')
		cdef int l = len(kargs['int2tube'])+len(kargs['exh2tube'])
		
		cdef intvec int2tube = intvec_factory(0)
		cdef intvec exh2tube = intvec_factory(0)
		
		for i in range(len(kargs['int2tube'])):
			int2tube.push_back(kargs['int2tube'][i])
		
		for i in range(len(kargs['exh2tube'])):
			exh2tube.push_back(kargs['exh2tube'][i])
				
		cdef int extras	=  boolean(kargs,'extras','Tank',0)

		self.thisptr = new_Tank(nnod, ndof, nnod_input, implicit, state_ini, histo, label, Volume, mass, h_film, Area_wall, T_wall, type_end, Cd_ports, int2tube, exh2tube, extras)
        
	def __dealloc__(self):
		del_Tank(self.thisptr)
   
