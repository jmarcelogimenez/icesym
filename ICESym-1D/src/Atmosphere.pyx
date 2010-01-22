from numpy import array, zeros, sign, arange, equal
cdef extern from "atmosphere.h":
    
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
            
    #tipo de dato clase atmosphere
	ctypedef struct c_Atmosphere "Atmosphere":
    
		unsigned int nnod
		unsigned int ndof
		unsigned int nnod_input
		int implicit
		doublevec state_ini
		intvec histo
		char* label
			
	c_Atmosphere *new_Atmosphere "new Atmosphere" (unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, doublevec state_ini, intvec histo, char* label)
	c_Atmosphere copyAtmosphere "new Atmosphere" (c_Atmosphere* a)
	void del_Atmosphere "delete" (c_Atmosphere *atmosphere)

#defino la clase
cdef class Atmosphere:
	cdef c_Atmosphere *thisptr
    
	def __cinit__(self, **kargs):
    	
		cdef unsigned int nnod 			= assignOptional(kargs,'nnod_input',1)
		cdef unsigned int ndof 			= validatePositive(kargs,'ndof','Atmosphere')
		cdef unsigned int nnod_input 	= assignOptional(kargs,'nnod_input',0)
		kargs['implicit']				= assignOptional(kargs,'implicit',0) 
		cdef int implicit 				= boolean(kargs,'implicit','Atmosphere')
    	
		onlyAssert(kargs,'state_ini','Atmosphere')
		cdef doublevec state_ini = doublevec_factory(0)
		#siempre solo 1 nodo?? asi por la interfaz.. se podria cambiar
		#for i in range(nnod):
		#	for j in range(ndof):
		#		state_ini.push_back(kargs['state_ini'][i][j])
		if isinstance(kargs['state_ini'][0],list):
			kargs['state_ini'] = kargs['state_ini'][0]

		for j in range(ndof):
			state_ini.push_back(kargs['state_ini'][j])

		#ATMOSPHERE SIN HISTO!!		
		#onlyAssert(kargs,'histo','Atmosphere')
		cdef intvec histo = intvec_factory(0)
		#for i in range(len(kargs['histo'])):
		#	if(kargs['histo'][i] in range(nnod)):
		#		histo.push_back(kargs['histo'][i])
		#	else:
		#		print 'Fail inicialitation in [%s,%s], some node to histo not exists'
		#		sys.exit()	
		
		s =  assignOptional(kargs,'label','atm_default')
		cdef char* label =  s
				
		self.thisptr = new_Atmosphere(nnod, ndof, nnod_input, implicit, state_ini,histo, label)
        
        
	def __dealloc__(self):
		del_Atmosphere(self.thisptr)
   
