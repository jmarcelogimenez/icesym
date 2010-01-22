from numpy import array, zeros, sign, arange, equal
cdef extern from "junction.h":

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
	
    #tipo de dato clase junction
	ctypedef struct c_Junction "Junction":
		intvec type_end
		int modelo_junc
		doublevec state_ini
		intvec node2tube
		int implicit
		unsigned int nnod
		unsigned int ndof
		unsigned int nnod_input
		intvec histo
		char* label
	
	c_Junction *new_Junction "new Junction" (unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, doublevec state_ini,intvec histo, char* label, intvec type_end,int modelo_junc, intvec node2tube)
	c_Junction copyJunction "new Junction" (c_Junction* c)
	void del_Junction "delete" (c_Junction *junc)

#defino la clase
cdef class Junction:
	cdef c_Junction *thisptr
    
	def __cinit__(self, **kargs):
    	
		cdef unsigned int nnod 			= validatePositive(kargs,'nnod','Junction')
		cdef unsigned int ndof 			= validatePositive(kargs,'ndof','Junction')
		cdef unsigned int nnod_input 	= validatePositive(kargs,'nnod_input','Junction', nnod)
		kargs['implicit']				= assignOptional(kargs,'implicit',1)
		cdef int implicit 				= boolean(kargs,'implicit','Junction')

    	
		onlyAssert(kargs,'state_ini','Junction')
		cdef doublevec state_ini = doublevec_factory(0)
		for i in range(nnod):
			for j in range(ndof):
				state_ini.push_back(kargs['state_ini'][i][j])
		
		onlyAssert(kargs,'histo','Junction')
		cdef intvec histo = intvec_factory(0)
		for i in range(len(kargs['histo'])):
			if(kargs['histo'][i] in range(nnod)):
				histo.push_back(kargs['histo'][i])
			else:
				print 'Fail inicialitation in [Junction], some node to histo not exists'
				sys.exit()	

		s =  assignOptional(kargs,'label','junc_default')
		cdef char* label =  s
				
		cdef int modelo_junc 			= assignOptional(kargs,'modelo_jun',1)
    	
    	
		# asignacion que puede ser mejorada
    	    			
		cdef intvec type_end = intvec_factory(0)
		for i in range(len(kargs['type_end'])):
			type_end.push_back(kargs['type_end'][i])
    	
		validateSize(kargs,'node2tube','Junction',nnod)
		cdef intvec node2tube = intvec_factory(0)
		for i in range(nnod):
			node2tube.push_back(kargs['node2tube'][i])
			
		self.thisptr = new_Junction(nnod, ndof, nnod_input, implicit,state_ini, histo, label, type_end, modelo_junc, node2tube)
        
        
	def __dealloc__(self):
		del_Junction(self.thisptr)
   
