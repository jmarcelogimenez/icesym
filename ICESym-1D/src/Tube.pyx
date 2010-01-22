from numpy import array, zeros, sign, arange, equal
from validator import *
cdef extern from "tube.h":

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
            
    #tipo de dato clase tube
	ctypedef struct c_Tube "Tube":
		unsigned int nnod
		unsigned int ndof
		unsigned int nnod_input
		double dt_max
		int implicit
		doublevec state_ini
		intvec histo
		char* label
		double longitud
		doublevec xnod
		doublevec hele
		doublevec Area
		doublevec twall
		doublevec dAreax
		doublevec curvature
		char* tleft
		unsigned int nleft
		char* tright
		unsigned int nright
       	
	c_Tube *new_Tube "new Tube" (unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, doublevec state_ini, 
				     intvec histo, char* label, double longitud, doublevec xnod, doublevec Area, doublevec twall, 
				     doublevec curvature, doublevec dAreax, char* tleft, unsigned int nleft, char* tright, 
				     unsigned int nright)
	void del_Tube "delete" (c_Tube *tube)
	c_Tube copyTube "new Tube" (c_Tube* t)
	
#defino la clase
cdef class Tube:
	cdef c_Tube *thisptr
    
	def __cinit__(self, **kargs):
    	
		cdef unsigned int nnod 			= validatePositive(kargs,'nnod','Tube')
		cdef unsigned int ndof 			= validatePositive(kargs,'ndof','Tube')
		cdef unsigned int nnod_input 	= validatePositive(kargs,'nnod_input','Tube',2)
		kargs['implicit']				= assignOptional(kargs,'implicit',0)
		cdef int implicit 			= boolean(kargs,'implicit','Tube')
			
		onlyAssert(kargs,'state_ini','Tube')
		cdef doublevec state_ini = doublevec_factory(0)

		#print "nnod: ",nnod, "  ndof:",ndof, " - arreglo: ",len(kargs['state_ini'])," x ",len(kargs['state_ini'][0])
		for i in range(nnod):
				for j in range(ndof):
					state_ini.push_back(kargs['state_ini'][i][j])

		onlyAssert(kargs,'histo','Tube')
		cdef intvec histo = intvec_factory(0)
		for i in range(len(kargs['histo'])):
			if(kargs['histo'][i] in range(nnod)):
				histo.push_back(kargs['histo'][i])
			else:
				print 'Fail inicialitation in [Tube], some node to histo not exists'
				sys.exit()

		s =  assignOptional(kargs,'label','tube_default')
		cdef char* label =  s
				
		cdef double longitud = validatePositive(kargs,'longitud','Tube')
		
		validateSize(kargs,'xnod','Tube',nnod)
		cdef doublevec xnod = doublevec_factory(0)
		for i in range(nnod):
			xnod.push_back(kargs['xnod'][i])
    	
		if not('Area' in kargs.keys()):
			if not('diameter' in kargs.keys()):
				print 'Fail inicialitation in Tube, area and diameter not exists'
				sys.exit()
			else:
				validateSize(kargs,'diameter','Tube',nnod)
				kargs['Area'] = []
				for d in kargs['diameter']:
					kargs['Area'].append(3.14159*(d/2)**2)
		validateSize(kargs,'Area','Tube',nnod)
		cdef doublevec Area = doublevec_factory(0)
		for i in range(nnod):
			Area.push_back(kargs['Area'][i])
			
		validateSize(kargs,'twall','Tube',nnod)
		cdef doublevec twall = doublevec_factory(0)
		for i in range(nnod):
			twall.push_back(kargs['twall'][i])
			
		cdef doublevec dAreax = doublevec_factory(0)
		kargs['dAreax'] = assignOptional(kargs, 'dAreax', []);
		if(len(kargs['dAreax'])>0):
			validateSize(kargs,'dAreax','Tube',nnod)
			for i in range(nnod):
				dAreax.push_back(kargs['dAreax'][i])
			
		if (not(hasattr(kargs,'curvature'))):
			kargs['curvature'] = zeros(nnod)
		cdef doublevec curvature = doublevec_factory(0)
		for i in range(nnod):
			curvature.push_back(kargs['curvature'][i])
		
		s = onlyAssert(kargs,'tleft','Tubes')
		cdef char* tleft = s
		s = onlyAssert(kargs,'tright','Tubes')
		cdef char* tright = s
		cdef unsigned int nleft 	= onlyAssert(kargs,'nleft','Tube')
		cdef unsigned int nright  	= onlyAssert(kargs,'nright','Tube')
		
		self.thisptr = new_Tube(nnod, ndof, nnod_input, implicit, state_ini,histo, label, longitud, 
					xnod, Area, twall, curvature, dAreax, tleft, nleft, 
					tright, nright)
        
	def __dealloc__(self):
		del_Tube(self.thisptr)
   
