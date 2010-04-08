from numpy import array, asarray, zeros, sign, arange, equal
from myfunction import *
cdef extern from "simulator.h":

	# templates de la stl
	ctypedef struct intvec "std::vector<int>":
		void (* push_back)(int elem)
		int (* at)(unsigned int index)
		int acc "operator[]"(unsigned int index)
	intvec intvec_factory "std::vector<int>"(int len)
    
	ctypedef struct doublevec "std::vector<double>":
		void (* push_back)(int elem)
		double (* at)(unsigned int index)
		double acc "operator[]"(unsigned int index)
		int (* size)()
	doublevec doublevec_factory "std::vector<double>"(int len)
    
	ctypedef struct cylindervec "std::vector<Cylinder>":
		void (* push_back)(c_Cylinder elem)
		c_Cylinder (* at)(unsigned int index)
		c_Cylinder acc "operator[]"(unsigned int index)
	cylindervec cylindervec_factory "std::vector<Cylinder>"(int len)
    
	ctypedef struct tankvec "std::vector<Tank>":
		void (* push_back)(c_Tank elem)
		c_Tank (* at)(unsigned int index)
		c_Tank acc "operator[]"(unsigned int index)
	tankvec tankvec_factory "std::vector<Tank>"(int len)
    
	ctypedef struct junctionvec "std::vector<Junction>":
		void (* push_back)(c_Junction elem)
		c_Junction (* at)(unsigned int index)
		c_Junction acc "operator[]"(unsigned int index)
	junctionvec junctionvec_factory "std::vector<Junction>"(int len)
    
	ctypedef struct tubevec "std::vector<Tube>":
		void (* push_back)(c_Tube elem)
		c_Tube (* at)(unsigned int index)
		c_Tube acc "operator[]"(unsigned int index)
	tubevec tubevec_factory "std::vector<Tube>"(int len)
	
	ctypedef struct atmvec "std::vector<Atmosphere>":
		void (* push_back)(c_Atmosphere elem)
		c_Atmosphere (* at)(unsigned int index)
		c_Atmosphere acc "operator[]"(unsigned int index)
	atmvec atmvec_factory "std::vector<Atmosphere>"(int len)
	
	#tipo de dato clase simulator
	ctypedef struct c_Simulator "Simulator":
       	
		double dt
		double tf
		int nrpms
		int irpm
		doublevec rpms
		double time 
		doublevec Xn  
		doublevec Xn1
		int ntubes
		int ncyl
		int ntank
		int njunc
		int natm
		int iter_sim1d
		int nsave
		int nappend
		int get_state
		double dtheta_rpm
		int inicia
		double Courant
		double ga
		int viscous_flow
		int heat_flow
		double R_gas
		double theta
		double omega 
		int nstroke
		int ncycles
		int icycle
		int engine_type
		char* filein_state
		char* filesave_state
		char* filein_spd
		char* filesave_spd
		char* folder_name
		intvec ig_order
		
		cylindervec cylinders
		tubevec tubes
		junctionvec junctions
		tankvec tanks
		atmvec atmospheres
       	
		int calc_engine_data
	
		void solver()
		void solveStep()
		void printData()
		void actualizeState()
		void actualizeDt()
		double GetTimeStep()
		void SetTimeStep(double dT)
		void SetTime(double t, int irpm)
		doublevec GetNewState()
		doublevec GetOldState()
		int GetNewStateSize()
		int GetOldStateSize()
		void SetStateValue(int i, double val)

		doublevec getCylinderMass(int icyl)
		void setCylinderMass(int icyl, int i, double mass)
       	
	c_Simulator *new_Simulator "new Simulator" (double dt, double tf, int nrpms, doublevec rpms, doublevec Xn, doublevec Xn1,
						    int ntubes, int ncyl, int ntank, int njunc, int iter_sim1d, int nsave, int nappend,
						    double dtheta_rpm, int inicia, double Courant, double ga, int viscous_flow, 
						    int heat_flow, double R_gas, int nstroke, int ncycles, int engine_type, 
						    char* filein_state, char* filesave_state, char* filein_spd, char* filesave_spd,char* folder_name,
						    intvec ig_order, cylindervec cylinders, tubevec tubes, junctionvec junctions,
						    tankvec tanks, int natm, atmvec atmospheres, int get_state,int calc_engine_data)
	void del_Simulator "delete" (c_Simulator *sim)

#defino la clase
cdef class Simulator:
	cdef c_Simulator *thisptr
    
	def __cinit__(self,**kargs):
	
		dataCylinders 	= kargs['Cylinders']
		dataTubes 	= kargs['Tubes']
		dataJunctions 	= kargs['Junctions']
		dataTanks	= kargs['Tanks']
		dataAtmospheres	= kargs['Atmospheres']
		sargs		= kargs['Simulator']
			
		cdef double dt	= assignOptional(sargs,'dt',0.001)
		cdef double tf	= assignOptional(sargs,'tf',1.0)
		
		cdef int ntubes	= len(dataTubes)
		cdef int ncyl	= len(dataCylinders)
		cdef int ntank	= len(dataTanks)
		cdef int njunc	= len(dataJunctions)
		cdef int natm	= len(dataAtmospheres)
		
		cdef int iter_sim1d = assignOptional(sargs,'iter_sim1d',1)
		cdef int nsave	    = validatePositive(sargs,'nsave','Simulator',0)
		cdef int nappend    = validatePositive(sargs,'nappend','Simulator',0)
		cdef int get_state  = assignOptional(sargs,'get_state',0)
		

		cdef doublevec Xn = doublevec_factory(0)
		# ahora hecho en c++!!!
		# el estado global se lee de un archivo
		#if(get_state==1):
		#	state_file = onlyAssert(sargs,'filein_state','Simulator')
		#	stateAnt   = readState(state_file)
		#	for i in range(len(stateAnt)):
		#		Xn.push_back(stateAnt[i])
			
		cdef doublevec Xn1 = doublevec_factory(0)
		
		onlyAssert(sargs,'rpms','Simulator')
		cdef doublevec rpms = doublevec_factory(0)
		for i in range(len(sargs['rpms'])):
			rpms.push_back(sargs['rpms'][i])
		cdef int nrpms = len(sargs['rpms'])
		
		cdef double dtheta_rpm = validatePositive(sargs,'dtheta_rpm','Simulator')
		cdef int inicia	       = assignOptional(sargs, 'inicia', 1)
		cdef double Courant    = validatePositive(sargs,'Courant','Simulator')
		cdef double ga	       = validatePositive(sargs,'ga','Simulator')
		cdef int viscous_flow  = onlyAssert(sargs,'viscous_flow','Simulator')
		cdef int heat_flow     = onlyAssert(sargs,'heat_flow','Simulator')
		cdef double R_gas      = assignOptional(sargs,'R_gas',286.9)
		cdef int nstroke       = validatePositive(sargs,'nstroke','Simulator',4)
		cdef int ncycles       = validatePositive(sargs,'ncycles','Simulator')
		cdef intvec ig_order   = intvec_factory(0)
		for i in range(ncyl):
			ig_order.push_back(sargs['ig_order'][i])
			
		cdef int engine_type	= validateInList(sargs,'engine_type','Simulator',[0,1],0)
		
		s = assignOptional(sargs, 'filein_state', '')
		cdef char* filein_state = s
		s = assignOptional(sargs, 'filesave_state', '')
		cdef char* filesave_state = s
		s = assignOptional(sargs, 'filesin_spd', '')
		cdef char* filein_spd = s
		s = assignOptional(sargs, 'filesave_spd', '')
		cdef char* filesave_spd = s
		s = assignOptional(sargs, 'folder_name', 'testDefault')
		cdef char* folder_name  = s

		cdef int calc_engine_data	= boolean(sargs,'calc_engine_data','Simulator',1)
    	
		#creo los tubos
		cdef tubevec tubes = tubevec_factory(0)
		cdef Tube auxTube
		for k in range(ntubes):
			auxTube = Tube(**dataTubes[k])
			tubes.push_back(copyTube(auxTube.thisptr))
		
		#creo los cilindros
		cdef double delta_encendido
		if(ncyl == 0):
			delta_encendido = 0.0
		else:
			delta_encendido = nstroke*pi/ncyl
			
		cdef double theta_start    = 0.0
		cdef cylindervec cylinders = cylindervec_factory(0)
		cdef Cylinder auxCyl
		for k in range(ncyl):
			for l in range(ncyl):
				icyl = sargs['ig_order'][l]
				if(icyl==k): break
			theta_0			    = theta_start-l*delta_encendido
			theta_0             = mod(theta_0,nstroke*pi)
			dataCylinders[k]['theta_0'] = theta_0
			auxCyl			    = Cylinder(**dataCylinders[k])
			cylinders.push_back(copyCylinder(auxCyl.thisptr))
		    	
		#creo los tanques
		cdef tankvec tanks = tankvec_factory(0)
		cdef Tank auxTank
		for k in range(ntank):
			auxTank = Tank(**dataTanks[k])
			tanks.push_back(copyTank(auxTank.thisptr))
		
		#creo las uniones
		cdef junctionvec junctions = junctionvec_factory(0)
		cdef Junction auxJunc
		for k in range(njunc):
			auxJunc = Junction(**dataJunctions[k])
			junctions.push_back(copyJunction(auxJunc.thisptr))	
			
		#creo las atmospheras
		cdef atmvec atmospheres = atmvec_factory(0)
		cdef Atmosphere auxAtm
		for k in range(natm):
			auxAtm = Atmosphere(**dataAtmospheres[k])
			atmospheres.push_back(copyAtmosphere(auxAtm.thisptr))
			
		self.thisptr = new_Simulator(dt,tf,nrpms,rpms,Xn,Xn1,ntubes,ncyl,ntank,njunc,iter_sim1d,nsave,nappend,
					     dtheta_rpm,inicia,Courant,ga,viscous_flow,heat_flow,R_gas,nstroke,ncycles,
					     engine_type,filein_state,filesave_state,filein_spd,filesave_spd,folder_name,
					     ig_order,cylinders,tubes,junctions,tanks,natm,atmospheres, get_state, 
					     calc_engine_data)
		
	def solver(self):
		self.thisptr.solver()

	def solve_step(self):
		self.thisptr.solveStep()
	
	def printData(self):
		self.thisptr.printData()

	def ActualizeGlobalState(self):
		self.thisptr.actualizeState()
        
	def actualizeDt(self):
		self.thisptr.actualizeDt()
		
	def GetTimeStep(self):
		return float(self.thisptr.GetTimeStep())
	
	def SetTimeStep(self, dT):
		self.thisptr.SetTimeStep(dT)
		
	def SetTime(self, t, irpm=0):
		self.thisptr.SetTime(t, irpm)
		
	def GetNewState(self):
		cdef doublevec a = self.thisptr.GetNewState()
		cdef int l = self.thisptr.GetNewStateSize()
		aPython = []	
		for i in range(l):
			aPython.append(a.at(i))
		return aPython

	def GetOldState(self):
		cdef doublevec a = self.thisptr.GetOldState()
		cdef int l = self.thisptr.GetOldStateSize()
		aPython = []			
		for i in range(int(l)):
			aPython.append(a.at(i))
		return aPython

	def SetState(self, indx, values):
		for i in range(len(indx)):
			self.thisptr.SetStateValue(indx[i],values[i])

	def GetCylinderMass(self, icyl):
		cdef doublevec mass = self.thisptr.getCylinderMass(icyl)
		lmass = []	
		for i in range(6):
			lmass.append(mass.at(i))
		return lmass

	def SetCylinderMass(self, icyl, mass_vec):
		for i in range(len(mass_vec)):
			self.thisptr.setCylinderMass(icyl, i, mass_vec[i])

	def __dealloc__(self):
		del_Simulator(self.thisptr)
