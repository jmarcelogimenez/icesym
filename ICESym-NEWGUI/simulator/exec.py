from numpy import array
import validator
import fuel_lib
import myfunction
import time
import sys
from simCythonCPP import Simulator

sys.path.append(sys.argv[1])
data = __import__(sys.argv[2])
now = time.time()
Sim = Simulator(**data.kargs)
print "termina de inicializar"
#Sim.printData()
Sim.solver()
now2 = time.time()
print now2-now
