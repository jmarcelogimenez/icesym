from numpy import array
import time
import sys
import os
from simCythonCPP import Simulator

sys.path.append("/home/etekken/Cimec/icesym/NewGUI/test2/icesym/ICESym-NEWGUI/saves")
data = __import__("multicyl4TSI")
now = time.time()
Sim = Simulator(**data.kargs)
print "termina de inicializar"
Sim.printData()
Sim.solver()
now2 = time.time()
print now2-now
