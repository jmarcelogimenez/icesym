import default_values
import os, sys
import string
from calcs import *
from numpy import trapz,pi,zeros
from unitsData import *

#def upGrid(grid):
#	rows = grid.GetNumberRows()
#	cols = grid.GetNumberCols()
#	ret = []
#	for i in range(cols):
#			aux = []
#			for j in range(rows):
#			if not(grid.GetCellValue(j,i) == ''):
#					dat = float(grid.GetCellValue(j,i))
#					aux.append(dat)
#			if not(aux==[]):
#				if cols>1:
#					ret.append(aux)
#				else:
#					ret = aux
#	return ret

def upGrid(grid):
	rows = grid.GetNumberRows()
	cols = grid.GetNumberCols()
	ret = []
	for i in range(rows):
			aux = []
			for j in range(cols):
				if not(grid.GetCellValue(i,j) == ''):
					dat = float(grid.GetCellValue(i,j))
					aux.append(dat)
			if not(aux==[]):
				if cols>1:
					ret.append(aux)
				else:
					ret = ret + aux
	return ret

def setGrid(data,grid):
	rows = grid.GetNumberRows()
	cols = grid.GetNumberCols()
	if isinstance(data,list) and not(data==[]):
		ldata = len(data)
		if isinstance(data[0],list): #caso matriz 2D
			ldata0 = len(data[0])
			if ldata0-cols > 0:
				grid.AppendCols(ldata0-cols)
			else:
				grid.DeleteCols(0,cols-ldata0)
			if ldata-rows > 0:
				grid.AppendRows(ldata-rows)
			else:
				grid.DeleteRows(0,rows-ldata) 
		else: #caso matriz 1D
			if ldata-rows > 0:
				grid.AppendRows(ldata-rows)
			else:
				grid.DeleteRows(0,rows-ldata)

		for i in range(ldata):
			if isinstance(data[i],list):
				for j in range(len(data[i])):
					grid.SetCellValue(i,j,str(data[i][j]))
			else:
				grid.SetCellValue(i,0,str(data[i]))
	return


def getWriteList(objectName, objectData, dict_generic, numCyl = -1):
	lines = []
	line = "\n#--------- Inicializacion de " + objectName + "\n\n"
	lines.append(line)
	if objectName in ['Cylinders','Junctions','Tubes','Tanks', 'Valves', 'Atmospheres']:
		line = objectName + " = []\n\n"
		lines.append(line)
	#if objectName in ['combustion','injection','fuel']:
	#print objectData
	for i in range(len(objectData)):
		line = objectName + str(i) + " = dict()\n"
		lines.append(line)
  	
		if objectName == "Atmospheres":
			line = objectName + str(i) + "['nnod'] = 1\n"
			lines.append(line)
			line = objectName + str(i) + "['ndof'] = 3\n"
			lines.append(line)
			line = objectName + str(i) + "['state_ini'] = [" +  str(objectData[i]['state_ini'][0]) + ", " +  str(objectData[i]['state_ini'][1]) + ", " +  str(objectData[i]['state_ini'][2]) + "]\n"
			lines.append(line)
			line = objectName + str(i) + "['position'] = (" + str(objectData[i]['position'][0]) + ","+ str(objectData[i]['position'][1]) +")\n" #guardo la posicion en el diagrama
			lines.append(line)
		else:
			for key in dict_generic:
				#print defaultNames[j], ": ", objectData[i][defaultNames[j]]
            	# imprime cada uno de los atributos
				if not(objectData[i][key]=="<none>"):
					line = objectName + str(i) + "['" + key + "'] = " + explodeAtribute(objectData[i][key]) + "\n"
					lines.append(line)

			if (objectName == "Cylinders"):
				f = [objectData[i]['fuel']]
				c = [objectData[i]['combustion']]
				inj = [objectData[i]['injection']]
				lines.append("\n")
				lines = lines + getWriteList("fuel", f, default_values.dict_fuel, i)
				lines = lines + getWriteList("combustion", c, default_values.dict_combustion, i)
				lines = lines + getWriteList("injection", inj, default_values.dict_injection, i)
			
            # imprime posicion en la grilla (no sirve en la simulacion pero si al levantar)
			if objectName not in ['Simulator','fuel','combustion','injection']:
	 			line = objectName + str(i) + "['position'] = (" + str(objectData[i]['position'][0]) + ","+ str(objectData[i]['position'][1]) +")\n" #guardo la posicion en el diagrama
				lines.append(line)

            # ACA IMPRIMIR CONEXIONES!!!
			if objectName == "Tubes":
				line = objectName + str(i) + "['nleft'] = " + str(objectData[i]['nleft']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['tleft'] = '" + objectData[i]['tleft'] + "'\n" #guardo la posicion en el diagrama
				lines.append(line)
				line = objectName + str(i) + "['nright'] = " + str(objectData[i]['nright']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['tright'] = '" + objectData[i]['tright'] + "'\n" #guardo la posicion en el diagrama
				lines.append(line)

			if objectName == "Junctions":
				line = objectName + str(i) + "['type_end'] = " + explodeAtribute(objectData[i]['type_end']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['node2tube'] = " + explodeAtribute(objectData[i]['node2tube']) +"\n" #
				lines.append(line)

			if objectName == "Tanks":
				line = objectName + str(i) + "['int2tube'] = " + explodeAtribute(objectData[i]['int2tube']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['exh2tube'] = " + explodeAtribute(objectData[i]['exh2tube']) +"\n" #
				lines.append(line)
 
			if objectName == "Cylinders":
				line = objectName + str(i) + "['exhaust_valves'] = []\n"
				lines.append(line)
				line = objectName + str(i) + "['intake_valves'] = []\n"
				lines.append(line)

			if objectName == "Valves":
				line = objectName + str(i) + "['tube'] = " + explodeAtribute(objectData[i]['tube']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['ncyl'] = " + explodeAtribute(objectData[i]['ncyl']) +"\n" #
				lines.append(line)
				line = objectName + str(i) + "['typeVal'] = " + explodeAtribute(objectData[i]['typeVal']) +"\n" #
				lines.append(line)
				if(objectData[i]['typeVal'] == "int" and not(objectData[i]['ncyl'] == -1)):
					line = "Cylinders" + str(objectData[i]['ncyl']) + "['intake_valves'].append(" + objectName + str(i) +")\n"
					lines.append(line)
				if(objectData[i]['typeVal'] == "exh" and not(objectData[i]['ncyl'] == -1)):
					line = "Cylinders" + str(objectData[i]['ncyl']) + "['exhaust_valves'].append(" + objectName + str(i) +")\n"
					lines.append(line)
				line = "Valves.append(" + objectName + str(i) +")\n"
				lines.append(line)
				line = "\n"
				lines.append(line)

			if objectName in ["Junctions"]:
				line = objectName + str(i) + "['state_ini'] = " + explodeAtribute(objectData[i]['state_ini'])
				lines.append(line)

		#print "object name: ", objectName
		if objectName in ['Cylinders','Junctions','Tubes','Tanks', 'Atmospheres']:
			line = "\n" + objectName + ".append(" + objectName +  str(i) + ")\n\n"
			lines.append(line)
		elif objectName in ["fuel", "injection","combustion"]:
				line = "\nCylinders" + str(numCyl) + "['" + objectName + "'] = " + objectName +  str(i) + "\n\n"
				lines.append(line)
		elif objectName in ['Simulator']:
				line = "\nSimulator = "	 + objectName +  str(i) + "\n\n"
				lines.append(line)

	line = "\n#--------- FIN Inicializacion de " + objectName + "\n\n"
	lines.append(line)

	return lines

def explodeAtribute(atribute):
	line = ''
	if not(isinstance(atribute,list)):
		if isinstance(atribute,str) or isinstance(atribute,unicode):
			line = line + "'" + atribute + "'"		
		else:
			line = line + str(atribute)
		return line
	else:	
		line = line + "["
		for l in range(len(atribute)):
			line = line + explodeAtribute(atribute[l])			
			if l+1 < len(atribute):
				line = line + ", "
		line = line + "]"
	return line
			

#def loadData(namefile,nRows,nCols):
#	archi = open(namefile, "r")
# data = []
#	linea = archi.readline()
#	lista = string.split(linea)
#	for i in range(len(lista)):
#		data.append([])
#		data[i].append(float(lista[i]))
#	for linea in archi.readlines():	
#		lista = string.split(linea)
#		# uso esta forma para armar la matriz ya transpuesta!
#		for i in range(len(lista)):
#			data[i].append(float(lista[i]))
#		#data.append(lista)
#	archi.close()
#	if not(nRows=="n"):
#		if not(nRows==len(data[0])):
#			return -1
#	if not(nCols=="n"):
#		if not(nCols==len(data)):
#			return -1
#	if len(data)==1:
#		return data[0]
#	else:
#		return data

def loadData(namefile,nRows,nCols):
	archi = open(namefile, "r")
	data = []
	for linea in archi.readlines():	
		lista = string.split(linea)
		data.append(lista)
	archi.close()
	if not(nCols=="n"):
		if not(nCols==len(data[0])):
			return -1
	if not(nRows=="n"):
		if not(nRows==len(data)):
			return -1
	if len(data)==1:
		return data[0]
	else:
		return data

def data2tuple(grid):
	rows = grid.GetNumberRows()
	cols = grid.GetNumberCols()
	data = []
	for i in range(rows):
		try:
			t = ( float(grid.GetCellValue(i,0)),float(grid.GetCellValue(i,1)) )
		except ValueError:
			wx.MessageBox("Incorrect Data", "Error")
			return false
		data.append(t)
	return data

def String2List(string):
	data = []
	aux = ""
	if not(string=="0") and not(string==""):
		for i in range(len(string)):
			if string[i] == ",":
				data.append(int(aux))
				aux = ""
			else:
				aux = aux + string[i]
		data.append(int(aux))
	return data

def List2String(lista):
	string = ""
	for l in range(len(lista)):
		if l < len(lista)-1:
			string = string + str(lista[l]) + ","
		else:
			string = string + str(lista[l])
	return string

def List2StringList(lista):
	string = []	
	for l in range(len(lista)):
		string.append(str(lista[l]))	
	return string

# este algoritmo permite generar el state de los cylinders (usado antes de guardar)
# agrega al state los estados de las valvulas
def couplingCode(home):
	genStates(home)
	for i in range(len(home.Cylinders)):
		n = 0
		histo = []
		if home.Cylinders[i]['ownState'] == 1:
			histo.append(0)
		home.Cylinders[i]['nvi'] = 0
		home.Cylinders[i]['nve'] = 0
		if home.Cylinders[i]['state_ini'] == '':
			home.Cylinders[i]['state_ini'] = [[0,0,0]]
		for v in home.Valves:
			if v['ncyl'] == i and v['typeVal']=='int':
				n = n+1
				#for j in range(home.Cylinders[i]['ndof']):
				#	home.Cylinders[i]['state_ini'][j] = home.Cylinders[i]['state_ini'][j] + v['state_ini'][j]
				home.Cylinders[i]['state_ini'] = home.Cylinders[i]['state_ini'] + v['state_ini']
				home.Cylinders[i]['nvi'] = home.Cylinders[i]['nvi'] + 1
				#print "para: ",n," valor: ",bool(v['histo']==1)
				if bool(v['histo']==1) == True:
					histo.append(n)
		for v in home.Valves:
			if v['ncyl'] == i and v['typeVal']=='exh':
				n = n+1
				#for j in range(home.Cylinders[i]['ndof']):
				#	home.Cylinders[i]['state_ini'][j] = home.Cylinders[i]['state_ini'][j] + v['state_ini'][j]
				home.Cylinders[i]['state_ini'] = home.Cylinders[i]['state_ini'] + v['state_ini']
				home.Cylinders[i]['nve'] = home.Cylinders[i]['nve'] + 1
				if bool(v['histo']==1) == True:
					histo.append(n)
		home.Cylinders[i]['nnod'] = home.Cylinders[i]['nnod'] + n
		home.Cylinders[i]['histo'] = histo

# este algoritmo permite generar el state de los cylinders para separarlos de las valves
# usado para ver en GUI
def decouplingCode(home):
	for i in range(len(home.Cylinders)):
		ns = home.Cylinders[i]['nnod'] - home.Cylinders[i]['nvi'] - home.Cylinders[i]['nve']
		new = []		
		for j in range(ns):
			new.append([])
			for k in range(home.Cylinders[i]['ndof']):
				new[j].append(home.Cylinders[i]['state_ini'][j][k])
		home.Cylinders[i]['state_ini'] = new
		home.Cylinders[i]['nnod'] = ns

# carga los datos respecto a theta
# extras = -1 .... dentro de state
# extras = 1 ... extras de cylinder
# extras = 2 ... extras de tank
def loadAngleTXT(namefile,cycle,nnod,ndof,extras=-1):
	namecopy = decompress(namefile)
	archi = open(namecopy, "r")
	data = []
	if extras == -1:
		for linea in archi.readlines():	
			lista = string.split(linea)
			if int(cycle) == int(lista[1]) and int(nnod) == int(lista[0]):
				t = (float(lista[2]),float(lista[4+ndof]))
				data.append(t)
		archi.close()
	elif extras == 1:
		ndof = ndof - 3		#los 3 primeros son del estado
		#print "ndof buscado: ",ndof
		row = 3
		col = ndof - 6
		if ndof in range(4):
			row = 1 #segunda fila de cada paso de tiempo
			col = ndof
		elif ndof in [4,5]:
			row = 2 #tercera fila de cada paso de tiempo
			col = ndof - 4
		#print "row: ", row
		#print "col: ", col
		#print "cycle: ", cycle
		line = 0
		read = False
		theta = 0
		for linea in archi.readlines():	
			if line%4 == 0: #leo cada 4 lineas el ciclo (esta en la primer linea)
				lista = string.split(linea)
				#print "linea leida: ", lista
				if int(cycle) == int(lista[0]):
					#print "linea leida: ", lista
					read = True
					theta = float(lista[1])
				else:
					read = False
			if line%4 == row: #leo cada 4 lineas la variable
				if read: 
					lista = string.split(linea)
					t = (theta,float(lista[col]))
					#print t
					data.append(t)
			line = line + 1
		archi.close()
	elif extras == 2:
		ndof = ndof - 3
		for linea in archi.readlines():	
			lista = string.split(linea)
			if int(cycle) == int(lista[0]):
				t = (float(lista[1]),float(lista[3+ndof]))
				data.append(t)
		archi.close()
	os.system("rm -r tests")
	return data


def loadCycleTXT(namefile,nnod,ndof,extras=-1):
	namecopy = decompress(namefile)
	archi = open(namecopy, "r")
	data = []
	cycles = []
	if extras == -1:
		for linea in archi.readlines():	
			lista = string.split(linea)
			c = int(lista[1])
			if len(cycles) < c:
				cycles.append([])
			if int(nnod) == int(lista[0]):
				cycles[c-1].append((float(lista[2]),float(lista[4+ndof])))
	elif extras == 1:
		ndof = ndof - 3		#los 3 primeros son del estado
		#print "ndof buscado: ",ndof
		row = 3 # por defecto cuarta fila
		col = ndof - 6
		if ndof in range(4):
			row = 1 #segunda fila de cada paso de tiempo
			col = ndof
		elif ndof in [4,5]:
			row = 2 #tercera fila de cada paso de tiempo
			col = ndof - 4
		line = 0
		theta = 0
		for linea in archi.readlines():	
			#print "linea: ", line
			if line%4 == 0: #leo cada 4 lineas el ciclo (esta en la primer linea)
				lista = string.split(linea)
				c = int(lista[0])
				theta = float(lista[1])
				if len(cycles) < c:
					cycles.append([])
			if line%4 == row: #leo cada 4 lineas la variable
				lista = string.split(linea)
				t = (theta,float(lista[col]))
				#print t
				cycles[c-1].append(t)
			line = line + 1
	elif extras == 2:
		ndof = ndof - 3
		for linea in archi.readlines():	
			lista = string.split(linea)
			c = int(lista[0])
			if len(cycles) < c:
				cycles.append([])
			t = (float(lista[1]),float(lista[3+ndof]))
			cycles[c-1].append(t)

	archi.close()
	#data = meanCycles(cycles)
	data = trapCycles(cycles)
	os.system("rm -r tests")
	return data

def loadSpaceTXT(namefile,time,ndof):
	data = []
	nearest = -1000
	tn = 0

	namecopy = decompress(namefile)
	archi = open(namecopy, "r")
	for linea in archi.readlines():	
		lista = string.split(linea)
		if abs(float(time)-float(lista[3])) < nearest:
			nearest = abs(float(time)-float(lista[3]))
			tn = float(lista[3])
 	archi.close()

	archi = open(namefile, "r")
	for linea in archi.readlines():	
		lista = string.split(linea)
		if float(lista[3]) == tn:
			t = (float(lista[7]),float(lista[4+ndof]))
			data.append(t)
	archi.close()
	os.system("rm -r tests")
	return data

def calcGlobalsData(path, dataRPMS, pd):
	rpms = dataRPMS['rpm']
	cyls = dataRPMS['Cylinders']
	nstroke = dataRPMS['nstroke']
	rho = dataRPMS['rho']
	nr = nstroke/2
	Q_fuel = cyls['Q_fuel']
	ncyls = len(cyls['labels'])
	#obtengo el trabajo en cada cilindro promediando los last cycles de cada rpm		
	#obtengo el volumen desplazado por ciclo (volMax - volMin)
	ncycles = len(dataRPMS['cycle'])
	w = [] #(ncyls x nrpms x ncycles)
	volDesp = [] #(ncyls x nrpms x ncycles)
	pMax = []  #(ncyls x nrpms x ncycles)
	valPD = 0
	maxPD = ncyls * len(rpms) + 1

	for j in range(ncyls):
		w.append([])	
		volDesp.append([])	
		pMax.append([])
		for irpm in rpms:
			w[-1].append([])	
			volDesp[-1].append([])	
			pMax[-1].append([])
			for icycle in range(ncycles):
				#volumen
				dataFile = path + "RPM_" + str(irpm) + "/cyl_extras_" + str(j) + ".txt"	
				ndofVol = 11
				volData = loadAngleTXT(dataFile,icycle+1,0,ndofVol,1)
				#presion
				dataFile = path + "RPM_" + str(irpm) + "/cyl_" + str(j) + ".txt"	
				ndofp = 1
				pData = loadAngleTXT(dataFile,icycle+1,0,ndofp,-1)
				
				work = []
				vd = []
				pMaxAux = []	
			
				#print volData	
			
				#primeros valores
				firstV = (volData[0][1] + volData[-1][1])/2
				firstp = (pData[0][1] + pData[-1][1])/2
				dV = [firstV]
				p = [firstp]
			
				#valores intermedios
				for k in range(len(volData)):
					dV.append(volData[k][1])
					p.append(pData[k][1])
				
				#ultimos valores
				dV.append(firstV)
				p.append(firstp)

				vMax = max(dV)
				vMin = min(dV)
			
				vd = (vMax-vMin)#volumen desplazado en el ciclo	    
				work = trapz(p,dV) #trabajo efecutuado en el ciclo
				pMaxAux = max(p) #presion maxima en el ciclo
				
				w[-1][-1].append(work)
				volDesp[-1][-1].append(vd)
				pMax[-1][-1].append(pMaxAux)

			valPD = valPD + 1
			#pd.Update(valPD,"Loading RPM " + str(irpm) + " for cylinder " + str(j))
			pd.Update(valPD,"Loading ...")


	#ordenes de magnitud bastante bien	
	#print "w: ", w
	#print "volDesp: ", volDesp
	#print "pMax: ", pMax
	#sys.exit()

	#calculo de IMEP_per_cylinder, FMEP_per_cylinder y BMEP_per_cylinder (ncyls x ncycles x nrpms)
	IMEP_per_cylinder = []
	c1 = 0.4 #bar
	c2 = 0.005 
	c3 = 0.09 #bar/(m/s)
	c4 = 9e-4 #bar/(m/s)^2
	FMEP_per_cylinder = []
	BMEP_per_cylinder = []
	Bar2Pa = 1e6
	for j in range(ncyls):
		IMEP_per_cylinder.append([])
		FMEP_per_cylinder.append([])
		BMEP_per_cylinder.append([])
		for icycle in range(ncycles):
			IMEP_per_cylinder[-1].append([])
			FMEP_per_cylinder[-1].append([])
			BMEP_per_cylinder[-1].append([])
		for irpm in range(len(rpms)):
			sm = float(rpms[irpm])*2*cyls['crank_radius'][j]/30
			for icycle in range(ncycles):
				aux1 = w[j][irpm][icycle]/volDesp[j][irpm][icycle]
				aux2 = (c1 + c2*pMax[j][irpm][icycle]/Bar2Pa + c3*sm + c4*sm*sm) * Bar2Pa
				aux3 = aux1 - aux2
				IMEP_per_cylinder[-1][icycle].append((float(rpms[irpm]),aux1))
				FMEP_per_cylinder[-1][icycle].append((float(rpms[irpm]),aux2))		
				BMEP_per_cylinder[-1][icycle].append((float(rpms[irpm]),aux3))		

	#print "IMEP_per_cylinder: ", IMEP_per_cylinder
	#print "FMEP_per_cylinder: ", FMEP_per_cylinder
	#print "BMEP_per_cylinder: ", BMEP_per_cylinder
	
	#sys.exit()

	#calculo de Potencia y Torque (efectivos e indicados) (ncycles x nrpms)
	power_indicated = []
	torque_indicated = []
	power_effective = []
	torque_effective = []
	mechanical_efficiency = []

	for icycle in range(ncycles):
			power_indicated.append([])
			torque_indicated.append([])
			power_effective.append([])
			torque_effective.append([])
			mechanical_efficiency.append([])

	
	for icycle in range(ncycles):
		for irpm in range(len(rpms)):
			N = float(rpms[irpm])/60
			aux1 = 0
			aux3 = 0
			for j in range(ncyls):
				Vd = volDesp[j][irpm][icycle]
				aux1 = aux1 + IMEP_per_cylinder[j][icycle][irpm][1]*Vd*N/nr
				aux3 = aux3 + BMEP_per_cylinder[j][icycle][irpm][1]*Vd*N/nr
			aux4 = aux3/(2*pi*N)	
			aux2 = aux1/(2*pi*N)
			aux5 = aux3/aux1			
			power_indicated[icycle].append((float(rpms[irpm]),aux1))
			torque_indicated[icycle].append((float(rpms[irpm]),aux2))
			power_effective[icycle].append((float(rpms[irpm]),aux3))
			torque_effective[icycle].append((float(rpms[irpm]),aux4))
			mechanical_efficiency[icycle].append((float(rpms[irpm]),aux5))
		#print "por ciclo: ", power_indicated[icycle]

	#print "final: ", power_indicated


	#calculo de consumo especifico y rendimientos mecanicos, volumetricos y de conversion de combustible
	(mfc,mair)  = getMasses(path, dataRPMS)
	#(mfc,mair)  = (1,1)
		
	SFC_indicated = []
	SFC_effective = []
	fuel_conversion_efficiency_indicated = []
	fuel_conversion_efficiency_effective = []

	for icycle in range(ncycles):
		SFC_indicated.append([])
		SFC_effective.append([])
		fuel_conversion_efficiency_indicated.append([])
		fuel_conversion_efficiency_effective.append([])

	for icycle in range(ncycles):
		for irpm in range(len(rpms)):
			mfcTotal = 0
			for j in range(ncyls):
				mfcTotal = mfcTotal + mfc[j][irpm][icycle]
			aux1 = mfcTotal*N/(nr*power_indicated[icycle][irpm][1])
			aux2 = mfcTotal*N/(nr*power_effective[icycle][irpm][1])
			aux4 = 1/(aux1*Q_fuel)
			aux5 = 1/(aux2*Q_fuel)
			SFC_indicated[icycle].append((float(rpms[irpm]),aux1))
			SFC_effective[icycle].append((float(rpms[irpm]),aux2))
			fuel_conversion_efficiency_indicated[icycle].append((float(rpms[irpm]),aux4))
			fuel_conversion_efficiency_effective[icycle].append((float(rpms[irpm]),aux5))

	
	volumetric_efficiency = []
	for j in range(ncyls):
		volumetric_efficiency.append([])
		for icycle in range(ncycles):
			volumetric_efficiency[-1].append([])
		for irpm in range(len(rpms)):
			for icycle in range(ncycles):
				Vd = volDesp[j][irpm][icycle]
				aux = mair[j][irpm][icycle]/(rho*Vd)
				volumetric_efficiency[-1][icycle].append((float(rpms[irpm]),aux))

	volumetric_efficiency_global = []
	IMEP_global = []
	BMEP_global = []
	FMEP_global = []

	for icycle in range(ncycles):
		volumetric_efficiency_global.append([])
		SFC_effective.append([])
		IMEP_global.append([])
		BMEP_global.append([])
		FMEP_global.append([])

	for icycle in range(ncycles):
		for irpm in range(len(rpms)):
			aux1 = 0
			aux2 = 0
			aux3 = 0
			aux4 = 0
			for j in range(ncyls):
				aux1 = aux1 + volumetric_efficiency[j][icycle][irpm][1]
				aux2 = aux2 + IMEP_per_cylinder[j][icycle][irpm][1]
				aux3 = aux3 + BMEP_per_cylinder[j][icycle][irpm][1]
				aux4 = aux4 + FMEP_per_cylinder[j][icycle][irpm][1]
			aux1 = aux1/ncyls
			aux2 = aux2/ncyls
			aux3 = aux3/ncyls
			aux4 = aux4/ncyls
			volumetric_efficiency_global[icycle].append((float(rpms[irpm]),aux1))
			IMEP_global[icycle].append((float(rpms[irpm]),aux2))
			BMEP_global[icycle].append((float(rpms[irpm]),aux3))
			FMEP_global[icycle].append((float(rpms[irpm]),aux4))

	valPD = maxPD
	pd.Update(valPD,"Loading RPM " + str(irpm) + " for cylinder " + str(j))
	data = dict()
	data['IMEP_per_cylinder'] = IMEP_per_cylinder
	data['FMEP_per_cylinder'] = FMEP_per_cylinder
	data['BMEP_per_cylinder'] = BMEP_per_cylinder
	data['IMEP_global'] = IMEP_global
	data['FMEP_global'] = FMEP_global
	data['BMEP_global'] = BMEP_global
	data['power_indicated'] = power_indicated
	data['torque_indicated'] = torque_indicated
	data['power_effective'] = power_effective
	data['torque_effective'] = torque_effective
	data['mechanical_efficiency'] = mechanical_efficiency
	data['SFC_indicated'] = SFC_indicated
	data['SFC_effective'] = SFC_effective
	data['volumetric_efficiency'] = volumetric_efficiency
	data['volumetric_efficiency_global'] = volumetric_efficiency_global
	data['fuel_conversion_efficiency_indicated'] = fuel_conversion_efficiency_indicated
	data['fuel_conversion_efficiency_effective'] = fuel_conversion_efficiency_effective
	
	return data
		

def getMasses(path, dataRPMS):
	rpms = dataRPMS['rpm']
	cyls = dataRPMS['Cylinders']
	ncyls = len(dataRPMS['Cylinders']['labels'])
	nstroke = dataRPMS['nstroke']
	ncycles = len(dataRPMS['cycle'])
	mfc = []	
	mair = []
	for j in range(ncyls):	
		mfc.append([])		
		mair.append([])
		angleClose = cyls['angleClose'][j]
		for irpm in rpms:
			mfc[-1].append([])		
			mair[-1].append([])
			for icycle in range(ncycles):
				namefile = path + "RPM_" + str(irpm) + "/cyl_extras_" + str(j) + ".txt"	
				ndofmfc = 3 #cuarto de la fila 4
				ndofmair = 4 #quinto de la fila 4
				row = 3
				namecopy = decompress(namefile)
				archi = open(namecopy, "r")
				data = []
				mfcAux = 0
				mairAux = 0
				line = 0
				next = True
				interpola = True
				for linea in archi.readlines():	
					if line%4 == 0: #leo cada 4 lineas el ciclo (esta en la primer linea)
						lista = string.split(linea)
						if int(icycle) == int(lista[0]) and angleClose < float(lista[1]):
							next = True
						elif int(icycle) == int(lista[0]) and angleClose >= float(lista[1]):
							next = False
							if angleClose == float(lista[1]):
								interpola = False
						
					if line%4 == row: #leo cada 4 lineas la variable
						lista = string.split(linea)
						if next == True:
							mfcAux = float(lista[ndofmfc])
							mairAux = float(lista[ndofmair])
						else: 
							if interpola == True:
								mfcAux = (mfcAux+float(lista[ndofmfc]))/2
								mairAux = (mairAux+float(lista[ndofmair]))/2
							else:
								mfcAux = float(lista[ndofmfc])
								mairAux = float(lista[ndofmair])
							break	
					line = line + 1			
		
				mfc[-1][-1].append(mfcAux)		
				mair[-1][-1].append(mairAux)
				archi.close()
				
				os.system("rm -r tests")

	return (mfc,mair)
	

# si trabajo dentro de la misma carpeta (sim-c/src/) no habria problemas en esto
def decompress(namefile):
	i1 = namefile.rfind("tests")
	namecopy = ''	
	for j in range(len(namefile) - i1):
		namecopy = namecopy + namefile[j+i1]
	os.system("tar jvxf " + namefile + ".tar.bz2")
	return namecopy


#genera los estados iniciales de los elementos Tank, Valve, Junction, buscando en los tubos conectados
#si no estan conectados, asigna 0 a cada ndof
def genStates(home):
	for i in range(len(home.Junctions)):
		home.Junctions[i]['state_ini'] = []
		for j in range(len(home.Junctions[i]['node2tube'])):
			if home.Junctions[i]['type_end'] == 1:
				home.Junctions[i]['state_ini'].append(home.Tubes[home.Junctions[i]['node2tube'][j]]['state_ini'][-1])
			else:
				home.Junctions[i]['state_ini'].append(home.Tubes[home.Junctions[i]['node2tube'][j]]['state_ini'][0])
	for i in range(len(home.Tanks)):
		home.Tanks[i]['state_ini'] = [home.Tanks[i]['state_ini'][0]]
		for j in range(len(home.Tanks[i]['int2tube'])):
			home.Tanks[i]['state_ini'].append(home.Tubes[home.Tanks[i]['int2tube'][j]]['state_ini'][-1])
		for j in range(len(home.Tanks[i]['exh2tube'])):
			home.Tanks[i]['state_ini'].append(home.Tubes[home.Tanks[i]['exh2tube'][j]]['state_ini'][0])
			
	for i in range(len(home.Valves)):
		if not(home.Valves[i]['tube'] == -1):
			if home.Valves[i]['typeVal'] == 'exh':
				home.Valves[i]['state_ini'] = [home.Tubes[home.Valves[i]['tube']]['state_ini'][0]]
			else:
				home.Valves[i]['state_ini'] = [home.Tubes[home.Valves[i]['tube']]['state_ini'][-1]]
		else:
			home.Valves[i]['state_ini'] = [zeros(3)]

def getWriteListPost(objectName, objectData):
	lines = []
	#print objectData
	line = objectName + " = []\n"
	lines.append(line)
	for i in range(len(objectData)):
		line = objectName + str(i) + " = []\n"
		lines.append(line)
		for j in range(len(objectData[i])):
			line = objectName + str(i) + str(j) + " = []\n"
			lines.append(line)
			for key in objectData[i][j]:
				line = objectName + str(i) + str(j) + "['" + key + "'] = " + explodeAtribute(objectData[i][j][key]) + "\n"
				lines.append(line)
			line = objectName + str(i) + ".append(" + objectName + str(i) + str(j) + ")\n\n"
			lines.append(line)
		line = objectName + ".append(" + objectName + str(i) + ")\n\n"
		lines.append(line)
		
	line = "\n#--------- FIN Inicializacion de " + objectName + "\n\n"
	lines.append(line)

	return lines


