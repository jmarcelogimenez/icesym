from numpy import trapz

#cycles en una matriz de pares
def meanCycles(cycles):
	data = []	
	for i in range(len(cycles)):
		s = 0
		for v in cycles[i]: 
			s = s + v[1]
		s = s/len(cycles[i])
		data.append((i+1,s))
	return data		

def trapCycles(cycles):
	data = []	
	for i in range(len(cycles)):
		x = [0]
		value0 = (cycles[i][0][1]+cycles[i][-1][1]) / 2 #considerando convergido el ciclo
		y = [value0]
		for v in cycles[i]: 
			x.append(v[0])
			y.append(v[1])
		maxAngle = 720		
		if cycles[i][-1][0]<360: #nstroke = 2
			maxAngle = 360
		x.append(maxAngle)	
		y.append(value0)
		res = trapz(y,x)	
		res = res/maxAngle
		data.append((i+1,res))
	return data		

#rpms en una matriz de pares
# lastCycles indica el numero de cyclos a promediar (de atras para adelante)
def meanRPMs(cdata,rpms,lastCycles = 1):
	data = []	
	for i in range(len(cdata)):
		s = 0
		nCycle = 1
		for v in cdata[i]: 
			if (nCycle+lastCycles) > len(cdata[i]): #es alguno de los last cycles
				s = s + v[1]
			nCycle = nCycle + 1
		s = s/lastCycles
		data.append((rpms[i],s))
	return data		


#con los datos separados en rpms que vienen en cdata, se saca un promedio por rpm de cada icycle
#obteniendo un para (rmp, valor_prom_icycle)
def getRPMfromCycle(cdata,rpms,icycle):
	data = []	
	for i in range(len(cdata)):
		s = 0
		nCycle = 1
		for v in cdata[i]: 
			if int(nCycle) == int(icycle):
				data.append((int(rpms[i]),v[1]))
				break
			nCycle = nCycle + 1
	return data	
