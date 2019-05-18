def readState(state_file):
	f = open(state_file, 'r')
	data = []	
	for line in f:
		data.append(float(line))
	return data
