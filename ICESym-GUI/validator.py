import sys

def validatePositive(kargs,arg,nameClass,opt=-11):
	if (arg in kargs.keys()):
		if(kargs[arg]>0):
			return kargs[arg]
	else:
		if(opt!=-11):
			return opt
	print 'Fail inicialitation in [%s,%s], number must be positive, or argument not exists' % (arg,nameClass)	
	sys.exit()	
	
def onlyAssert(kargs,arg,nameClass):
	if not(arg in kargs.keys()):
		print 'Fail inicialitation in [%s,%s], argument not exists' % (arg,nameClass)	
		sys.exit()
	else:
		return kargs[arg]
		
def assignOptional(kargs,arg,opt):
	if not(arg in kargs.keys()):
		return opt
	else:
		return kargs[arg]
		
def validateSize(kargs, arg,nameClass,size):
	if(arg in kargs.keys()):
		if(len(kargs[arg])==size):
			return 1
	print 'Fail inicialitation in [%s,%s], argument have not the requered size' % (arg,nameClass)
	sys.exit()
	
def boolean(kargs,arg,nameClass,opt=-11):
	if(arg in kargs.keys()):
		if(kargs[arg]==0 or kargs[arg]==1):
			return kargs[arg]
	else:
		if(opt!=-11):
			return opt
	print 'Fail inicialitation in [%s,%s], argument must be 0 or 1' % (arg,nameClass)
	sys.exit()			
	
def validateInList(kargs,arg,nameClass,lista,opt):
	if(arg in kargs.keys()):
		if(kargs[arg] in lista):
			return kargs[arg]
		else:
			print 'Fail inicialitation in [%s,%s], argument must be in [%s]' % (arg,nameClass,lista)
	else:
		return opt
	sys.exit()	
