import operator

UNITS = {}
UNITS['Density']        = ['kg/m^3']
UNITS['Pressure']       = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['Velocity']       = ['m/s']
UNITS['Temperature']    = ['K',"C"]

UNITS['Convective Heat-Transfer Coeff'] = ['W/(K.m^2)']
UNITS['Radiactive Heat-Transfer Coeff'] = ['W/(K.m^2)']
UNITS['Convective Heat-Transfer Rate'] = ['W','kW']
UNITS['Radiactive Heat-Transfer Rate'] = ['W','kW']
UNITS['Burned Mass Fraction'] = []
UNITS['Burned Mass Fraction Rate'] = []
UNITS['Mass Flow Rate trought Intake Port'] = ['kg/s','g/s']
UNITS['Mass Flow Rate trought Exhaust Port'] = ['kg/s','g/s']
UNITS['Volume'] = ['m^3','cm^3']
UNITS['Mass of Fuel'] = ['kg','g']
UNITS['Mass of Air'] = ['kg','g']
UNITS['Mass of Residual Gas'] = ['kg','g']
UNITS['Total Heat-Transfer Rate'] = ['W','kW']
UNITS['Fuel Chemical Energy Release'] = ['W','kW']
UNITS['Power Indicated'] = ['W', 'kW', 'hp', 'cv']
UNITS['Power Effective'] = ['W', 'kW', 'hp', 'cv']
UNITS['Torque'] = ['N.m','kgf']
UNITS['Torque Indicated'] = ['N.m','kgf']
UNITS['Torque Effective'] = ['N.m','kgf']
UNITS['IMEP per Cylinder'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['IMEP Global'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['FMEP per Cylinder'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['FMEP Global'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['BMEP per Cylinder'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['BMEP Global'] = ['Pa', 'kPa', 'MPa','Bar','Atm']
UNITS['SFC Indicated'] = ['kg/(W.s)', 'g/(kW.h)','lbm/(hp.h)']
UNITS['SFC Effective'] = ['kg/(W.s)', 'g/(kW.h)','lbm/(hp.h)']
UNITS['Mechanical Efficiency'] = []
UNITS['Volumetric Efficiency per Cylinder'] = []
UNITS['Volumetric Efficiency Global'] = []
UNITS['Fuel Conversion Efficiency Indicated'] = []
UNITS['Fuel Conversion Efficiency Effective'] = []

UNITS['Mass Flow Rate'] = ['kg/s','g/s']
UNITS['Enthalpy Flow Rate'] = []
UNITS['Mass'] = ['kg','g']

# Las conversiones son desde el sistema internacional a la unidad requerida
CONVERSIONS = {}
# Presion
CONVERSIONS['Pa'] = [1, operator.mul]
CONVERSIONS['Bar'] = [1e-5, operator.mul]
CONVERSIONS['kPa'] = [1e-3, operator.mul]
CONVERSIONS['MPa'] = [1e-6, operator.mul]
CONVERSIONS['Atm'] = [9.87e-6, operator.mul]
# Energia
CONVERSIONS['W'] = [1, operator.mul]
CONVERSIONS['kW'] = [1e-3, operator.mul]
CONVERSIONS['hp'] = [1.34e-3, operator.mul]
CONVERSIONS['cv'] = [1.36e-3, operator.mul]
# Eficiencia
CONVERSIONS['kg/(W.s)'] = [1, operator.mul]
CONVERSIONS['g/(kW.h)'] = [1e3*3600.0/1e-3, operator.mul]
CONVERSIONS['lbm/(hp.h)'] = [2.2*3600.0/1.34e-3, operator.mul]
# Temperatura
CONVERSIONS['K'] = [1, operator.mul]
CONVERSIONS['C'] = [273.15, operator.sub]
# Volumen
CONVERSIONS['m^3'] = [1, operator.mul]
CONVERSIONS['cm^3'] = [1e6, operator.mul]
# Otras
CONVERSIONS['kg/m^3'] = [1, operator.mul]
CONVERSIONS['m/s'] = [1, operator.mul]
CONVERSIONS['W/(K.m^2)'] = [1, operator.mul]
# Caudal Masico
CONVERSIONS['kg/s'] = [1, operator.mul]
CONVERSIONS['g/s'] = [1e3, operator.mul]
# Masa
CONVERSIONS['kg'] = [1, operator.mul]
CONVERSIONS['g']  = [1e3, operator.mul]
# Momento de Fuerza
CONVERSIONS['N.m'] = [1, operator.mul]
CONVERSIONS['kgf'] = [0.1, operator.mul]
# Ninguna
CONVERSIONS[''] = [1.0, operator.mul]

def conversion(data,units):
	if not(units == ''):
		dataRet = []
		if units in ['C']:
			for i in range(len(data)):
				aux = (data[i][0],data[i][1]-CONVERSIONS[units])
				dataRet.append(aux)
		else:
			for i in range(len(data)):
				aux = (data[i][0],data[i][1]*CONVERSIONS[units])
				dataRet.append(aux)
		data = dataRet
	return data