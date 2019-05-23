#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:40:34 2019

@author: etekken
"""

from numpy import trapz
from math import pi
from utils import show_message
from units import UNITS, CONVERSIONS

GROUP_A = ['IMEP per Cylinder','FMEP per Cylinder','BMEP per Cylinder']
GROUP_B = ['Power Indicated','Power Effective','Torque Indicated','Torque Effective','Mechanical Efficiency']

class GeneralAttributes():
    def __init__(self, run_attributes, read_normal_txt, read_extras_txt, cylinders, rho, current_test_dir):
        self.general_atributes = {}
        self.run_attributes = run_attributes
        self.read_extras_txt = read_extras_txt
        self.read_normal_txt = read_normal_txt

        self.rpms = [int(irpm) for irpm in self.run_attributes['rpms']]
        self.ncycles = self.run_attributes['ncycles']
        self.ncyls = len(cylinders)
        self.cylinders = cylinders
        self.rho = rho
        self.current_test_dir = current_test_dir
        return
    
    def return_calculated_variable(self, variable, icycle, label, unit):        
        if variable not in self.general_atributes.keys():
            if variable in GROUP_A:
                self.calculate_IMEP_FMEP_BMEP()
            elif variable in GROUP_B:
                self.calculate_torque_and_power()
            else:
                self.calculate_SFC_fuel_conv_vol_eff()

        datas   = []
        legends = []
        scale   = CONVERSIONS[unit]
        if variable in self.general_atributes.keys():
            if 'per Cylinder' in variable:
                for icylinder in range(self.ncyls):
                    data = [ [float(irpm),self.general_atributes[variable][icylinder][irpm][icycle]*scale] for irpm in self.rpms ]
                    legends.append(label+"_"+variable+"_Cylinder_"+str(icylinder)+"_Cycle_"+str(icycle))
                    datas.append(data)
            else:
                data = [ [float(irpm),self.general_atributes[variable][icycle][irpm]*scale] for irpm in self.rpms ]
                legends.append(label+"_"+variable+"_Cycle_"+str(icycle))
                datas.append(data)
        else:
            show_message('%s not yet calculated. Aborting'%variable)
        return [datas,legends]
    
    def calculate_work_vd_pMax(self):
        w       = {}    # (ncyls x nrpms x ncycles)
        volDesp = {}    # (ncyls x nrpms x ncycles)
        pMax    = {}    # (ncyls x nrpms x ncycles)

        for icylinder in range(self.ncyls):
            w[icylinder]        = {}
            volDesp[icylinder]  = {}
            pMax[icylinder]     = {}
            for irpm in self.rpms:
                w[icylinder][irpm]          = {}
                volDesp[icylinder][irpm]    = {}
                pMax[icylinder][irpm]       = {}
                for icycle in range(1,self.ncycles+1):
    				# Volumen
                    archive = "%s/RPM_%s/cyl_extras_%s.txt"%(self.current_test_dir,irpm,icylinder)
                    volData = self.read_extras_txt(archive,icycle,'Volume','Cylinders','m^3')
    				# Presion
                    archive = "%s/RPM_%s/cyl_%s.txt"%(self.current_test_dir,irpm,icylinder)
                    pData = self.read_normal_txt(archive,0,icycle,1,'Pa')

                    work    = []
                    vd      = []
                    pMaxAux = []

                    if not len(volData) or not len(pData):
                        show_message('There was an error loading the volume and presurre data')
                        return False

                    # primeros valores
                    firstV = (volData[0][1] + volData[-1][1])/2
                    firstp = (pData[0][1] + pData[-1][1])/2
                    dV = [firstV]
                    p = [firstp]

                    # valores intermedios
                    kk = len(volData)
                    if len(pData) < len(volData):
                        kk = len(pData)
                    for k in range(kk):
                        dV.append(volData[k][1])
                        p.append(pData[k][1])

                    # ultimos valores
                    dV.append(firstV)
                    p.append(firstp)
                    vMax = max(dV)
                    vMin = min(dV)

                    vd = (vMax-vMin)    # volumen desplazado en el ciclo	    
                    work = trapz(p,dV)  # trabajo efecutuado en el ciclo
                    pMaxAux = max(p)    # presion maxima en el ciclo
    				
                    w[icylinder][irpm][icycle] = work
                    volDesp[icylinder][irpm][icycle] = vd
                    pMax[icylinder][irpm][icycle] = pMaxAux

        self.general_atributes['w']         = w
        self.general_atributes['volDesp']   = volDesp
        self.general_atributes['pMax']      = pMax
        return True

    def calculate_IMEP_FMEP_BMEP(self):
        # calculo primero estas variables auxiliares
        if (not self.general_atributes.keys()) or\
        not (all(i in ['w','volDesp','pMax'] for i in self.general_atributes.keys())):
            if not self.calculate_work_vd_pMax():
                show_message('There was an error calculating the variable')
                return False

        	# calculo de IMEP_per_cylinder, FMEP_per_cylinder y BMEP_per_cylinder (ncyls x ncycles x nrpms)
        c1 = 0.4       # bar
        c2 = 0.005
        c3 = 0.09      # bar/(m/s)
        c4 = 9e-4      # bar/(m/s)^2
        Bar2Pa = 1e5
        IMEP_per_cylinder = {}
        FMEP_per_cylinder = {}
        BMEP_per_cylinder = {}

        for icylinder in range(self.ncyls):
            IMEP_per_cylinder[icylinder] = {}
            FMEP_per_cylinder[icylinder] = {}
            BMEP_per_cylinder[icylinder] = {}
            for irpm in self.rpms:
                IMEP_per_cylinder[icylinder][irpm] = {}
                FMEP_per_cylinder[icylinder][irpm] = {}
                BMEP_per_cylinder[icylinder][irpm] = {}
                sm = float(irpm)*2.0*self.cylinders[icylinder].object['crank_radius']/30.0
                for icycle in range(1,self.ncycles+1):
                    aux1 = self.general_atributes['w'][icylinder][irpm][icycle]/self.general_atributes['volDesp'][icylinder][irpm][icycle]
                    aux2 = (c1 + c2*self.general_atributes['pMax'][icylinder][irpm][icycle]/Bar2Pa + c3*sm + c4*sm*sm) * Bar2Pa
                    aux3 = aux1 - aux2
                    IMEP_per_cylinder[icylinder][irpm][icycle] = aux1
                    FMEP_per_cylinder[icylinder][irpm][icycle] = aux2
                    BMEP_per_cylinder[icylinder][irpm][icycle] = aux3

        self.general_atributes['IMEP per Cylinder'] = IMEP_per_cylinder
        self.general_atributes['FMEP per Cylinder'] = FMEP_per_cylinder
        self.general_atributes['BMEP per Cylinder'] = BMEP_per_cylinder
        return True
    
    def calculate_torque_and_power(self):
        # calculo de Potencia y Torque (efectivos e indicados) (ncycles x nrpms)

        if (not self.general_atributes.keys()) or\
        not (all(i in ['IMEP per Cylinder','FMEP per Cylinder','BMEP per Cylinder'] for i in self.general_atributes.keys())):
            if not self.calculate_IMEP_FMEP_BMEP():
                show_message('There was an error calculating the variable')
                return False

        power_indicated         = {}
        torque_indicated        = {}
        power_effective         = {}
        torque_effective        = {}
        mechanical_efficiency   = {}
        nstroke = self.run_attributes['nstroke']
        nr = nstroke/2

        for icycle in range(1,self.ncycles+1):
            power_indicated[icycle]         = {}
            torque_indicated[icycle]        = {}
            power_effective[icycle]         = {}
            torque_effective[icycle]        = {}
            mechanical_efficiency[icycle]   = {}
            for irpm in self.rpms:
                N = irpm/60.0
                aux1 = 0
                aux3 = 0
                for icylinder in range(self.ncyls):
                    Vd = self.general_atributes['volDesp'][icylinder][irpm][icycle]
                    aux1 = aux1 + self.general_atributes['IMEP per Cylinder'][icylinder][irpm][icycle]*Vd*N/nr
                    aux3 = aux3 + self.general_atributes['BMEP per Cylinder'][icylinder][irpm][icycle]*Vd*N/nr
                aux4 = aux3/(2*pi*N)
                aux2 = aux1/(2*pi*N)
                aux5 = aux3/aux1
                power_indicated[icycle][irpm]       = aux1
                torque_indicated[icycle][irpm]      = aux2
                power_effective[icycle][irpm]       = aux3
                torque_effective[icycle][irpm]      = aux4
                mechanical_efficiency[icycle][irpm] = aux5

        self.general_atributes['Power Indicated']       = power_indicated
        self.general_atributes['Power Effective']       = power_effective
        self.general_atributes['Torque Indicated']      = torque_indicated
        self.general_atributes['Torque Effective']      = torque_effective
        self.general_atributes['Mechanical Efficiency'] = mechanical_efficiency
        return True
    
    def calculate_SFC_fuel_conv_vol_eff(self):
        # Calculo de consumo especifico y rendimientos mecanicos, volumetricos y de conversion de combustible
        
        if (not self.general_atributes.keys()) or\
        not (all(i in ['Power Indicated','Power Effective','Torque Indicated','Torque Effective','Mechanical Efficiency'] for i in self.general_atributes.keys())):
            if not self.calculate_torque_and_power():
                show_message('There was an error calculating the variable')
                return False
        
        (mfc,mair)  = self.getMasses() if 'mfc' not in self.general_atributes.keys() else \
                    (self.general_atributes['mfc'] ,self.general_atributes['mair'])
        nstroke = self.run_attributes['nstroke']
        nr = nstroke/2
        Q_fuel = self.cylinders[0].object['Q_fuel']

        SFC_indicated                           = {}
        SFC_effective                           = {}
        fuel_conversion_efficiency_indicated    = {}
        fuel_conversion_efficiency_effective    = {}

        for icycle in range(1,self.ncycles+1):
            SFC_indicated[icycle] = {}
            SFC_effective[icycle] = {}
            fuel_conversion_efficiency_indicated[icycle] = {}
            fuel_conversion_efficiency_effective[icycle] = {}
            for irpm in self.rpms:
                N = irpm/60.0
                mfcTotal = 0
                for icylinder in range(self.ncyls):
                    mfcTotal = mfcTotal + mfc[icylinder][irpm][icycle]
                aux1 = mfcTotal*N/(nr*self.general_atributes['Power Indicated'][icycle][irpm])
                aux2 = mfcTotal*N/(nr*self.general_atributes['Power Effective'][icycle][irpm])
                aux4 = 1/(aux1*Q_fuel)
                aux5 = 1/(aux2*Q_fuel)
                SFC_indicated[icycle][irpm] = aux1
                SFC_effective[icycle][irpm] = aux2
                fuel_conversion_efficiency_indicated[icycle][irpm] = aux4
                fuel_conversion_efficiency_effective[icycle][irpm] = aux5

        self.general_atributes['SFC Indicated'] = SFC_indicated
        self.general_atributes['SFC Effective'] = SFC_effective
        self.general_atributes['Fuel Conversion Efficiency Indicated'] = fuel_conversion_efficiency_indicated
        self.general_atributes['Fuel Conversion Efficiency Effective'] = fuel_conversion_efficiency_effective

        volumetric_efficiency = {}
        for icylinder in range(self.ncyls):
            volumetric_efficiency[icylinder] = {}
            for irpm in self.rpms:
                volumetric_efficiency[icylinder][irpm] = {}
                for icycle in range(1,self.ncycles+1):
                    Vd = self.general_atributes['volDesp'][icylinder][irpm][icycle]
                    aux = mair[icylinder][irpm][icycle]/(self.rho*Vd)
                    volumetric_efficiency[icylinder][irpm][icycle] = aux

        self.general_atributes['Volumetric Efficiency per Cylinder'] = volumetric_efficiency

        volumetric_efficiency_global = {}
        IMEP_global = {}
        BMEP_global = {}
        FMEP_global = {}

        for icycle in range(1,self.ncycles+1):
            volumetric_efficiency_global[icycle]    = {}
            IMEP_global[icycle]                     = {}
            BMEP_global[icycle]                     = {}
            FMEP_global[icycle]                     = {}
            for irpm in self.rpms:
                aux1 = 0
                aux2 = 0
                aux3 = 0
                aux4 = 0
                active_cyls = self.ncyls
                for icylinder in range(self.ncyls):
                    aux1 = aux1 + self.general_atributes['Volumetric Efficiency per Cylinder'][icylinder][irpm][icycle]
                    aux2 = aux2 + self.general_atributes['IMEP per Cylinder'][icylinder][irpm][icycle]
                    aux3 = aux3 + self.general_atributes['BMEP per Cylinder'][icylinder][irpm][icycle]
                    aux4 = aux4 + self.general_atributes['FMEP per Cylinder'][icylinder][irpm][icycle]
                    if self.general_atributes['IMEP per Cylinder'][icylinder][irpm][icycle] <= 0.0:
    					active_cyls -= 1
                aux1 = aux1/active_cyls
                aux2 = aux2/active_cyls
                aux3 = aux3/active_cyls
                aux4 = aux4/active_cyls
                volumetric_efficiency_global[icycle][irpm] = aux1
                IMEP_global[icycle][irpm] = aux2
                BMEP_global[icycle][irpm] = aux3
                FMEP_global[icycle][irpm] = aux4

        self.general_atributes['Volumetric Efficiency Global'] = volumetric_efficiency_global
        self.general_atributes['IMEP Global'] = IMEP_global
        self.general_atributes['FMEP Global'] = FMEP_global
        self.general_atributes['BMEP Global'] = BMEP_global
        return
    
    def getMasses(self):
        mfc  = {}
        mair = {}
        for icylinder in range(self.ncyls):
            mfc[icylinder]  = {}
            mair[icylinder] = {}
            angleClose = self.cylinders[icylinder].object['angleClose']
            for irpm in self.rpms:
                mfc[icylinder][irpm]  = {}
                mair[icylinder][irpm] = {}
                rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
                archive = rpm_folder + "/cyl_extras_" + str(icylinder) + ".txt"
                for icycle in range(1,self.ncycles+1):
                    ms = []
                    for ivariable in ('Mass of Fuel','Mass of Air'):
                        data = self.read_extras_txt(archive,icycle,ivariable,'Cylinders','kg')
                        ndata = len(data)
                        mAux = [ data[i][1] for i in range(0,ndata) if (data[i][0]==angleClose)]
                        # Si no hay uno que sea igual al angulo, interpolar con el anterior y proximo
                        if not len(mAux):
                            ia = 0
                            # En algunos casos, el angulo del extras de cilindro no comienza
                            # desde cero, por lo que es necesario llevar a ia al cero
                            while data[ia][0]>angleClose:
                                ia+=1
                            while data[ia][0]<angleClose:
                                mAux_ant = data[ia][1]
                                ia+=1
                            mAux_pos = data[ia][1]
                            ms.append((mAux_ant+mAux_pos)/2.0)
                        else:
                            ms.append(mAux[0])
                    mfc[icylinder][irpm][icycle]    = ms[0]
                    mair[icylinder][irpm][icycle]   = ms[1]

        self.general_atributes['mfc']   = mfc
        self.general_atributes['mair']  = mair
        return (mfc,mair)