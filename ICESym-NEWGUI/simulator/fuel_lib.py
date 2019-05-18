import sys

def data_fuel(fuel):
    #
    if fuel == 'gasoline': # Ferguson
        alpha, beta, gamma, delta = 7.0, 17.0, 0.0, 0.0
        Q_LHV = 44.5e6
        h_vap = 310e3
        Afuel=[4.0652, 6.0977e-02, -1.8801e-05, 0.0, 0.0]
    elif fuel == 'diesel': # Ferguson
        alpha, beta, gamma, delta = 14.4, 24.9, 0.0, 0.0
        Q_LHV = 42.5e6
        h_vap = 230e3
        Afuel=[7.971, 1.1954e-01, -3.6858e-05, 0.0, 0.0]
    elif fuel == 'methane': # Ferguson
        alpha, beta, gamma, delta = 1.0, 4.0, 0.0, 0.0
        Q_LHV = 50.0e6
        h_vap = 509e3
        Afuel=[1.971324, 7.871586e-03, -1.048592e-06, 0.0, 0.0]
    elif fuel == 'methanol': # Ferguson
        alpha, beta, gamma, delta = 1.0, 4.0, 1.0, 0.0
        Q_LHV = 19.9e6
        h_vap = 1215e3
        Afuel=[1.779819, 1.262503e-02, -3.624890e-06, 0.0, 0.0]
    elif fuel == 'nitromethane': # Ferguson
        alpha, beta, gamma, delta = 1.0, 3.0, 2.0, 1.0
        Q_LHV = 0.0
        h_vap = 0.0
        Afuel=[1.412633, 2.087101e-02, -8.142134e-06, 0.0, 0.0]
    elif fuel == 'benzene': # Ferguson
        alpha, beta, gamma, delta = 6.0, 6.0, 0.0, 0.0
        Q_LHV = 40.2e6
        h_vap = 433e3
        Afuel=[-2.545087, 4.79554e-02, -2.030765e-05, 0.0, 0.0]
    elif fuel == 'toluene': # Raine
        alpha, beta, gamma, delta = 7.0, 8.0, 0.0, 0.0
        Q_LHV = 40.6e6
        h_vap = 412e3
        Afuel=[-2.09053, 5.654331e-2, -2.350992e-5, 0.0, 0.0]
    elif fuel == 'isooctane': # Raine
        alpha, beta, gamma, delta = 8.0, 18.0, 0.0, 0.0
        Q_LHV = 44.3e6
        h_vap = 308e3
        Afuel=[6.678e-1, 8.398e-2, -3.334e-5, 0.0, 0.0]
    elif fuel == 'methane_h': # Heywood
        alpha, beta, gamma, delta = 1.0, 4.0, 0.0, 0.0
        Q_LHV = 50.0e6
        h_vap = 509e3
        Afuel=[-1.4669e-01, 1.3248e-02, -5.3392e-06, 7.8785e-10, 8.3400e+04]
    elif fuel == 'propane': # Heywood
        alpha, beta, gamma, delta = 3.0, 8.0, 0.0, 0.0
        Q_LHV = 46.4e6
        h_vap = 426e3
        Afuel=[-7.4815e-01, 3.7409e-02, -1.9659e-05, 4.0531e-09, 6.1343e+03]
    elif fuel == 'hexane': # Heywood
        alpha, beta, gamma, delta = 6.0, 14.0, 0.0, 0.0
        Q_LHV = 44.14e6
        h_vap = 0.0
        Afuel=[-1.0456e+01, 1.0592e-01, -8.2592e-05, 2.6586e-08, 2.8500e+05]
    elif fuel == 'isooctane_h': # Heywood
        alpha, beta, gamma, delta = 8.0, 18.0, 0.0, 0.0
        Q_LHV = 44.3e6
        h_vap = 308e3
        Afuel=[-2.7835e-01, 9.1396e-02, -4.9209e-05, 1.0267e-08, -1.5575e+04]
    elif fuel == 'methanol_h': # Heywood
        alpha, beta, gamma, delta = 1.0, 4.0, 1.0, 0.0
        Q_LHV = 20.0e6
        h_vap = 1103e3
        Afuel=[-1.3617, 2.2227e-02, -1.3839e-05, 3.6329e-09, 1.0215e+05]
    elif fuel == 'ethanol': # Heywood
        alpha, beta, gamma, delta = 2.0, 6.0, 1.0, 0.0
        Q_LHV = 26.9e6
        h_vap = 840e3
        Afuel=[3.51756, 0.02, -0.00001, 0.0, 0.0]
    elif fuel == 'gasoline_h1': # Heywood
        alpha, beta, gamma, delta = 8.26, 15.5, 0.0, 0.0
        Q_LHV = 44.0e6
        h_vap = 350e3
        Afuel=[-1.2117e+01, 1.2914e-01, -1.0149e-04, 3.2584e-08, 2.9227e+05]
    elif fuel == 'gasoline_h2': # Heywood
        alpha, beta, gamma, delta = 7.76, 13.1, 0.0, 0.0
        Q_LHV = 44.0e6
        h_vap = 350e3
        Afuel=[-1.1323e+01, 1.1473e-01, -8.9202e-05, 2.8205e-08, 2.4381e+05]
    elif fuel == 'diesel_h': # Heywood
        alpha, beta, gamma, delta = 10.8, 18.7, 0.0, 0.0
        Q_LHV = 43.0e6
        h_vap = 250e3
        Afuel=[-4.5825, 1.2428e-01, -7.2334e-05, 1.6269e-08, 2.6067e+04]
    else:
        print 'Fuel %s is not included in the fuel library' % (fuel)
    
    return alpha, beta, gamma, delta, Q_LHV, h_vap, Afuel
