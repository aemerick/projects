from yt.mods import *

TEMPERATURE_MULT_CONSTANT = 1.0
print "WARNING TEMPERATURE_MULT_CONSTANT SET TO ", TEMPERATURE_MULT_CONSTANT, " in derived fields"

ICM_num_density_cut = 0.1

print "USING AN ICM DENSITY CUT OF ", ICM_num_density_cut 

def _HI_NumberDensity_2(field,data):
    return 2.0*data['HI_NumberDensity']

def _HI_Analytical_HM12(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'] / 1.0E4
    Gamma_12 = 0.0228

    n_H = data['HI_NumberDensity'] + data['HII_NumberDensity']

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12

    return n_HI

def _HI_Analytical_HM01(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'] / 1.0E4
    Gamma_12 = 0.10

    n_H = data['HI_NumberDensity'] + data['HII_NumberDensity']

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12

    return n_HI
    
def _HI_Analytical_Dave(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'] / 1.0E4
    Gamma_12 = 0.10 / 0.67

    n_H = data['HI_NumberDensity'] + data['HII_NumberDensity']

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12

    return n_HI


def _HI_Column_Density(field, data):
    l = data['CellVolume'] ** (1.0/3.0)
    return l*data['HI_NumberDensity']
    
def _Warm_Gas_Density(field, data):
    T   = data['Temperature']
    rho = data['Density']

    not_warm = (T < 2.0E4) #* (T > 1.0E6)

    rho[not_warm] = 0.0

    not_warm = T > 1.0E6
    
    rho[not_warm] = 0.0

    not_ICM = data['NumberDensity'] > ICM_num_density_cut

    rho[not_ICM] = 0.0

    return rho

def _Cold_Gas_Density(field,data):
    T = data['Temperature']
    rho = data['Density']

    not_cold = (T > 2.0E4)

    rho[not_cold] = 0.0

    not_ICM = data['NumberDensity'] > ICM_num_density_cut

    rho[not_ICM] = 0.0

    return rho

def _All_Gas_Density(field,data):
    not_ICM = data['NumberDensity'] > ICM_num_density_cut
    rho = data['Density']
    rho[not_ICM] = 0.0

    return rho

def _Hot_Gas_Density(field,data):
    T = data['Temperature']
    rho = data['Density']

    not_hot = T < 1.0E6

    rho[not_hot] = 0.0

    not_ICM = data['NumberDensity'] > ICM_num_density_cut

    rho[not_ICM] = 0.0

    return rho
    

def add(field):

    fDict = {'HI_Column_Density': _HI_Column_Density,
	      'Warm_Gas_Density': _Warm_Gas_Density,
               'Hot_Gas_Density': _Hot_Gas_Density,
              'Cold_Gas_Density': _Cold_Gas_Density,
              'All_Gas_Density': _All_Gas_Density, 
              'HI_Analytical_HM12': _HI_Analytical_HM12, 
              'HI_Analytical_HM01': _HI_Analytical_HM01,
              'HI_NumberDensity_2': _HI_NumberDensity_2}

    
    uDict = {'HI_Column_Density': r"\rm{cm}^{2}",
             'Warm_Gas_Density':  r"\rm{g} \rm{cm}^{-3}",
             'Cold_Gas_Density':  r"\rm{g} \rm{cm}^{-3}",
             'Hot_Gas_Density':   r"\rm{g} \rm{cm}^{-3}",
             'All_Gas_Density': r"\rm{g} \rm{cm}^{-3}",
             'HI_Analytical_HM12': r"\rm{cm}^{-3}",
             'HI_Analytical_HM01': r"\rm{cm}^{-3}",
             'HI_NumberDensity_2': r"\rm{cm}^{-3}"}
    
    
    add_field(field, function=fDict[field], units = uDict[field])
