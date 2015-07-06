import yt
import numpy as np

TEMPERATURE_MULT_CONSTANT = 1.0
print "WARNING TEMPERATURE_MULT_CONSTANT SET TO ", TEMPERATURE_MULT_CONSTANT, " in derived fields"

ICM_num_density_cut = 0.1

print "USING AN ICM DENSITY CUT OF ", ICM_num_density_cut 

def _CellVolume(field,data):
    return (data['dx'] * data['dy'] * data['dz']).convert_to_units('cm**3')

def _HI_NumberDensity_2(field,data):
    return 2.0*data['HI_NumberDensity']

def _HI_Analytical_HM96(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'].convert_to_units('K').value / 1.0E4
    Gamma_12 = 0.03

    n_H = (data['H_number_density'] + data['H_p1_number_density']).convert_to_units('cm**-3').value
    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12
    # do a fix for collisional ionization
    hot_gas = T4>=10.0
    n_HI[hot_gas] = n_H[hot_gas] *\
        (10.0**((13.9 - 5.4*np.log10(T4[hot_gas]*1.0E4) + 0.33*np.log10(T4[hot_gas]*1.0E4)**2)))


    return n_HI

def _HI_Analytical_HM12(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'].convert_to_units('K').value / 1.0E4
    Gamma_12 = 0.0228

    n_H = (data['H_number_density'] + data['H_p1_number_density']).convert_to_units('cm**-3').value
   

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12
    # do a fix for collisional ionization
    hot_gas = T4>=10.0
    n_HI[hot_gas] = n_H[hot_gas] *\
        (10.0**((13.9 - 5.4*np.log10(T4[hot_gas]*1.0E4) + 0.33*np.log10(T4[hot_gas]*1.0E4)**2)))
    return n_HI

def _HI_Analytical_HM01(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'].convert_to_units('K').value / 1.0E4
    Gamma_12 = 0.10

    n_H = (data['H_number_density'] + data['H_p1_number_density']).convert_to_units('cm**-3').value

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12
    # do a fix for collisional ionization
    hot_gas = T4>=10.0
    n_HI[hot_gas] = n_H[hot_gas] *\
        (10.0**((13.9 - 5.4*np.log10(T4[hot_gas]*1.0E4) + 0.33*np.log10(T4[hot_gas]*1.0E4)**2)))

    return n_HI
    
def _HI_Analytical_Dave(field, data):
    T4 = TEMPERATURE_MULT_CONSTANT*data['Temperature'].convert_to_units('K').value / 1.0E4
    Gamma_12 = 0.10 / 0.67

    n_H = (data['H_number_density'] + data['H_p1_number_density']).convert_to_units('cm**-3').value

    n_HI = 0.46 * n_H**2 * T4**(-0.76) / Gamma_12

    return n_HI


def _HI_Column_Density(field, data):
    l = data['dx'].convert_to_units('cm')
    return l*data['H_number_density']
    
def _Warm_Gas_Density(field, data):
    T   = data['Temperature']
    rho = data['Density'].convert_to_units('g/cm**3')

    not_warm = (T.value < 2.0E4) #* (T > 1.0E6)

    rho[not_warm] = 0.0 * yt.units.g / yt.units.cm**3

    not_warm = T.value > 1.0E6
    
    rho[not_warm] = 0.0 * yt.units.g / yt.units.cm**3

    not_ICM = data['number_density'].convert_to_units('cm**(-3)').value > ICM_num_density_cut

    rho[not_ICM] = 0.0 * yt.units.g / yt.units.cm**3

    return rho

def _Cold_Gas_Density(field,data):
    T = data['Temperature']
    rho = data['Density'].convert_to_units('g/cm**3')

    not_cold = (T.value > 2.0E4)

    rho[not_cold] = 0.0 * yt.units.g / yt.units.cm**3

    not_ICM = data['number_density'].convert_to_units('cm**(-3)').value > ICM_num_density_cut

    rho[not_ICM] = 0.0 * yt.units.g / yt.units.cm**3

    return rho

def _All_Gas_Density(field,data):
    not_ICM = data[('gas','number_density')].convert_to_units('cm**(-3)').value > ICM_num_density_cut
    rho = data['Density'].convert_to_units('g/cm**3')
    rho[not_ICM] = 0.0 * yt.units.g / yt.units.cm**3

    return rho

def _Hot_Gas_Density(field,data):
    T = data['Temperature']
    rho = data['Density'].convert_to_units('g/cm**3')

    not_hot = T.value < 1.0E6

    rho[not_hot] = 0.0 * yt.units.g / yt.units.cm**3

    not_ICM = data['number_density'].convert_to_units('cm**(-3)').value > ICM_num_density_cut

    rho[not_ICM] = 0.0 * yt.units.g / yt.units.cm**3

    return rho
    

def add(field):

    fDict = {'HI_Column_Density': _HI_Column_Density,
	      'Warm_Gas_Density': _Warm_Gas_Density,
               'Hot_Gas_Density': _Hot_Gas_Density,
              'Cold_Gas_Density': _Cold_Gas_Density,
              'All_Gas_Density': _All_Gas_Density, 
              'HI_Analytical_HM96': _HI_Analytical_HM96,
              'HI_Analytical_HM12': _HI_Analytical_HM12, 
              'HI_Analytical_HM01': _HI_Analytical_HM01,
              'HI_NumberDensity_2': _HI_NumberDensity_2,
              'CellVolume'        : _CellVolume}

    
    display_name_dict = {'HI_Column_Density': r"N$_{\rm{HI}}$",
             'Warm_Gas_Density':  r"\rm{g} \rm{cm}^{-3}",
             'Cold_Gas_Density':  r"\rm{g} \rm{cm}^{-3}",
             'Hot_Gas_Density':   r"\rm{g} \rm{cm}^{-3}",
             'All_Gas_Density': r"\rm{g} \rm{cm}^{-3}",
             'HI_Analytical_HM96': r"\rm{cm}^{-3}",
             'HI_Analytical_HM12': r"\rm{cm}^{-3}",
             'HI_Analytical_HM01': r"\rm{cm}^{-3}",
             'HI_NumberDensity_2': r"\rm{cm}^{-3}",
             'CellVolume': 'Volume'}
    
    units_dict = {'HI_Column_Density': "cm**(-2)",
                  "Warm_Gas_Density" : "g/cm**3",
                  "Cold_Gas_Density" : "g/cm**3",
                  "Hot_Gas_Density"  : "g/cm**3",
                  "All_Gas_Density"  : "g/cm**3",
                  "HI_Analytical_HM96" : "cm**(-3)",
                  "HI_Analytical_HM12" : "cm**(-3)",
                  "HI_Analytical_HM01" : "cm**(-3)",
                  "HI_NumberDensity_2" : "cm**(-3)",
                  "CellVolume"         : "cm**(3)"}

    
    yt.add_field(field, function=fDict[field], display_name = display_name_dict[field], units=units_dict[field])
