#
# merging grackles cooling_rate.py and run_from_yt.py
# to hopefully calculate cooling rate in each cell in an Enzo simulations
# this will be fun?
#
import yt
import numpy as np

from itertools import chain

from pygrackle.api import *
from pygrackle.fluid_container import _yt_to_grackle

from pygrackle.grackle_wrapper import \
     calculate_cooling_time, \
     calculate_temperature, \
     chemistry_data, \
     solve_chemistry

from pygrackle.fluid_container import FluidContainer

from utilities.api import \
     setup_fluid_container, \
     set_cosmology_units, \
     get_cooling_units

from utilities.physical_constants import \
     mass_hydrogen_cgs, \
     sec_per_Gyr, \
     cm_per_mpc

def define_cooling(ds, grackle_data_file = "CloudyData_UVB=HM2012.h5",check=1):
    """
    Function loads grackle chemistry settings in order to compute the cooling 
    time and rate from a dataset. Function defines two fields, 'cooling_rate' and
    'cooling_time'.
    """

    # do some grackle-y things
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1

    # try to guess based on fields present within data set
    if 'H2II_Density' in ds.field_list:
        my_chemistry.primordial_chemistry = 2 # 9 species.. Joung 2012
    else:
        my_chemistry.primordial_chemistry = 1

    my_chemistry.metal_cooling = 1
    my_chemistry.grackle_data_file = grackle_data_file;

    # potentially useful parameters
#    my_chemistry.UVbackgroun
    my_chemistry.photoelectric_heating = 1

    # Set the units for the data set... this should be done carefully
    # and may not work universally
    my_chemistry.comoving_coordinates = 0
    my_chemistry.length_units         = ds.length_unit.convert_to_units('cm').value
    my_chemistry.density_units        = ds.mass_unit.convert_to_units('g').value / my_chemistry.length_units**3
    my_chemistry.time_units           = 1.0
    my_chemistry.velocity_units       = ds.velocity_unit.convert_to_units('cm/s').value
    my_chemistry.a_units              = 1.0

    # constant in cgs
    gravitational_constant = 6.67259E-8 * my_chemistry.density_units * my_chemistry.time_units**2

    a_value = 1.0

    my_chemistry.initialize(a_value)


    #
    # now that that is all setup, load in the data from yt
    #

    field_list = ds.field_list + ds.derived_field_list

    def _cooling_time(field, data):

        data_shape = np.shape(data['density'].value)

        fc_list = data_to_grackle(my_chemistry, field_list, data, update=False)

        cooling_time = [None]*len(fc_list)
        i = 0
        for fc in fc_list:

            # don't need the temperature calculation
            calculate_temperature(fc, a_value)

            # compute the cooling time
            calculate_cooling_time(fc, a_value)

            # get the cooling rate
            cooling_time[i] = np.abs(fc['cooling_time'])

            i = i + 1

        # list to array, and reshape
        cooling_time = list(chain.from_iterable(cooling_time))
        cooling_time = np.array(cooling_time).flatten()
        cooling_time = cooling_time.reshape(data_shape)


        return cooling_time * my_chemistry.time_units * yt.units.s

    def _cooling_rate(field, data):
        """
        Computes the cooling rate from the cooling time. Defined as energy loss rate
        density.
        """
        
        cooling_time = data['cooling_time']
    
        units = 'erg/s * cm**3'

        cooling_rate = data['thermal_energy'] / cooling_time * (data['Density']/(data['number_density'])**2)
    
        return cooling_rate.convert_to_units(units)

               
    def _cooling_time_dynamical_time_ratio(field, data):
        """
        returns ratio of cooling time to dynamical time
        """
 
        cooling_time = data['cooling_time'].convert_to_units('s')
        dynamical_time = data['dynamical_time'].convert_to_units('s')

        return cooling_time / dynamical_time

    # add the fields to the data set
    ds.add_field('cooling_time', function=_cooling_time, units='s',force_override=True)
    ds.add_field('cooling_rate', function=_cooling_rate, units='erg/s *cm**3',force_override=True)
    ds.add_field('CoolingTimeDynamicalTimeRatio',
                 function=_cooling_time_dynamical_time_ratio, units='', force_override=True)

##
