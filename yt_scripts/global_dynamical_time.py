import yt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad



def define_dynamical_time(ds, r, m, halo_center):
    """
    Given the radial profile of the dark matter halo density, defines
    the dynamical time.
    """
       
    # allow user to give arbitrary units... if none given assume cgs
    if hasattr(r, 'value'):
        r = r.value
        r = r.units
    else:
        r_units = yt.units.cm

    if hasattr(m, 'value'):
        m = m.value
        m_units = m.units
    else:
        m_units = yt.units.g

#    density_profile_function = UnivariateSpline(r, rho)

#    integral = lambda x : 4.0 * np.pi * x * x * density_profile_function(x)       

#    mass_profile_integral_function    = lambda x : quad(integral, 0.0, x)

#    sample_points = np.linspace(np.min(r), np.max(r), 1000.0)

#    mass = np.zeros(np.size(sample_points))
    
#    for i in np.arange(np.size(sample_points)):
#        mass[i] = mass_profile_integral_function(sample_points[i])[0]

#    mass_profile_function = UnivariateSpline(sample_points, mass)
    mass_profile_function = UnivariateSpline(r, m)

    def _global_dynamical_time(field,data):
        # get radial distance from center
        x = data['x']
        y = data['y']
        z = data['z']

        R = ( (halo_center[0]-x)**2 + (halo_center[1]-y)**2 +\
              (halo_center[2]-z)**2)**0.5
       
#        R = R.convert_to_units(r_units)
        R = R.convert_to_units('cm')
 
        M = mass_profile_function(R.value) * m_units
        average_density = M / (4.0 * np.pi * R**3 / 3.0)
        average_density = average_density.value * m_units / r_units**3

        tau = (3.0 * np.pi / (32.0 * yt.physical_constants.G * average_density))**0.5
 
        return tau.convert_to_units('s')


    ds.add_field(('gas','global_dynamical_time'),\
                   function=_global_dynamical_time, units='s',force_override=True)

