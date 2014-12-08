from yt.mods import *
import yt.analysis_modules.halo_profiler.api as HP
from yt.analysis_modules.halo_mass_function.api import *
import os.path


def load_hmf(pfname, pfdir = '.', type='co-moving'):
    """
    type = co-moving or volume 
    """
    
    pf = load(pfdir + "/" + pfname)
    
    if type == 'co-moving': # calculate the comoving volume cubed
        conv = pf['Mpc'] * pf.hubble_constant * (1.0 + pf.current_redshift)
    else:
        conv = 1.0

    hmf_data = np.genfromtxt(pfdir + "/" + pfname + "_hmf-haloes.dat", names = ['log10mass','mass','HMF'],skip_header=4,dtype=np.dtype(float))


    if os.path.isfile(pfdir + "/" + pfname + "_hmf-fit.dat"):
        hmf_fit = np.genfromtxt(pfdir + "/" + pfname + "_hmf-fit.dat", skip_header=5, names=['log10mass','mass','dHMF_fit','HMF_fit'], dtype=np.dtype(float))

#    print hmf_data['log10mass']

    hmf_dict = {'log10mass': hmf_data['log10mass'], 'mass': hmf_data['mass'],
                'HMF': hmf_data['HMF']*conv, 'dHMF_fit': hmf_fit['dHMF_fit'],
                'HMF_fit': hmf_fit['HMF_fit']*conv}    


    return hmf_dict


def full_halo_hmf(pfname, vo=200.0 , pfdir='.', overwrite=False):
    """ 
    Function operates in three stages:
    1) Runs the halo finder on (for now) the entire supplied dataset
    2) Runs the halo profiler a la the default 
    3) Calculates and does analytical fit for the halo mass function

    pfname => name of data file (e.x. 'RD0020')
    vo     => virial overdensity to use  Default 200.0 
    pfdir  => directory that contains pfname data file. Default '.' (current dir)
    overwrite => If the halo file exists already, the default is to load the halo
                 file instead of re-running the halo finder (since this takes time).
                 overwrite = True ignores this and recalculates, writing over prev.
    """

    pf = load(pfdir + "/" + pfname)
    halo_filename = "%s_halos.txt"%(pfname)

#    if os.path.isfile(halo_filename):
#        halos = LoadHaloes(pf, "%s_HopAnalysis"%(pfname))
    if overwrite or not os.path.isfile(halo_filename):
        halos = HaloFinder(pf)
        halos.write_out(halo_filename)
        halos.dump("%s_HopAnalysis"%(pfname))
    else:
        halos = LoadHaloes(pf, "%s_HopAnalysis"%(pfname))
    
    print "Now starting to run the halo profiler"
    # run the halo profiler and output
    hp = HP.HaloProfiler(pfdir + "/" + pfname, halo_list_file=halo_filename)
    print "adding filter"
    hp.add_halo_filter(HP.VirialFilter, must_be_virialized=True,
                       overdensity_field = 'ActualOverdensity',
                       virial_overdensity = vo,
                       virial_filters     = [['TotalMassMsun','>=','1e8']],
                       virial_quantities  = ['TotalMassMsun','RadiusMpc'])
    print "making profile"
    hp.make_profiles(filename="%s_VirialHalos.out"%(pfname))

    print "HMF"
    # run the halo mass function
    
    hmf = HaloMassFcn(pf, halo_file="%s_VirialHalos.out"%(pfname),
                      sigma8input=0.9, primordial_index=0.962, omega_baryon0=0.0490,
                      fitting_function=4, mass_column=6, num_sigma_bins=200)

    print "HMF out"
    hmf.write_out(prefix='%s_hmf'%(pfname))

    
def calculate_hmf(pf, pfname):
    hmf = HaloMassFcn(pf, halo_file="%s_VirialHalos.out"%(pfname),
                      sigma8input=0.9, primordial_index=0.962, omega_baryon0=0.0490,
                      fitting_function=4, mass_column=6, num_sigma_bins=200)

    print "HMF out"
    hmf.write_out(prefix='%s_hmf'%(pfname))

    return hmf
