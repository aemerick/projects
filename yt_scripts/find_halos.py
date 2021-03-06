from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
import os.path
import numpy as np

def find_halos(pf, center=[0.5,0.5,0.5], rc=[1.,1.,1.], lc=[0.,0.,0.],
               sub_vol=None,print_n_halos=5,output_dir=None,
               virial_overdensity=100.0):
    """
    Find halos and returns things. Defaults to using the hole box unless
    center, rc, and lc coordinates are provided, OR unless a specific
    subvolume is provided
    """

    if sub_vol == None:
        reg = pf.h.region(center,lc,rc)
    else:
        reg = sub_vol
        
    # generate list of halos
    halo_list = HaloFinder(pf,subvolume=reg)
    
    #
    if np.size(halo_list) < print_n_halos: 
        print_n_halos = np.size(halo_list)
       
    outDir = os.getcwd() + "/" +str(pf) +"_halos/"
    if not os.path.exists(outDir): os.makedirs(outDir)
  
    print "# using virial overdensity of ", virial_overdensity 
    print "# COM Total_mass maximum_radius virial_mass virial_radius"
    for halo in halo_list[0 : print_n_halos - 1]:
	    print halo.center_of_mass(), halo.total_mass(), halo.maximum_radius(), halo.virial_mass(virial_overdensity=virial_overdensity), pf['Mpc']*halo.virial_radius(virial_overdensity=virial_overdensity)
	    
	 
    halo_list.write_out(outDir + str(pf) + "_HopHalos.out")
    halo_list.dump(outDir + str(pf) + "_HopDump")
    
    return halo_list


def load_halos(filename):
    """
        Load a previously generated halos list from the given file
    """
    

   
