import yt
from absorbers.gas_plot_3d import *
import numpy as np

# stuff to make calculations from the absorber pickle....

# to do's:
# write the procedure to calculate the surface density of the absorbers



def loadAbsorbers(filepath):
    """
    load the absorbers pickle given the filepath to the directory containing
    the 'absorbers' folder that the pickle is located in
    """
    
    abs_pickle = filepath + "/absorbers/absorberList.pickle"
    
    if not os.path.isfile(abs_pickle):
        print "file does not exist, please generate the absorbers first"
    else:
        alist = loadAbsorberList(abs_pickle)
        
        
        
    return alist


def generateAbsorbers(pfname, pfdir, filepath, i = 0):
    """
    pfname -> name of data dump.... like "RD0020"
    pfdir  -> parent directory to dump ==>  pfdir + "/" + pfname + "/" pfname
                                            should be the data dump itself
    filepath -> path to absorption results... should contain observation result
                text files.
    """
    
    #i = 0
    halos = np.genfromtxt(pfdir + "/" + pfname + "_halos.txt", names=True)
    center = np.array([halos['x'][i],halos['y'][i],halos['z'][i]])
    R_vir  = halos['virial_radius'][i]
    
    ds = yt.load(pfdir + "/" + pfname + "/" + pfname)
    
    cluster = GalCluster(ds, center, R_vir)
    
    line_file = filepath + "/QSO_data_absorbers.out"
    data = np.genfromtxt(line_file, names=True)
    all_points = plot_gas_3D(cluster, filepath, data)
    
    if not os.path.isdir(filepath + "/absorbers"):
        os.mkdir(filepath + "/absorbers")
    
    saveAbsorberList(all_points, filepath + "/absorbers/absorberList.pickle")

    return all_points

def get_params_as_array(alist, xname, check_plot=False):
    """
    Given a list of absorbers, returns the desired parameter as an array
    
    if check_plot = True, this checks the 'plot' parameter of the given absorber.
    if false, it ignores the absorber
    """

    x = []
    
    for a in alist:   
        if not check_plot or (check_plot and a.plot):
    
            if xname=="pos":
                x.append(a.pos)
            else:
                x.append(a.param[xname])

    x = np.array(x)
    return x
    

def generate_histogram(alist, bin_field, bins, cut_field = None, cut_min = None, cut_max = None):
    """
    bins data by bin field into supplied bins. cut_field and cuts are allowed to be made to data
    """

    if bin_field == 'R':
        print 'when using R'
        print 'make sure units of bins are in cm!!!!'

        pos = get_params_as_array(alist, "pos")

        r = np.zeros(np.size(pos) / 3.0)
    
        #print 'pos'
        #print pos
        

        i = 0
        for p in pos:
            r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
            i = i + 1

        #print r
        #a = b
        binned_field_data = r 

    else:
        binned_field_data = get_params_as_array(alist, binned_field)

    
    if not cut_field == None:

        cut_field_data = get_params_as_array(alist, cut_field)

        binned_field_data = binned_field_data[ (cut_field_data<cut_max)*(cut_field_data>=cut_min) ]

    
    hist, bins = np.histogram(binned_field_data, bins=bins)

    return hist, bins
    
    
    
def calculateAbsorberDensity(pf, alist, R_vir, density_type, units = 'Mpc', bins = None, volume=None,
                                                  cut_field=None,
                                                  cut_min = None, cut_max = None):
    """
    density_type = 'area' or 'volume'
    length units are fixed to be Mpc...
    
    cut_field -> remove data points in the specified range 
    cmin, cmax -> lower and upper bound of range to cut points
    ---- can supply a list of cut_fields and cmin and cmax's for some more 
         complicated selections of the data set
    
    """
    for a in alist:
        a.plot = True    
    
    # make the cuts:
    if not cut_field == None:
        if not isinstance(cut_field, list):
            cut_field = [cut_field]
            cut_min,cut_max = [cut_min],[cut_max]
            
        for (cfield, cmin, cmax) in zip(cut_field,cut_min,cut_max):
            for a in alist:
                if a.param[cut_field] > cmin and a.param[cut_field] < cmax:
                    a.plot = False
    
    
    # pick the default bins depending on the density tupe
    if bins == None:
        if density_type == 'area':
            bins = np.arange(0.0, 3.6, 0.5)
        elif density_type == 'volume':
            bins = np.arange(0.0, 5.1, 0.5)
        
    
    # get position and calculate impact parameter
    pos = get_params_as_array(alist, "pos", check_plot = True)
     

    r = np.zeros(np.size(pos) / 3.0)
    b = np.zeros(np.size(pos) / 3.0)
    i = 0
    for p in pos:
        r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
        b[i] = (p[0]**2 + p[1]**2)**0.5
        i = i + 1
    
    # compute the histogram

    if units == 'Rvir':
        conv = pf['Mpc'] / pf['cm'] / R_vir
    elif units == 'Mpc':
        conv = pf['Mpc'] / pf['cm']
    elif units == 'cm':
        conv = 1.0
        
    r = r * conv
    b = b * conv
    
    center = (bins[:-1] + bins[1:])*0.5

    
    
    # normalize to number density of desired type
    if density_type == 'volume':
        hist, bins = np.histogram(r, bins = bins, density = False)
#        vol = 4.0/3.0 * np.pi * ((conv*bins[1:])**3 - (conv*bins[:-1])**3)
    
        print np.size(hist), np.size(volume)
        print hist
        print "center", center
        print "volume", volume
    
        num_density = hist / volume
    elif density_type == 'area':
        hist, bins = np.histogram(b, bins = bins, density = False)
        area = np.pi*((conv*bins[1:])**2 - (conv*bins[:-1])**2)
        num_density = hist / area
    
    
       
    
    # reset from the cuts
    for a in alist:
        a.plot = True
    
    return num_density, center
        
      
        
def obs_absorberDensity(RA,DEC, R_vir, density_type ='area', units = 'Mpc', bins = None,
                                                  cluster = 'virgo', cut_field=None,
                                                  cut_min = None, cut_max = None):
    """
    pass array of distances (projected) as 'r'. Use 'density_type'='volume'
    if 'r' is 3D distances... 
    
    pass an array to 'cut_field' to select points based upon cut_min and
    cut_max intervals. NOTE: this does not work like the function for the 
    simulation absorbers.
    """
    
    if cluster == 'virgo':
        clusterRA = 187.69708
        clusterDEC = 12.33694
        z          = 0.003600
#        R_vir = 1.6
    
    
    c = 3.0E5
    H_o = 70.0
    ang_sep = ((RA - clusterRA)**2 + (DEC - clusterDEC)**2)**0.5
    r       = (np.pi*ang_sep/180.0) * (1.0/(c*z*H_o))
    
    
    
    if not cut_field == None:
        r = r[ (cut_field < cut_min) * (cut_field > cut_max) ]


    if bins == None:
        if density_type == 'area':
            bins = np.arange(0.0, 3.6, 0.5)
        elif density_type == 'volume':
            bins = np.arange(0.0, 4.6, 0.5)
    
    if units == 'Rvir':
        conv = 1.0/R_vir
    elif units == 'Mpc':
        conv = 1.0


    center = (bins[:-1] + bins[1:])*0.5
    
    if density_type == 'volume':
        hist, bins = np.histogram(r, bins = bins, density = False)
        vol = 4.0/3.0 * np.pi * ((conv*bins[1:])**3 - (conv*bins[:-1])**3)
        num_density = hist / vol
    elif density_type == 'area':
        hist, bins = np.histogram(b, bins = bins, density = False)
        area = np.pi*((conv*bins[1:])**2 - (conv*bins[:-1])**2)
        num_density = hist / area
        
        
    return num_density, center
    
