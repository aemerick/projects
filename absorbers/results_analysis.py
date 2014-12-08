from yt.mods import *
import numpy as np
from scipy import optimize
import pyfits

CONST_c = 2.99e5 # km/s
m_H = 1.6733E-24
k_b = 1.380658E-16

powerlaw = lambda x,amp,index: amp * (x**index)

MULT_FACTOR = 1.0
print "WARNING USING ", MULT_FACTOR, " MULTIPLIER ON INTEGRATED HI"

class resultsObj:
    """
        asdf
    """
    
    def __init__(self, pfname, pfdir = '', i = 0, res_type='cluster'):
        self.pfname = pfname
        self.pf = load(pfdir + pfname + "/" + pfname)
        self.res_type = res_type

        if res_type == 'cluster':
            halos = np.genfromtxt(pfdir + pfname +"_halos.txt", names=True)
            self.cluster_center = np.array([halos['x'][i],halos['y'][i],halos['z'][i]])
            self.mass   = halos['virial_mass'][i]
            self.R_vir  = halos['virial_radius'][i] 
            self.ray_length = 5.0*self.R_vir # *****************************************

            

                
        elif res_type == 'igm':
            self.ray_length = self.pf['Mpc'] * 0.99
            
            
        self.z_center = self.ray_length/2.0 * self.pf.hubble_constant*100.0/CONST_c
        

        
        


        self.defineFieldLabels()

    def defineFieldLabels(self):
        self.fieldLabels = {'HI_N': r'N$_{\rm{HI}}$ (cm$^{-2}$)',
                            'integrated_HI_N': r'N$_{\rm{HI}}$ (cm$^{-2}$)'}



    def loadResults(self, loadDir):
    
        data_file = loadDir + "/QSO_data_fitted_lines.out"
        data      = np.genfromtxt(data_file,missing_values='nan',\
                                   filling_values=float('nan'), names=True)
                                   

        
        # make cuts for nonexistent lines
        
        isnan_array = np.array(np.where(\
                               np.array([not x for x in np.isnan(\
                               data['HI_N'])]) == True )).tolist()
               
        self.data = data[ isnan_array[0] ]
        self.total_abs = int(np.size(data['QSO_id']))
        
        # load the raw data from file
        raw_data_file = loadDir + "/QSO_data.out"
        self.raw_data      = np.genfromtxt(raw_data_file, names=True)
        self.total_lines     = int(np.size(self.raw_data['QSO_id']))
        self.data_dir = loadDir        
        ###
    
    
    def makeCuts(self, field, ltgt, val):
        """
        ltgt = 'lt' to remove things less than cut_val
        ltgt = 'gt' to remove things greater than cut_val
        """
        
        keep_array = np.array(np.where(self.data[field] >= val)).tolist()


        
        if ltgt == 'gt':
            keep_array = [not i for i in keep_array]
            
        keep_array = np.sort(keep_array[0])
        
        self.cut_data       = self.data[keep_array]
        
        print "Cuts made, removed all lines with " + field + " " + ltgt + " " + str(val)
    
    
        
    def listPossibleFields(self):
        try:
            print self.data.dtype.names
            print self.raw_data.dtype.names
        except:
            print "Need to load data first"
            
            
            
def cf(results, cftype, bins=None, rmin=0, rmax = np.inf, return_center = False):
    """
    Computes the covering fraction as a function of NHI for the provided bins.
    If no bins are provided, some default is used. Can also provide a rmin and
    rmax that will bin by radius as well.
    
    cftype = 'spectra'
    cftype = 'integrated'
    
    be warned... the way i coded this is suuuuper ineffecient...
    
    return_center = True returns centers of bins instead of bins
    """
    



    if bins == None:
        # make some default HIN bins
        nmin, nmax = 11.5, 21.0
        n = 25.0
        
        bins = np.linspace(nmin,nmax,(nmax-nmin)*n)
    
    if cftype == 'spectra':
        data = results.data
        NHI = data['HI_N']
    elif cftype == 'integrated':
        data = results.raw_data
        NHI = data['integrated_HI_N']*MULT_FACTOR
    
    impact_param = data['impact_parameter'] / results.R_vir
    
    
    NHI = np.log10(NHI)
    
    numbins = np.size(bins) - 1
    cf = np.zeros(numbins)
    prev_id = '0000'
    for j in np.arange(numbins):
        
        for i in np.arange(np.size(data['QSO_id'])):
            rho   = impact_param[i]
            qsoid = data['QSO_id'][i]
        
            if qsoid == prev_id:
                i = i + 1
        
            else:
            
                if any(NHI[ data['QSO_id'] == qsoid ] >= bins[j]) and rho <= rmax and rho >= rmin:
                
                    cf[j] += 1 # increment counter

                prev_id = qsoid   
                
        
    if return_center:
        bins = (bins[:-1] + bins[1:])*0.5
        

    ip = results.raw_data['impact_parameter'] / results.R_vir 
    numlines = np.size( ip[(ip >= rmin) * (ip <=rmax) ])


    cf = cf / (1.0 * numlines)

    return cf, bins

def mean_flux_decriment(results,use_fit=True):
    """
    spectra_type is either 'noise' or 'noiseless' depending on which one
    is desired to compute the mean flux decriment. If 'all' is supplied,
    it is computed for both noisy and noiseless spectra and saved to file
    with the given outDir
    """
    lmin = 1215.67 * (1.0+results.pf.current_redshift)
    lmax = ((results.ray_length*results.pf.hubble_constant*100.0)/CONST_c) * 1215.67 + lmin                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

    lmax = (((results.ray_length * results.pf.hubble_constant * (1.0 + results.pf.current_redshift))\
             *results.pf.hubble_constant*100.0)/CONST_c)*1215.67 + lmin
    
    lmin = lmin - ((500.0 / CONST_c)*1215.67)
    lmax = lmax + ((500.0 / CONST_c)*1215.67)
    #lmin, lmax = 3644, 3652.5
    
    print "using comoving distance!!!"
 #   print results.pf.current_redshift
 #   print results.pf.hubble_constant
 #   print results.ray_length
 #   print results.pf['Mpc']
    print lmin, lmax

    qso_id = results.raw_data['QSO_id']


    f = open(results.data_dir + "/flux_decriment.out", 'w')
    
    if use_fit:
        outstring = "%0004d %5.4e %5.4e %5.4e"
        f.write("#QSO_id noisy fit raw")
    else:
        outstring = "%0004d %5.4e"
        f.write("#QSO_id raw")
 
    D_noisy = 0.0
    D_fit   = 0.0
    D_raw   = 0.0

    counter = 0
    for id in qso_id:
        
        fullid = "_QSO-" + "%0004d"%(id)

        if counter % 1000 == 0:
            print fullid

        hdulist = pyfits.open(results.data_dir + "/" + fullid + "/spectrum"+fullid+".fits")
        tbdata = hdulist[1].data


        

           

     
        # load in the noisy, fitted, raw flux
   #     wavelength = data['wavelength']

        


        wavelength    = tbdata['wavelength']
        raw   = tbdata['flux']
 
        hdulist.close()

        wavelength_cut = (wavelength >= lmin)*(wavelength <= lmax)        



        raw   = raw[wavelength_cut]



        d_raw   = np.average(1.0 - raw)

        if use_fit:
            data = np.genfromtxt(results.data_dir + "/" + fullid + "/" + fullid + "_fit_data.out",
                                 names=True)
        
           
            fit   = data['fitted_flux']
            noisy = data['flux']
    
            fit   = fit[wavelength_cut]
            noisy = noisy[wavelength_cut]

            d_noisy = np.average(1.0 - noisy)
            d_fit   = np.average(1.0 - fit)

            f.write(outstring % (id, d_noisy, d_fit, d_raw))
            D_fit   += d_fit
            D_noisy += d_noisy
        else:
            f.write(outstring % (id, d_raw))



        D_raw   += d_raw

        counter += 1

    D_noisy = D_noisy / (1.0*results.total_lines)
    D_fit   = D_fit   / (1.0*results.total_lines)
    D_raw   = D_raw   / (1.0*results.total_lines)        

    
    print D_noisy, D_fit, D_raw

    if use_fit:
        f.write("#"+ outstring % (9999, D_noisy, D_fit, D_raw))
    else:
        f.write("#" + outstring % (9999, D_raw))
    
    f.close()
    
def calculateProbedVolume(results, outDir, radius_bins=None):
    """
    This function computes the volume along each line of sight for the cluster
    and sums it to calculate the total probed volume of the cluster. Binned
    by 0.5 Rvir. Outputs results to a text file
    """
    
    if radius_bins == None:
        radius_bins = np.arange(0.0,5.1,0.5)
    
    # make the array to store the results
    output = np.array(
                      [np.zeros(3 + 1 + np.size(radius_bins))
                                 for i in np.arange(results.total_lines)])
    
    # loop over the lines
    for j in np.arange(results.total_lines):
        
        R = results.raw_data['impact_parameter'][j] / results.R_vir
    
 #       print R
#        print R*results.R_vir
  #      apple = banana

        # load start and end points 
        start = np.array([results.raw_data['x0'][j], results.raw_data['y0'][j],
                                                     results.raw_data['z0'][j]])
        end   = np.array([results.raw_data['x1'][j], results.raw_data['y1'][j],
                                                     results.raw_data['z1'][j]])
                  
        los = end - start                                   
        # draw the ray and calculate total volume
        ray = results.pf.h.ray(start,end)
        cellVolume = ray['CellVolume'] * results.pf['Mpc']**3 / results.pf['cm']**3
        
        # in units of Mpc^3
        totalVolume = np.sum(cellVolume)
        
        output[j][3] = totalVolume
        
        first_time_through = True
        for i in np.arange(np.size(radius_bins) - 1):
            rmin = radius_bins[i]
            rmax = radius_bins[i+1]
            
            
            # if impact is in right range calculate violumes... otherwise skip
            if 0==0:
                if not first_time_through:
                    output[j][3 + i + 1] = np.sum( cellVolume[(dist<rmax)*(dist>=rmin)] )
                else:
                              
                    dist = np.zeros(np.size(ray['t']))
                    k = 0
                    for t in ray['t']:
                         
                        # in code units
                        pos = los * t  + start
                        
    #                    print pos
   #                     print results.cluster_center
  #                      print los
 #                       print start
#                        apple = pie
        
                        # in code units
                        dist[k] = ( (pos[0] - results.cluster_center[0])**2 +\
                                    (pos[1] - results.cluster_center[1])**2 +\
                                    (pos[2] - results.cluster_center[2])**2 )**0.5
                              
                        # now in units of R_vir
                        dist[k] = dist[k] * results.pf['Mpc'] / results.R_vir
                        
                        k = k + 1
                        
                    output[j][3 + i + 1] = np.sum( cellVolume[(dist<rmax)*(dist>=rmin)] )
                    first_time_through = False

                # calculate fractional volume
                output[j][3 + i + 1] = output[j][3 + i + 1] / totalVolume
                
        # do the area stuff here
        cellArea = cellVolume**(2.0/3.0)
        output[j][0] = np.min(cellArea)
        output[j][1] = np.max(cellArea)
        output[j][2] = np.average(cellArea)
        
        
    # now write out to a file until there is a better way to store...
    
    
    f = open(outDir + "/area_volume_sampled.txt",'w')
    
    f.write("# minArea maxArea avgArea totalVolume 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0\n")
    
    formatString = "%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n"
    
#    for j in np.arange(results.total_lines):
        
 #       f.write(formatString%(output[j][0], output[j][1], output[j][2],
  #                            output[j][3], output[j][4], output[j][5],
   #                           output[j][6], output[j][7], output[j][8],
    #                          output[j][9], output[j][10], output[j][11],
     #                         output[j][12], output[j][13]))
                              

    f.close()
    
    centers = 0.5*(radius_bins[1:] + radius_bins[:-1])
    sampledVolume = np.zeros(np.size(centers) + 1)
    
    for i in np.arange(np.size(centers)):
        for j in np.arange(results.total_lines):
            sampledVolume[i] += output[j][3]*output[j][3+1+i]
        
    sampledVolume[-1] = totalVolume
    # sampled volume in units of cubic MPC
    return sampledVolume, radius_bins, centers 

        

def f_N(results, Narray=None, dz=None, bins=None, fit=False):
    """
    pulls data from the results passed unless Narray and dz is specified.
    
    dz can either be an array or a single number. if single number, should be
    total redshift interval to normalize over.
    
    Narray is array of column densities
    
    use_dz = True computs dCounts / dz
    use_dN = True computs dCounts / dlog(N)
    both true gives    dCounts / (dz dlog(N) )
    
    if fit=True, returns slope of the distribution (beta)
    """
    
    if bins == None:
#        bins = np.log10( np.logspace(12,22,num=(22-12)*3.0) ) # default bins
        bins = np.logspace(12,22,num=(22-12)*3.0)   
    # calculate dz interval
    if dz == None:
        # assuming constant redshift interval lines of sight
        dz_line =  results.ray_length * results.pf.hubble_constant*100.0/CONST_c
    
        dz_total = results.total_lines * dz_line

    elif np.size(dz) == 1:
        dz_total = dz
    else:
        dz_total = np.sum(dz)
        
    # load column unless already loaded
    if Narray == None:
        Narray = results.data['HI_N']
        
    if np.average(Narray) < 100.0:
#        Narray = np.log10(Narray)
        Narray = 10**Narray
   
    # compute the plain old histogram
    hist, bins = np.histogram(Narray, bins = bins, density = False)
    
    # divide by total redshift interval   
    hist = hist / (1.0 * dz_total)

    # calculate bin centers and sizes
    centers = 0.5*(bins[:-1] + bins[1:])
    dN      = (bins[1:] - bins[:-1]) 
      
    # divide by bin sizes to get in terms of dN
    hist = hist / dN

    if fit:
        fitfunc = lambda p, x: p[0] + p[1]*x
        errfunc = lambda p, x, y: (y - fitfunc(p, x))
        
        sim_nmin, sim_nmax = 13.3, 14.5 # from Dave and pendelton
        sim_nmin, sim_nmax = 13.0, 14.0 #
        

#        obs_nmin, obs_nmax = 12.2, 17.5 # from Danforth et. al. 2014
        
        danforth_beta = 1.68 # +/- 0.02
        dave_beta = 1.70
        # adding lines from dave and 
        
        cut_array = (centers > 10**sim_nmin) * (centers < 10**sim_nmax) * (hist > 0)
#        print hist
#        print centers
#        print cut_array        

        # make the fit to dN/dzdN using unlogged bins
       
#        hfit, bins = np.histogram(Narray, bins = bins, density = False)
        
        
        hfit      = np.log10(hist[cut_array])
        cfit      = np.log10(centers[cut_array])
        
        # starting guesses        
        pinit = [0.01, - 1.7]

#        print hfit
 #       print cfit

        # make the fit        
        out = optimize.leastsq(errfunc, pinit, args=(cfit,hfit), full_output=1)
        
        pfinal = out[0]
        covar  = out[1]
        
        index  = pfinal[1]
        amp    = 10**pfinal[0]
        
        plaw = powerlaw(10**cfit, amp, index)

        # beta is defined as -1.0 times the power on the 
        # slope of the power law
 
        my_beta = -1.0*index
       
        return hist, bins, centers, amp, my_beta
    
    else:
        return hist, bins, centers
            
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    



