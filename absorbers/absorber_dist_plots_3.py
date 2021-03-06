from absorbers.gas_plot_3d import *
import numpy as np
import yt
from matplotlib import rc
from plotting import plotTools as my_pt

from scipy.stats import ks_2samp

fsize = 17.5
linewidth = 2.0

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fsize)

pfname = "RD0020"
ORIENTATION = 'z'
#center     = np.array([0.49730353,  0.50146712,  0.5017492]) # in code units
#R_vir  = 1.85 # Mpc
overwrite_list = False
scheme_vh_vpec_1D = 'N_bins'
reverse_lists = True
#to_dir = "SP/10x10/b_30/"
#to_dir = "SP/50x50/z_5Rvir/"
#to_dir = "SP/50x50/absorber_testing/"
to_dir = "SP/50x50/new_abs/"
abs_dir = to_dir + "absorbers/"
pec_vel=True


halos = np.genfromtxt(pfname + "_halos.txt", names=True)
i = 0
center = np.array([halos['x'][i],halos['y'][i],halos['z'][i]])
mass   = halos['virial_mass'][i]
R_vir  = halos['virial_radius'][i]

CONSTANT_m_H = 1.67493E-24
CONSTANT_k   = 1.38066E-16



#if pfname == "RD0020":

#    center     = np.array([0.49730353,  0.50146712,  0.5017492] )# in code units
#    R_vir      = 1.85 # in MPC
#elif pfname == "RD0013":
#    center = np.array([ 0.49571984,  0.51244465 , 0.51293379])
#    R_vir = 0.663835303192

# load pf
ds = yt.load(pfname +"/"+pfname)


cluster = GalCluster(ds,center,R_vir)
# load generated list of absorbers
if not os.path.isfile(abs_dir + "absorberList.pickle") or overwrite_list:
    line_file_name = to_dir + "QSO_data_absorbers.out"

    data = np.genfromtxt(line_file_name,names=True)
    all_points = plot_gas_3D(cluster,to_dir,data,pec_vel=pec_vel)
    saveAbsorberList(all_points, abs_dir + "absorberList.pickle")

alist = loadAbsorberList(abs_dir+"absorberList.pickle")

# some globals (ew.. I know)
point_size = 15
cmap       = 'spectral'


def get_params_as_array(alist,xname):

    x = []
    
    for a in alist:
        if xname=="pos":
            x.append(a.pos)
        else:
            x.append(a.param[xname])

    x = np.array(x)
    return x


def labeldict():
    """
    """
    
    labels = {"HI_Column_Density": r'log N$_{HI}$ cm$^{-2}$',
              "Metallicity": r'log Z/Z$_{sun}$',
              "Temperature": r'log T K',
              "vh": r"v$_{Hubble}$ km/s",
              "vpec": r"v$_{pec}$ km/s",
              "r" : r"R (R$_{vir}$)"}
              
    return labels

def plot_N_Z(alist,abs_dir):
    """
    """
    labels = labeldict()
    
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    N, Z = [],[]
    for a in alist:
        N.append(a.param["HI_Column_Density"])
        Z.append(a.param["Metallicity"])
    
    pos = get_params_as_array(alist,"pos")
    N   = get_params_as_array(alist,"HI_Column_Density")
    
    
    r = np.zeros(np.size(N))
    i = 0
    for p in pos:
        r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
        i = i + 1
    
    N = np.array(N)
    Z = np.array(Z)
#    N = np.log10(N)
 #   Z = np.log10(Z)
    
    ax1.scatter(N,Z, s = point_size)
    
    ax1.set_ylim(np.min(Z),np.max(Z))
    
    ax1.set_xlabel(labels["HI_Column_Density"])
    ax1.set_ylabel(labels["Metallicity"])
    ax1.loglog()
    #plt.tight_layout()
    fig.savefig(abs_dir + "N_Z_dist.png")
    plt.close(fig)

    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    Nbins = [ [12.0,13.0] , [13.0,14.0], [14.0,15.0] , [15.0,17.0] ,[17.0,25.0] ]

    colors = ['black','blue', 'green','orange', 'red']
    #colors = ['black'] * 4
    psizes  = [15.0, 30.0, 75.0,  115.0, 165.0]
    alphas  = [1.0, 0.9, 0.8, 0.70, .6]
       
    if reverse_lists:
        Nbins.reverse()
        colors.reverse()
        psizes.reverse()
        #alphas.reverse()

    j = 0


     
    for Nlim in Nbins:

        selection = ( np.log10(N)  >= Nlim[0])*( np.log10(N) <=Nlim[1])
        ax1.scatter((r[selection]*yt.units.cm).convert_to_units('Mpc').value /R_vir, np.log10(Z[selection]), s = psizes[j], c=colors[j],alpha=alphas[j])
        j = j + 1
    
    #ax1.scatter(r*pf['Mpc']/pf['cm']/R_vir, Z, s = point_size)
    #ax1.semilogy()
    ax1.set_ylim(-9,1)
    ax1.minorticks_on()
    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(labels['Metallicity'])
    fig.savefig(abs_dir + "metallicity_position.png")
    plt.close(fig)


    ##################################################
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)

    mappable = ax1.scatter(N,Z, s = point_size,c=(r*yt.units.cm).convert_to_units('Mpc').value / R_vir, cmap='Blues')
    
    ax1.set_ylim(np.min(Z),np.max(Z))
    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["r"])
    ax1.set_xlabel(labels["HI_Column_Density"])
    ax1.set_ylabel(labels["Metallicity"])
    ax1.loglog()
#    ax1.set_ylim(1E-9,10)
    #plt.tight_layout()
    fig.savefig(abs_dir + "N_Z_r_dist.png")
    plt.close(fig)
    
    
def plot_vH_N(alist,abs_dir):
    """
    """
    labels=labeldict()
    
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    vh = get_params_as_array(alist,"vh")
    N  = get_params_as_array(alist,"HI_Column_Density") 
   
   
    # make sure vH is normalized...
    if np.min(vh) >= 0.0:
        mid = 0.5*(np.min(vh) + np.max(vh))
        
        vh = vh - mid
        
    
    ax1.scatter(vh,N,s=point_size)
    ax1.semilogy()
    ax1.set_xlabel(labels['vh'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    fig.savefig(abs_dir + "vh_N_dist.png")
    plt.close(fig)
    
def plot_vel_dir(alist, abs_dir):
    velocity = get_params_as_array(alist,"Velocity")
    pos = get_params_as_array(alist,"pos")
    
    v = np.zeros(np.size(N))
    i=0
    
    
    
    for (vel,p) in zip(velocity,pos):
        r = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
        v[i] = vel*p/r
        
        i = i + 1
    
    v = (v * yt.units.cm).convert_to_units('km').value #pf['km'] / pf['cm'] # now in km/s
    
    hist, bins = np.hist(v, bins = 25, density =False)
    
    centers = 0.5*(bins[:-1] + bins[1:])
    
    plt.plot(centers, hist)
    plt.xlabel('vel projected towards cluster center')
    plt.close()
    
    
def plot_vH_vpec(alist,abs_dir):
    """
    """
    labels=labeldict()
    

    
    xname = 'vh'
    yname = 'vpec'
    
    vh = get_params_as_array(alist,xname)
    vpec = get_params_as_array(alist,yname)
    velocity = get_params_as_array(alist,"Velocity")
    N = get_params_as_array(alist,"HI_Column_Density")
    
    # make sure vH is normalized...
    if np.min(vh) >= 0.0:
        mid = 0.5*(np.min(vh) + np.max(vh))
        
        vh = vh - mid
      
   # print velocity  
    v = np.zeros(np.size(N))
    i=0
    for vel in velocity:
        v[i] = ((vel[0]**2 + vel[1]**2 + vel[2]**2))**(0.5)
        v[i] = (v[i] * yt.units.cm).convert_to_units('km').value# pf['km'] / pf['cm']
        i = i + 1
        
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    ax1.scatter(vh,vpec,s=point_size)

    ax1.set_xlabel(labels[xname])
    ax1.set_ylabel(labels[yname])
    plt.tight_layout()
    fig.savefig(abs_dir + xname +"_"+yname+"_dist.png")
    plt.close(fig)
    
    fig = plt.figure(figsize = [6.5,6.5])
    ax1 = fig.add_subplot(111)

    if scheme_vh_vpec_1D == 'N_bins':

        Nbins = [ [12.0,13.0] , [13.0,14.0], [14.0,15.0] , [15.0,17.0] ,[17.0,25.0] ]

        colors = ['black','blue', 'green','orange', 'red']
        #colors = ['black'] * 4
        psizes  = [15.0, 30.0, 75.0,  115.0, 165.0]
        alphas  = [1.0, 0.9, 0.8, 0.70, .6]
       
        if reverse_lists:
            Nbins.reverse()
            colors.reverse()
            psizes.reverse()
            #alphas.reverse()

        j = 0

        for Nlim in Nbins:

            selection = ( np.log10(N)  >= Nlim[0])*( np.log10(N) <=Nlim[1])
            ax1.scatter(vh[selection],vpec[selection], s = psizes[j], c=colors[j],alpha=alphas[j])
            j = j + 1

    else:
        ax1.scatter(vh,vpec,s=point_size)

    ax1.set_xlabel(r'v$_{\rm{H}}$ (km/s)')
    ax1.set_ylabel(r'v$_{\rm{pec}}$ (km/s)')

#    ax1.set_aspect(1.)
    ax1.set_ylim(-2250,2250)
    ax1.set_xlim(-600,600)
    ax1.minorticks_on()
#    ax1.set_aspect(1.)

    ax1.plot([-600,600],[600,-600],color='black')
    ax1.plot([0,0],[-2250,2250],color='black')


#    ax1.set_aspect(1.)
    divider = make_axes_locatable(ax1)
    axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax1)
    axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=ax1)

    plt.setp(axHistx.get_xticklabels()+axHisty.get_yticklabels(), visible=False)

    ybins = np.linspace(-2250,2250,75)
    xbins = np.linspace(-600,600,75)

    axHistx_n, bins, patches = axHistx.hist(vh, bins=xbins, color = 'black')
    axHisty_n, bins, patches = axHisty.hist(vpec, bins = ybins, orientation='horizontal', color='black')

    for tl in axHistx.get_xticklabels():
        tl.set_visible(False)

    
    if np.max(axHistx_n) > 100.0:
        yticks = [30,60,90,120]
    else:
        yticks = [20,40,60,80]

    axHistx.set_yticks(yticks)
#    axHistx.set_yticks([30,60,90,120])
    for tl in axHisty.get_yticklabels():
        tl.set_visible(False)
    
    if np.max(axHisty_n) > 100.0:
        yticks = [30,60,90,120]
    else:
        yticks = [20,40,60,80]

    axHisty.set_xticks(yticks)
    
    plt.tight_layout()
    ax1.set_ylim(-2250,2250)
#    plt.tight_layout()
    fig.savefig(abs_dir + "vH_vpec_1D-hist.png")
    plt.close(fig)
    ###############
    # color by HI column
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)

    
    mappable = ax1.scatter(vh,vpec,s=point_size,c=np.log10(N),cmap=cmap)

  #  x0,x1 = ax1.get_xlim()
  #  y0,y1 = ax1.get_ylim()
   # ax1.set_aspect((x1-x0)/(y1-y0))    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right",size="2.5%",pad=0.001)

   # cbar = fig.colorbar(mappable)
   # cbar.set_label(labels["HI_Column_Density"])
    
    ax1.set_ylim(-2250,2250)
    ax1.set_xlim(-600,600)
    
    
    ax1.plot([-600,600],[600,-600],color='black')
    ax1.plot([0,0],[-2250,2250],color='black')
    
    ax1.minorticks_on()
 #   plt.minorticks_on()
    plt.colorbar(mappable,label=labels["HI_Column_Density"],cax=cax)#,cax=cax)
    ax1.set_xlabel(labels[xname])
    ax1.set_ylabel(labels[yname])
    plt.tight_layout()
    ax1.set_ylim(-2250,2250)
    fig.savefig(abs_dir + xname +"_"+yname+"_Ncolor_dist.png")
    plt.close(fig)
    
    # printing some statistics for this plot
    total_number_absorbers = np.size(vh) * 1.0
    logN = np.log10(N)
    frac_all = 1.0 - np.size(vh[np.sign(vh+vpec) == np.sign(vh)]) / total_number_absorbers
    frac_12  = np.size(vh[(np.sign(vh+vpec) == np.sign(vh)) * (logN>=12.0)*(logN<13.0)])
    frac_12  = 1.0 - frac_12 / (1.0 * np.size(vh[(logN>=12.0)*(logN<13.0)]))
    frac_13  = np.size(vh[(np.sign(vh+vpec) == np.sign(vh)) * (logN>=13.0)*(logN<15.0)])
    frac_13  = 1.0 - frac_13 / (1.0 * np.size(vh[(logN>=13.0)*(logN<15.0)]))
    frac_15  = np.size(vh[(np.sign(vh+vpec) == np.sign(vh)) * (logN>=15.0)*(logN<17.0)])
    frac_15  = 1.0 - frac_15 / (1.0 * np.size(vh[(logN>=15.0)*(logN<17.0)]))
    frac_17  = np.size(vh[(np.sign(vh+vpec) == np.sign(vh)) * (logN>=17.0)            ])
    frac_17  = 1.0 - frac_17 / (1.0 * np.size(vh[(logN>=17.0)            ]))


    print "Statistics on vh vpec distribution - fraction of abs in ul and lr quadrants"
    print "all", frac_all
    print "12-13", frac_12
    print "13-15", frac_13
    print "15-17", frac_15
    print "17+", frac_17 
    
    ###loooooppppp#######
    Nmin, Nmax = 10.0,11.0
    while Nmax < 22.0:
    
        frac = np.size( vh[ (np.sign(vh+vpec) == np.sign(vh)) * (logN>=Nmin)*(logN<Nmax) ] )
        total = (1.0 * np.size(vh[(logN>=Nmin)*(logN<Nmax)]))
        if total > 0:
            frac = 1.0 - frac / total
        else:
           total = np.float('nan')
        
     
    
        print Nmin, Nmax, frac
        Nmin += 1.0
        Nmax += 1.0
    
        ###############
    # color by Metals
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)

    Z = get_params_as_array(alist,"Metallicity")
    mappable = ax1.scatter(vh,vpec,s=point_size,c=np.log10(Z),cmap=cmap)

    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["Metallicity"])
    ax1.set_xlabel(labels[xname])
    ax1.set_ylabel(labels[yname])
    plt.tight_layout()
    fig.savefig(abs_dir + xname +"_"+yname+"_Zcolor_dist.png")
    plt.close(fig)
    
            ###############
    # color by r
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)

    pos = get_params_as_array(alist,"pos")
    r = np.zeros(np.size(N))
    i = 0
    for p in pos:
        r[i] = ((p[0]**2 + p[1]**2 + p[2]**2)**0.5 * yt.units.cm).convert_to_units('Mpc').value / R_vir# pf['Mpc']/pf['cm'] /R_vir
        i = i + 1


    mappable = ax1.scatter(vh,vpec,s=point_size,c=r,cmap=cmap)
    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["r"])
    ax1.set_xlabel(labels[xname])
    ax1.set_ylabel(labels[yname])
    plt.tight_layout()
    fig.savefig(abs_dir + xname +"_"+yname+"_rcolor_dist.png")
    plt.close(fig)
    
def plot_vtot_N(alist,abs_dir):
    """
    """
    labels=labeldict()
    
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    vh = get_params_as_array(alist,"vh")
    N = get_params_as_array(alist,"HI_Column_Density")
    Z = get_params_as_array(alist,"Metallicity")
    vpec = get_params_as_array(alist,"vpec")
    # make sure vH is normalized...
    if np.min(vh) >= 0.0:
        mid = 0.5*(np.min(vh) + np.max(vh))
        
        vh = vh - mid
        
    
    
    mappable = ax1.scatter(vh+vpec,N,c=np.log10(Z),s=point_size,cmap=cmap)
    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["Metallicity"])
    ax1.semilogy()
    ax1.set_xlabel(labels['vh'] + " + " + labels['vpec'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    fig.savefig(abs_dir + "vtot_N_dist.png")
    plt.close(fig)
    
def plot_vpec_N(alist,abs_dir):
    """
    """
    labels=labeldict()
    
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    vpec = get_params_as_array(alist,"vpec")
    N = get_params_as_array(alist,"HI_Column_Density")
    #Z,vpec = get_params_as_array(alist,"Metallicity","vpec")
    # make sure vH is normalized...
    
    
    mappable = ax1.scatter(vpec,N,s=point_size,cmap=cmap)
   # cbar = fig.colorbar(mappable)
    #cbar.set_label(labels["Metallicity"])
    ax1.semilogy()
    ax1.set_xlabel(labels['vpec'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    fig.savefig(abs_dir + "vpec_N_dist.png")
    plt.close(fig)

def plot_dist_N(alist,abs_dir):
    """
    """
    
    labels=labeldict()

    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    pos = get_params_as_array(alist,"pos")
    N = get_params_as_array(alist,"HI_Column_Density")
    Z = get_params_as_array(alist,"Metallicity")
    T = get_params_as_array(alist,"Temperature")
    r,b = np.zeros(np.size(N)), np.zeros(np.size(N))
    i = 0

    if ORIENTATION == 'z':
        x_b = 0
        y_b = 1
    elif ORIENTATION == 'x':
        x_b = 1
        y_b = 2
    elif ORIENTATION == 'y':
        x_b = 0
        y_b = 2


    for p in pos:
        r[i] = ((p[0]**2 + p[1]**2 + p[2]**2)**0.5 * yt.units.cm).convert_to_units('Mpc').value / R_vir #/ pf['cm'] * pf['Mpc'] / R_vir
        b[i] = ((p[x_b]**2 + p[y_b]**2)**0.5 *yt.units.cm).convert_to_units('Mpc').value / R_vir #pf['cm'] * pf['Mpc'] / R_vir
        i = i + 1
    T= np.log10(T)
    Z = np.log10(Z)
    
    print 'correlation coefficient for projected and actual distance and column'
    print np.corrcoef(r,N)
    print np.corrcoef(b,N)

    xmin, xmax = 0.0,4.5
    ymin, ymax = 12.0,21.0

    mappable = ax1.scatter(r,np.log10(N),s=point_size,cmap=cmap)

    rbins = np.arange(0.0,4.9,0.2)
    rbins_2 = np.arange(0.0,5.1,0.2)
    rbin_cent = 0.5 * (rbins_2[1:] + rbins_2[0:-1])
    Navg  = np.zeros(np.size(rbins))

    i = 0
    print "average column in radial bins of absorbers"
    for rmin in rbins:
        rmax = rmin + (rbins[1] - rbins[0])
        select = (r >= rmin)*(r < rmax)
        Navg[i] = np.average(np.log10(N[select]))

        i = i + 1 
        print rmin, rmax, np.average(np.log10(N[select]))

#    ax1.semilogy()
    ax1.plot(rbin_cent, Navg, lw = linewidth, label = 'Average')
    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
#    xmin,xmax = ax1.get_xlim()
 #   ymin,ymax = ax1.get_ylim()
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)

    fig.savefig(abs_dir + "dist_N_dist.png")
    plt.close(fig)
    
    

#    ymin,ymax = np.log10(ymin), np.log10(ymax)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    logN = np.log10(N)
    xbinsize = 0.2
    ybinsize = 0.2
    xbins = np.arange(0.0,4.5+xbinsize,xbinsize)
    ybins = np.arange(12,21+ybinsize,ybinsize)
    
    rmesh, Nmesh, H2d = my_pt.my_histogram2d(r,logN,logscale=True)

    mappable = ax1.pcolormesh(rmesh, Nmesh, H2d, cmap=plt.cm.gist_heat)
    plt.cm.gist_heat.set_bad('black',1.)

    divider = make_axes_locatable(ax1)
    cax     = divider.append_axes('bottom',size="2.5%",pad=0.6)
    plt.colorbar(mappable, label = 'log Counts', cax=cax, orientation='horizontal')

    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    fig.savefig(abs_dir + "dist_N_dist_heat.png")
    plt.close(fig)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    logN = np.log10(N)
    xbins = np.arange(xmin,xmax+0.5,0.5)
    ybins = np.arange(ymin,ymax+0.5,0.5)

    bmesh, Nmesh, H2d = my_pt.my_histogram2d(b,logN,logscale=True)

    mappable = ax1.pcolormesh(bmesh, Nmesh, H2d, cmap=plt.cm.gist_heat)
    plt.cm.gist_heat.set_bad('w',1.)

    divider = make_axes_locatable(ax1)
    cax     = divider.append_axes('bottom',size="2.5%",pad=0.6)
    plt.colorbar(mappable, label = 'log Counts', cax=cax, orientation='horizontal')

    ax1.set_xlabel('Projected' + labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    fig.savefig(abs_dir + "dist_N_projected_dist_heat.png")
    plt.close(fig)





#    mappable = ax1.scatter(r,
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    ax1.scatter(b,np.log10(N),s=point_size)
#    ax1.semilogy()
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    ax1.set_xlabel('Projected ' + labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    
    print "average column in radial bins (projected) of absorbers"
    i = 0
    Navg = np.zeros(np.size(rbins))
    for rmin in rbins:
        rmax = rmin + (rbins[1] - rbins[0])
        select = (b >= rmin)*(b < rmax)

        Navg[i] = np.average(np.log10(N[select]))
        print rmin, rmax, np.average(np.log10(N[select]))
        i = i + 1

    



    ax1.plot(rbin_cent, Navg, lw = linewidth, label = 'Average')
    fig.savefig(abs_dir + 'dist_N_projected_dist.png')
    plt.close(fig)


    #######
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)    
    mappable = ax1.scatter(r,N,c=T,s=point_size,cmap=cmap)
    ax1.semilogy()
    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    
    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["Temperature"])
    fig.savefig(abs_dir + "dist_N_T_dist.png")
    plt.close(fig)
    ####
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)    
    mappable = ax1.scatter(r,N,c=Z,s=point_size,cmap=cmap)
    ax1.semilogy()
    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(labels['HI_Column_Density'])
    
    cbar = fig.colorbar(mappable)
    cbar.set_label(labels["Metallicity"])
    fig.savefig(abs_dir + "dist_N_Z_dist.png")
    plt.close(fig)

def plot_dist_hist(alist,abs_dir,conv):

    labels=labeldict()

    fig = plt.figure(figsize = [10,8])
    ax1 = fig.add_subplot(111)
    
    pos = get_params_as_array(alist,"pos")
    N   = get_params_as_array(alist,"HI_Column_Density")
    r = np.zeros(np.size(N))
    b = np.zeros(np.size(N))
    i = 0

    if ORIENTATION == 'z':
        x_b = 0
        y_b = 1
    elif ORIENTATION == 'x':
        x_b = 1
        y_b = 2
    elif ORIENTATION == 'y':
        x_b = 0
        y_b = 2


    for p in pos:
        r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5 
        b[i] = (p[x_b]**2 + p[y_b]**2)**0.5
        i = i + 1
        
    r = (r * yt.units.cm).convert_to_units('Mpc').value/R_vir# pf['Mpc'] / pf['cm'] / R_vir
    b = (b * yt.units.cm).convert_to_units('Mpc').value/R_vir# pf['Mpc'] / pf['cm'] / R_vir
        
    NminCuts = np.array([1.0  , 1.0,    10.0**13, 10.0**15, 10.0**17])
    NmaxCuts = np.array([10.0**25, 10.0**13, 10.0**15, 10.0**17, 10.0**25]   )
        
    
    bins = np.arange(0.0,4.6,0.5)

    for Nmin,Nmax in zip(NminCuts,NmaxCuts):
        print Nmin,Nmax
        rcut = r[ (N > Nmin) * (N < Nmax) ]
        

        hist, bins = np.histogram(rcut, bins=bins, density=False)

        center = (bins[:-1] + bins[1:])*0.5
    

        vol = 4.0/3.0 *np.pi*((conv*bins[1:])**3 - (conv*bins[:-1])**3)
        hist = hist / vol
        logNmin = np.log10(Nmin)
        logNmax = np.log10(Nmax)
        ax1.plot(center, hist, label=r"%2d $<$ log(N$_{HI}$) $<$ %2d cm$^{-2}$"%(logNmin,logNmax))
        ax1.scatter(center,hist, color='black')
        
        if Nmin == 10**13 and Nmax == 10**15:
            print (hist*vol)[0], (hist*vol)[1]

    ax1.set_xlim(0,4.5)
    ax1.set_ylim(0.0,1.5)

    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(r'n$_{abs}$ Mpc$^{-3}$ ')
   
    ## Do this again for the projected distance@  
    #ax1.set_ylabel(r'Counts')
    ax1.legend(loc='best')
    fig.savefig(abs_dir + "dist_vol_histogram.png")
    plt.close(fig)
    fig = plt.figure(figsize = [10,8])
    ax1 = fig.add_subplot(111)  
    bins = np.arange(0.0,3.6,0.5)
    for Nmin,Nmax in zip(NminCuts,NmaxCuts):

        bcut = b[ (N > Nmin) * (N < Nmax) ]
        

        hist, bins = np.histogram(bcut, bins=bins, density=False)
        width = 0.85*(bins[1:]-bins[:-1])
        center = (bins[:-1] + bins[1:])*0.5
    

        area = np.pi*((conv*bins[1:])**2 - (conv*bins[:-1])**2)
        hist = hist / area
        logNmin = np.log10(Nmin)
        logNmax = np.log10(Nmax)
        ax1.plot(center, hist, label=r"%2d $<$ log(N$_{HI}$) $<$ %2d cm$^{-2}$"%(logNmin,logNmax))
        ax1.scatter(center,hist, color='black')
        
        if Nmin == 10**13 and Nmax == 10**15:
            print (hist*area)[0], (hist*area)[1]
        
    ax1.set_ylim(0.0,12.0)
    ax1.set_xlim(0,3.5)
    #ax1.set_ylim(0.0,1.0)
    ax1.axvspan(2.5,3.5,color='g',alpha=0.5)
    ax1.set_xlabel(r'$\rho$ (R$_{vir}$)')
    ax1.set_ylabel(r'n$_{abs}$ Mpc$^{-2}$ ')
   
    #ax1.set_ylabel(r'Counts')
    ax1.legend(loc='best')
    fig.savefig(abs_dir + "impact_area_histogram.png")
    plt.close(fig)
    
    
    
    
    
    
    
    
def plot_temp_hist(alist,abs_dir):

    labels=labeldict()

    fig = plt.figure(figsize = [8,6])
    ax1 = fig.add_subplot(111)
    
    pos = get_params_as_array(alist,"pos")
    T = get_params_as_array(alist,"Temperature")
   # r = np.zeros(np.size(N))
  #  i = 0
  ##  for p in pos:
  # #     r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
   #     i = i + 1
 #   bins = np.arange(0.25,3.5,0.25)
    bins = np.logspace(np.min(np.log10(T)),np.max(np.log10(T)),15)
   # print bins
    hist, bins = np.histogram(T, bins=bins)#, density=True)
    width = 0.85*(bins[1:]-bins[:-1])
    center = (bins[:-1] + bins[1:])*0.5
    #print np.sum(hist)
    total = np.sum(hist)*1.0
    ax1.bar(center, hist/total, align='center',width=width,facecolor='black',alpha=0.75)
#plt.semilogy()
    ##  ax1.set_ylim(0,1)
    ax1.semilogx()
    ax1.set_xlabel(labels['Temperature'])
    ax1.set_ylabel(r'Fraction')
    fig.savefig(abs_dir + "temperature_histogram.png")
    plt.close(fig)
    
    
def plot_metallicity_hist(alist,abs_dir):

    labels=labeldict()

    fig = plt.figure(figsize = [8,6])
    ax1 = fig.add_subplot(111)
    
    pos = get_params_as_array(alist,"pos")
    Z = get_params_as_array(alist,"Metallicity")
   # r = np.zeros(np.size(N))
  #  i = 0
  ##  for p in pos:
  # #     r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5
   #     i = i + 1
 #   bins = np.arange(0.25,3.5,0.25)
    bins = np.logspace(np.min(np.log10(Z)),np.max(np.log10(Z)),10)
   # print bins
    hist, bins = np.histogram(Z, bins=bins)#, density=True)
    width = 0.85*(bins[1:]-bins[:-1])
    center = (bins[:-1] + bins[1:])*0.5
    #print np.sum(hist)
    total = np.sum(hist)*1.0
    ax1.bar(center, hist/total, align='center',width=width,facecolor='black',alpha=0.75)
#plt.semilogy()
    ##  ax1.set_ylim(0,1)
    ax1.semilogx()
    ax1.set_xlabel(labels['Metallicity'])
    ax1.set_ylabel(r'Fraction')
    fig.savefig(abs_dir + "metallicity_histogram.png")
    plt.close(fig)   
    
    
    
def plot_dist_hist_cut(alist,abs_dir,conv):

    labels=labeldict()

    fig = plt.figure(figsize = [8,6])
    ax1 = fig.add_subplot(111)
    
    pos = get_params_as_array(alist,"pos")
    N = get_params_as_array(alist,"HI_Column_Density")
    Z = get_params_as_array(alist,"Metallicity")

    r = np.zeros(np.size(N))
    i = 0
    for p in pos:
        r[i] = (p[0]**2 + p[1]**2 + p[2]**2)**0.5    
        i = i + 1
        
    r =  (r * yt.units.cm).convert_to_units('Mpc').value/R_vir#pf['Mpc'] / pf['cm'] / R_vir
    #r = r[N >10**16]    
    r = r[Z < 0.01]   
    bins = np.arange(0.0,3.5,0.25)
    hist, bins = np.histogram(r, bins=bins, density=False)
    width = 0.85*(bins[1:]-bins[:-1])
    center = (bins[:-1] + bins[1:])*0.5
    

    vol = 4.0/3.0 *np.pi*((conv*bins[1:])**3 - (conv*bins[:-1])**3)
    hist = hist / vol
    ax1.bar(center, hist, align='center',width=width,facecolor='black',alpha=0.75)
#plt.semilogy()
    ax1.set_xlim(0,3.5)
#    ax1.set_ylim(0,1)
    ax1.set_xlabel(labels['r'])
    ax1.set_ylabel(r'n$_{abs}$ Mpc$^{-3}$ ')
    fig.savefig(abs_dir + "dist_histogram_cut.png")
    plt.close(fig)
    
    
def plot_dist_length(alist,abs_dir):

    labels=labeldict()

    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)
    
    n  = get_params_as_array(alist,"HI_NumberDensity")
    N     = get_params_as_array(alist,"HI_Column_Density")
    #n = n / 1.67E-27
    l = N / n
    l = (l *yt.units.cm).convert_to_units('Mpc').value # pf['Mpc'] / pf['cm']
    print np.min(N), np.average(N), np.max(N)
    print np.min(n), np.average(n), np.max(n)
    
    print np.min(l), np.average(l), np.max(l)
    
    ax1.scatter(N,l,s = point_size)
    ax1.semilogx()
    ax1.semilogy()
    ax1.set_ylim(np.min(l)*0.5,np.max(l)*1.5)
    ax1.set_xlabel(labels['HI_Column_Density'])
    ax1.set_ylabel(r'Size (Mpc)')
    plt.tight_layout()
    
    fig.savefig(abs_dir + "l_N.png")
    plt.close(fig)
    

def plot_vel_dist(alist,abs_dir):
    plt.close()
    
    binsize = 100.0
    color = {'All gas': 'black','Cold': 'blue', 'Warm': 'yellow', 'Hot': 'red',
             'All_particle': 'red','DM': 'black', 'Star': 'green'}

    symbolDict = {'All gas' : '-o', 'Cold': '-v', 'Warm': '-s', 'Hot': '-p',
                 'All_particle': '-o', 'DM': '-s', 'Star': '-*'}
                 
    labels = labeldict()
    
    T = get_params_as_array(alist, "Temperature")
    v3 = get_params_as_array(alist, "Velocity")
    vpec = get_params_as_array(alist,'vpec')
    losArray = get_params_as_array(alist,'los')
    
    vmag = np.zeros(np.size(T))
    vlos = np.zeros(np.size(T))
    
    i = 0
    for v in v3:
        vmag[i] = (v[0]**2 + v[1]**2 + v[2]**2)**0.5
        i = i +1 
    i = 0
    for (v,los) in zip(v3,losArray):
        vlos[i] = (v[0]*los[0] + v[1]*los[1] + v[2]*los[2])/(np.sum(los*los)**0.5)
        i = i + 1 

    #print "vmag"
    #print np.min(vmag), np.max(vmag), np.size(vmag)

    vlos = (vlos * yt.units.cm).convert_to_units('km').value# pf['km'] / pf['cm']
    vmag = (vmag * yt.units.cm).convert_to_units('km').value #pf['km'] / pf['cm'] 

    #print "vmag"
    #print np.min(vmag), np.max(vmag), np.size(vmag)
    ##print "vlos"
    #print np.min(vlos), np.max(vlos), np.size(vlos)
    #print "vpec"
    #print np.min(vpec), np.max(vpec), np.size(vpec)


    bins = np.arange(0.0, 3001.0, binsize)

    TminDict = {"All gas": -np.inf,
                "Cold": -np.inf,
                "Warm": 2.0 * 10.0**4,
                "Hot": 10.0**6
               }


    TmaxDict = {"All gas": np.inf,
                "Cold": 2.0*10.0**4,
                "Warm": 10.0**6,
                "Hot": np.inf
               }

#    types = ['All gas',"Cold","Warm","Hot"]
    types = ['All gas', 'Cold', 'Warm']
    for gas_type in types:
        Tmin, Tmax = TminDict[gas_type], TmaxDict[gas_type]
    
    
        vmag_reg = vmag[ (T > Tmin) * (T < Tmax) ]
    
        hist, bins = np.histogram(vmag_reg, bins = bins)
        center     = (bins[:-1] + bins[1:])*0.5
        
        tot = 1.0 * np.sum(hist)
        
        vdisp = np.std(vmag_reg)

        if tot == 0:
            vdisp = 0

        #print center, hist
        #print np.shape(center), np.shape(hist/tot)
        
        plt.plot(center, hist/tot, symbolDict[gas_type],
                 label = r" " + gas_type + " %3d - %3d"%(vdisp,tot),
                 color = color[gas_type])
                 

    plt.xlim(0.0,3000.0)
    plt.xlabel(r'$\mid$v$\mid$ (km/s)')
    plt.ylabel(r'Fraction - Absorbers')
    plt.legend(loc='best',fancybox=True)
    
    plt.savefig(abs_dir + "3D_velocity_distribution.png")
    plt.close()
    
    bins = np.arange(-2000.0,2000.0,binsize)
    for gas_type in types:
        Tmin, Tmax = TminDict[gas_type], TmaxDict[gas_type]
    
    
        vlos_reg = vlos[ (T > Tmin) * (T < Tmax) ]
    
        hist, bins = np.histogram(vlos_reg, bins = bins)
        center     = (bins[:-1] + bins[1:])*0.5
        
        tot = 1.0 * np.sum(hist)
        
        vdisp = np.std(vlos_reg)
        if tot == 0:
            vdisp = 0
        #print center, hist
        #print np.shape(center), np.shape(hist/tot)
        
        plt.plot(center, hist/tot, symbolDict[gas_type],
                 label = r" " + gas_type + " %3d - %3d"%(vdisp,tot),
                 color = color[gas_type])
                 
    plt.xlim(-2000.0,2000.0)
  #  plt.xlim(bins[np.where(bins[hist>0])[0]-1]-binsize,bins[np.where(bins[hist>0])[0]+1]+binsize)
    plt.xlabel(r'v$_{\rm{los}}$ (km/s)')
    plt.ylabel(r'Fraction - Absorbers')
    plt.legend(loc='best',fancybox=True)
    
    plt.savefig(abs_dir + "vlos_distribution.png")
    plt.close()
    
    
    bins = np.arange(0.0,2000.0,binsize)
    for gas_type in types:
        Tmin, Tmax = TminDict[gas_type], TmaxDict[gas_type]
    
    
        vlos_reg = np.abs(vlos[ (T > Tmin) * (T < Tmax) ])
    
        hist, bins = np.histogram(vlos_reg, bins = bins)
        center     = (bins[:-1] + bins[1:])*0.5
        
        tot = 1.0 * np.sum(hist)
        
        vdisp = np.std(vlos_reg)
        if tot == 0:
            vdisp = 0
        #print center, hist
        #print np.shape(center), np.shape(hist/tot)
        
        plt.step(center, hist/tot, symbolDict[gas_type],
                 label = r" " + gas_type + " %3d - %3d"%(vdisp,tot),
                 color = color[gas_type])
                 
    plt.xlim(0.0,2000.0)
  #  plt.xlim(bins[np.where(bins[hist>0])[0]-1]-binsize,bins[np.where(bins[hist>0])[0]+1]+binsize)
    plt.xlabel(r'$\mid$v$_{\rm{los}}\mid$ (km/s)')
    plt.ylabel(r'Fraction - Absorbers')
    plt.legend(loc='best',fancybox=True)
    
    plt.savefig(abs_dir + "abs_vellos_distribution.png")
    plt.close()

    bins = np.arange(-2000.0,2000.0,binsize)
    for gas_type in types:
        Tmin, Tmax = TminDict[gas_type], TmaxDict[gas_type]
    
    
        vpec_reg = vpec[ (T > Tmin) * (T < Tmax) ]
    
        hist, bins = np.histogram(vpec_reg, bins = bins)
        center     = (bins[:-1] + bins[1:])*0.5
        
        tot = 1.0 * np.sum(hist)
        
        vdisp = np.std(vpec_reg)
        if tot == 0:
            vdisp = 0
        #print center, hist
        #print np.shape(center), np.shape(hist/tot)
        
        plt.plot(center, hist/tot, symbolDict[gas_type],
                 label = r" " + gas_type + " %3d - %3d"%(vdisp,tot),
                 color = color[gas_type])
                 
    plt.xlim(-2000.0,2000.0)
    #plt.xlim(np.min(vpec)-binsize,np.max(vpec)+binsize)
    #plt.xlim(bins[np.where(bins[hist>0])[0]]-binsize,bins[np.where(bins[hist>0])[0]]+binsize)
    plt.xlabel(r'v$_{\rm{pec}}$ (km/s)')
    plt.ylabel(r'Fraction - Absorbers')
    plt.legend(loc='best',fancybox=True)
    
    plt.savefig(abs_dir + "velpec_distribution-test.png")
    plt.close()
    
    
def plot_b_T_errors(alist,abs_dir,cfield_name='error'):

    qso_id = get_params_as_array(alist, 'QSO_id')
    T = get_params_as_array(alist, 'Temperature')
    b = get_params_as_array(alist, 'b')  * 1000.0 * 100.0 # to get to cm/s
    
    
    T_b = 0.5* b * b * CONSTANT_m_H / CONSTANT_k 
    
    print np.average(T_b)
    
    error = (T-T_b)/T
    
    T   = np.log10(T)
    T_b = np.log10(T_b)
    

    
    fig = plt.figure(figsize = [6,6])
    ax1 = fig.add_subplot(111)


    if cfield_name == 'error':
        colorField = error
    elif cfield_name == 'HI_Column_Cell':
        cellN = get_params_as_array(alist, 'HI_NumberDensity')
        cellL = get_params_as_array(alist, 'CellVolume')**(1.0/3.0)
        cell_HIN = np.log10(cellN * cellL)
        
        colorField = cell_HIN
        
    
    mappable = ax1.scatter(T,T_b,s=25,c=colorField,cmap='algae')

    ax1.set_xlabel(r'b$^{2}$m$_{\rm{H}}$/(2k$_{\rm{B}}$) (K)')
    ax1.set_ylabel(r'T (k)')

    x0,x1 = ax1.get_xlim()
    y0,y1 = ax1.get_ylim()
    
    x = np.min([x0,y0])
    y = np.max([x1,y1])
    print x, y
    ax1.plot([x,y],[x,y],linewidth=2.0,color='black')
    
    ax1.set_xlim([x,y])
    ax1.set_ylim([x,y])


    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right",size="2.5%",pad=0.001)
    ax1.minorticks_on()
    plt.colorbar(mappable,label='Fractional Error',cax=cax)





    plt.tight_layout()
    fig.savefig(abs_dir + "b_T_error.png")
    plt.close(fig)
    
    spec_HIN = np.log10(    get_params_as_array(alist,'HI_Column_Density'))
    error = np.log10(np.abs(error))
    plt.scatter(cell_HIN, error, s =15, c=np.log10(T),cmap='algae')
    plt.xlabel('log cell HIN')
    plt.ylabel('log Fractional Error')
#    plt.ylim(100,-)
    plt.savefig(abs_dir + "cell_HIN_bt-error.png")
    plt.close()
    

    plt.scatter(spec_HIN, error, s =15, c=np.log10(T),cmap='algae')
    plt.xlabel('log spec HIN')
    plt.ylabel('log Fractional Error')
#    plt.ylim(100,-)
    plt.savefig(abs_dir + "spec_HIN_bt-error.png")
    plt.close()
    
    plt.scatter(cell_HIN,spec_HIN,
                s = 15, c = error, cmap='RdBu')
    plt.plot([12,22],[12,22],color='black')
    plt.xlabel('log cell HIN')
    plt.ylabel('log spec HIN')
    plt.savefig(abs_dir + "cellHIN_specHIN_bterror.png")
    plt.close()
    
    print '---- id ----'
    print qso_id[0][ (cell_HIN < 12.0)*(spec_HIN > 18.0) ]
    print spec_HIN[ (cell_HIN < 12.0)*(spec_HIN > 18.0) ]
    
 
def plot_radial_vel(alist,abs_dir):
    v3 = get_params_as_array(alist, "Velocity")
    pos = get_params_as_array(alist,"pos")
    N = get_params_as_array(alist,"HI_Column_Density")
    N = np.log10(N)
    print v3
    print pos
    i = 0
    v_radial = np.zeros(np.size(pos)/3.0)
    r = np.zeros(np.size(pos)/3.0)  
#  print np.size(pos)
    for p,v in zip(pos,v3):
        p = (-1.0*p * yt.units.cm).convert_to_units('km').value# pf['km'] / pf['cm']  # vector from absorber to cluster center (0,0,0)
        v = (v * yt.units.cm).convert_to_units('km').value #pf['km'] / pf['cm']
        r[i] = (np.sum(p*p)**0.5 * yt.units.cm).convert_to_units('Mpc').value / R_vir #pf['Mpc'] / pf['km'] / R_vir     

#   print p, v
        v_radial[i] = np.sum(p*v) / (1.0*np.sum(p*p)**0.5)
        i = i + 1
 
    
    print v_radial
    vmin,vmax = np.min(v_radial), np.max(v_radial)

    vmin,vmax = -500.0, 2500.0
    print vmin,vmax
    bins = np.linspace(vmin,vmax,50.0)
    
    
    
    
    hist, bins = np.histogram(v_radial,bins=bins)
    centers = 0.5*(bins[:-1] + bins[1:])
    fractional = False
    if fractional:
        plt.step(centers, hist/(1.0*np.sum(hist)), label = 'All', lw=linewidth)
    else:
        plt.step(centers, hist, label = 'All', lw=linewidth)
    
    Nlimits = [   [12.0,13.0], [13.0,14.0], [14.0,15.0], [16.0,25.0] ]
    labels  = ['12 - 13', '13 - 14', '14 - 15', '16 - inf']

    rlimits = [   [0.0,1.0], [1.0,2.0], [2.0,3.0], [3.0,4.0], [4.0,6.0] ]
    labels = ['0 - 1', '1 - 2', '2 - 3', '3 - 4', '4 - 6' ]
#    colors  - ['green', 
    i = 0
#    for Nlim in Nlimits:
    for rlim in rlimits:   
        nselect = (r >= rlim[0])*(r<rlim[1])
        vselect = v_radial[nselect]
        
        
        hist, bins = np.histogram(vselect, bins=bins)
        if fractional:
            plt.step(centers, hist/(1.0*np.sum(hist)), label = labels[i], lw = linewidth)
        else:
            plt.step(centers,hist, label = labels[i], lw = linewidth)
        i = i +1
        
    
    print "12-13 with 13-14,14-15,15 and up"
    print ks_2samp( v_radial[(N >= 12)*(N<13)] , v_radial[(N >= 13)*(N<14)]  )
    print ks_2samp( v_radial[(N >= 12)*(N<13)] , v_radial[(N >= 14)*(N<15)]  )
    print ks_2samp( v_radial[(N >= 12)*(N<13)] , v_radial[(N >= 15)*(N<25)]  )
    print "13-14 with 14-15, 15 and up"
    print ks_2samp( v_radial[(N >= 13)*(N<14)] , v_radial[(N >= 14)*(N<15)]  )
    print ks_2samp( v_radial[(N >= 13)*(N<14)] , v_radial[(N >= 15)*(N<25)]  )
    print "14-15 with 15 and up"
    print ks_2samp( v_radial[(N >= 14)*(N<15)] , v_radial[(N >= 15)*(N<25)]  )
            
    print "ENDING KS TESTING"
    plt.xlabel('Velocity Towards Cluster Center (km/s)')

    if fractional:
        plt.ylabel('Fraction')
    else:
        plt.ylabel('Counts')

    plt.legend()
    plt.minorticks_on()
    plt.xlim(vmin,vmax)
    
    plt.savefig(abs_dir + 'radial_velocity_rbins.png')
    plt.ylim(0.0,10.0)
    plt.savefig(abs_dir + 'radial_velocity_rbins_zoom.png')
    plt.close()




plot_radial_vel(alist,abs_dir)    
#    pos, Z = get_params_as_array(alist,"pos","Metallicity") 
plot_b_T_errors(alist,abs_dir,cfield_name = 'HI_Column_Cell')
plot_dist_length(alist,abs_dir)  
# make plots
plot_N_Z(alist,abs_dir)
plot_vH_N(alist,abs_dir)
plot_vH_vpec(alist,abs_dir)
plot_vtot_N(alist,abs_dir)
plot_dist_N(alist,abs_dir)
#plot_dist_hist(alist,abs_dir,R_vir)
#plot_dist_hist_cut(alist,abs_dir,R_vir)
plot_temp_hist(alist,abs_dir)
plot_metallicity_hist(alist,abs_dir)
plot_vel_dist(alist,abs_dir)
