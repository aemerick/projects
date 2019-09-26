### The purpose of this script is to 

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import plotting.scatter3D as s3d # for arrow drawing
import yt

# for saving objects
try:
    import cPickle as pickle
except:
    import pickle

## mencoder "mf://*.png" -mf w=800:h=600:fps=2:type=png -ovc copy -oac copy -o movie.avi


DEFAULT_ARROW_TIME = 5.0E16

class GalCluster:
    def __init__(self, ds, center, R_vir, inputfile=None):
        self.ds = ds
        self.center= center
        self.R_vir = R_vir *yt.units.Mpc# in Mpc!!
        
        if not inputfile == None:
            try:
                f = file(inputfile,'rb')
                c = pickle.load(f)
                self = c
                f.close()
            except:
                print inputfile , " does not exist"        
        
        
    def save(self, output):
        f = file(output,'wb')
        pickle.dump(self, f, -1)
        f.close()
    
    


def saveAbsorberList(absorbers, outputfile):
    """
        Saves the list of absorbers with pickle
    """
    
    #ta daaa
    f = file(outputfile,'wb')
    pickle.dump(absorbers, f, -1)
    f.close()
    
    print "Saved as pickle: ",outputfile
    
def loadAbsorberList(inputfile):
    """
        Loads the list of absorbers with pickle
    """
    
    f = file(inputfile,'rb')
    absorbers = pickle.load(f)
    f.close()
    
    return absorbers


class Absorber:
    """ making a class for the absorbers to store data better for
        post processing after initial render """
    absorberCount = 0
        
    
    def __init__(self, pos, param, plot = True):
        self.pos = pos          # 3d array giving position
        self.param = param
   
   #     self.vel = vel          # 3d array giving velocity
   #     self.N   = N            # Fitted HI column
        self.plot = plot        # bool ... should this be plotted?
        Absorber.absorberCount +=1
        
    # Functions to set attributes if needed
    def setPos(self, pos):
        self.pos = pos
        
  #  def setVel(self, vel):
  #      self.vel = vel
    
  #  def setN(self, N):
  #      self.N = N
        
    def setPlot(self, plot):
        self.plot = plot
        
    def addParams(self, new_params):
        """
        Add a dictionary of new parameters
        """
    
        for c in new_params:
            self.param[c] = new_params[c]
      
    def list_available_parameters(self):
        for c in self.param:
            print c  

    # Calculates arrow end point for velocity    
    def arrowEndPoint(self, t = DEFAULT_ARROW_TIME ):
        """
        Calculates arrow end point
        """
        
        return t*self.param["Velocity"] + self.pos
    
    
    
    
        
class AbsorberPlot:
    """
    Class for the absorber plot. Initialize with a list of absorbers to plot
    """

    def __init__(self, absorbers, cluster, cparam="HI_Column_Density", 
                 cscale="log"):
        """
        Must have a list of absorber class objects and a galaxy cluster
        class object to initialize
        """
        self.absorbers = absorbers
        self.cluster = cluster
        
        self.cscale = cscale
        self.cparam = cparam
        
    

    # assign figure and ax to pre-existing
    def set_plot(self, fig, ax, points=None):
        """
        Pass a figure, axis, and potentially the points mapper to set the 
        plot
        """
        self.fig = fig
        self.ax  = ax
        self.points = points

    # generate a new figure and ax
    def create_plot_base(self):
        """
        Generate plt.figure and 3D axes object for plotting purposes
        """
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d')
        
        self.fig = fig
        self.ax  = ax
#        mappable = ax.scatter
    def clear_plot(self):
        del self.cbar
        del self.ax
        del self.points
        del self.fig

    def set_N_cut(self, Nmin, Nmax):
        self.Nmin = Nmin
        self.Nmax = Nmax

    def set_velocity_cut(self, vmin, vmax, absvel=False):
        a = self.absorbers[0]
        pf = self.cluster.ds

        conversion = 1000.0*100.0*pf['Mpc']/pf['cm'] / self.cluster.R_vir

        vmin,vmax = vmin*conversion,vmax*conversion

        self.vmin = vmin
        self.vmax = vmax
        self.absvel = absvel

    def set_T_cut(self, Tmin, Tmax):
        """ 
        Temperature cutoffs
        """
        self.Tmin = Tmin
        self.Tmax = Tmax


    def make_cuts(self):
        """
        Function makes cuts with vmin, vmax, Nmin, Nmax attributes... if
        they do not exist, print a statement saying so as a warning to the
        user
        """
        try:
            vmin, vmax = self.vmin, self.vmax
            
            # right now just looking at velocity projected along its position
            # vector (i.e. towards or away from cluster center)
            if not self.absvel:
                for a in self.absorbers:
                    v = a.param['Velocity']
                    pos = a.pos
                    vproj = v * pos / (pos[0]**2 + pos[1]**2 +pos[2]**2)**0.5
                    vproj = np.sum(vproj)
                    
                    if vproj < vmin or vproj > vmax:
                        a.plot = False
                    else:
                        a.plot = True
            
            elif self.absvel and (vmin<0 or vmax<0):
                # throw exception... vmin vmax must be non-zero
                raise ValueError("vmin vmax sign error with absvel")
            else:
            
                for a in self.absorbers:
                    v = a.param['Velocity']
                    vmag = (v[0]**2 + v[1]**2 + v[2]**2)**0.5                       
                    
                    if vmag > vmax or vmag < vmin:
                        a.plot = False
                    else:
                        a.plot = True

                # make the cutoff in v_r
       
        except ValueError:
            print "vmin and vmax must be non-zero when absvel is True"
        
        except AttributeError:
            print "Using all velocities"
#            print "use set_velocity_cut to set cuts"
            
            
            
        # try to make the HI Column cuts
        try:
            Nmin, Nmax = self.Nmin, self.Nmax
            for a in self.absorbers:
                a.param['HI_Column_Density']
                if N > Nmax or N < Nmin:
                    a.plot = False
                else:
                    a.plot = True
                    
                    
            print "Made HI column cuts. Including points for:"
            print "%3.2e < N < %3.2e" % (Nmin, Nmax)
        except AttributeError:
            print "Using all column densities"
#            print "use set_N_cut to set cuts"

        try:
            Tmin, Tmax = self.Tmin, self.Tmax
            for a in self.absorbers:
                T = a.param['Temperature']
                if T > Tmax or T < Tmin:
                    a.plot = False
                else:
                    a.plot = True
            print "Made temperature cuts. Including points for:"
            print "%3.2e < T < %3.2e" % (Tmin, Tmax)
        except AttributeError:
            print "Using all temperatures"


        try:
            Zmin, Zmax = self.Zmin, Zmax

            for a in self.absorbers:
                Z = a.param['Metallicity']
                if Z > Zmax or Z < Zmin:
                    a.plot = False
                else:
                    a.plot = True
            print "Made Metallicity cuts. Including points for:"
            print "%3.2e < Z < %3.2e" %(Zmin,Zmax)
        except AttributeError:
            print "Using all Metallicities"


        self.plot_points()
                
    def remove_cuts(self):
        """
            Removes cuts made... aka sets plot to True for all points
        """
        for a in self.absorbers:
            a.plot = True


    def plot_points(self, cparam=None):
        """
        Nmin, Nmax, vmin, vmax are the HI column and velocity
        cutoffs. Velocity cutoff is vel magnitude AND is in km/s.
        absvel=True set the velocity cutoffs as absolute velocity
        """
        self.cparam = cparam 
        cmin,cmax = None,None
        if self.cparam == "HI_Column_Density":
            cmin,cmax = 1.0E9,1.0E19
        elif self.cparam == "Temperature":
            cmin,cmax = 1.0E3, 1.0E6
            
        xpos = []
        ypos = []
        zpos = []
        all_color   = []
        
        if not cparam == None:
            self.cparam = cparam
        
        try:
            self.points.remove()
            print "Cleared previous ponits"
        except AttributeError:
            print "Making new plot"

        for a in self.absorbers:
 
            if a.plot:
                pos = a.pos
                
                
                #conv = self.cluster.pf['Mpc'] / self.cluster.pf['cm'] / self.cluster.R_vir
                
                
                
                
                xpos.append((pos[0]*yt.units.cm).convert_to_units('Mpc').value / self.cluster.R_vir)
                ypos.append((pos[1]*yt.units.cm).convert_to_units('Mpc').value / self.cluster.R_vir )
                zpos.append((pos[2]*yt.units.cm).convert_to_units('Mpc').value / self.cluster.R_vir)

                # this is bad... need to fix to generalize to take an arbitrary
                # color value decided by user at any point.
                
                all_color.append(a.param[cparam])
                
        if self.cscale == "log":
            all_color = np.array(all_color)
            all_color = np.log10(all_color)
            if not cmin ==None and not cmax == None:
                cmin, cmax = np.log10(cmin), np.log10(cmax)
                                            
        # draw the points... make the points mappable
        points = self.ax.scatter(xpos,ypos,zpos, c = all_color,
                            s=50, cmap='spectral',
                            vmin=cmin,vmax=cmax)

        # set the points object
        self.points = points
        
        # color bar
       # if np.size(self.fig.axes) > 1:
       #     self.fig.delaxes(self.fig.axes[1])
     
       
        try:
            cbar = self.cbar
           
            if self.cscale == "log":
                cbar.set_label(self.cscale + self.cparam)
            else:
                cbar.set_label(self.cparam)
                
            cbar.set_clim(cmin,cmax)
            cbar.update_bruteforce(self.points)            


        except AttributeError:
       
            cbar = self.fig.colorbar(self.points)
            if self.cscale == "log":
                cbar.set_label(self.cscale + self.cparam)
            else:
                cbar.set_label(self.cparam)
      #  self.cbar = cbar
        # redraw
        
        self.cbar = cbar
        self.fig.canvas.draw()
        plt.draw()


    # draw a sphere of radius r in plot units and color        
    def drawSphere(self, r, color):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x=np.cos(u)*np.sin(v) 
        y=np.sin(u)*np.sin(v)
        z=np.cos(v)
        
        self.ax.plot_wireframe(r*x,r*y,r*z, color=color)
    
    def drawArrow(self, pos0, pos1, *args, **kwargs):
        """
        Takes 2 3D positions and draws the arrow. just a wrapper on
        the Arrow3D method. Takes in kwargs.        
        
        """
        
        s3d.drawArrow(np.array([pos0,pos1]).flatten(), self.ax, *args, **kwargs)
        
        
    def drawVelocityArrows(self, t=DEFAULT_ARROW_TIME):
        """
        Draw velocity arrows for all absorbers
        """
        # need to draw arrows all at once.... not one by one
        
        for a in self.absorbers:
            if a.plot:
                xi = a.pos
                xf = a.arrowEndPoint(t)
        
                coords = np.array([xi,xf]).flatten()
                s3d.drawArrow(coords, self.ax, color = "black")
            
    def clearArrows(self):
        """
        Need a function to remove the drawn arrows
        """
        
        print "THIS DOES NOTHING RIGHT NOW!!!!"
        
        
    def clearPoints(self):
        """
        Need a function to remove the drawn points
        """
        
        print "THIS DOES NOTHING RIGHT NOW!!!!"
    
    def rotate_plot(self, elevation, azim, outputStr,counter_persist=False):
        """
        set counter_persist to true in order to keep the output image number 
        counter going. This allows for more complex movie making. For example,
        changing color scale after a rotation in the axis without having
        to do some messy editing of the image frame names after generation. Set
        False (default) to start from 0 each time.
        """
        
        
        # check if counter exists already
        if counter_persist:
            try:
                i = self.frame_counter
            except:
                i = 0
        else:
            i = 0
        
        for az in azim:
            self.ax.view_init(elev=elevation,azim=az)
            self.fig.canvas.draw()
            self.fig.savefig(outputStr + "%06d.png"%(i))
            i = i+1
        
        if counter_persist:
            self.frame_counter = i
        
   # def animate(self, 
    


def losVel(los,ray):
    """
    Calculate the scalar line of sight velocity given the line of sight
    vector and the yt ray object. Does this for all cells along the ray
    """
    ray_sort = np.argsort(ray['t'])
    
    vx = ray['x-velocity'][ray_sort]
    vy = ray['y-velocity'][ray_sort]
    vz = ray['z-velocity'][ray_sort]
    
    vlos = (vx*los[0] + vy*los[1] + vz*los[2])/(np.sum(1.*los*los)**0.5)
    
    vlos = vlos.convert_to_units('km/s')
    
    return vlos#/(1000.0*100.0) # convert from cm/s to km/s

def find_absorber_position_velocity(pf, ray, cz, abs_N, pec_vel=True):
    """
        Given the pf, yt ray object of the line of sight, and the fitted
        velocity of the absorber, this function gives the most likely physical
        cell associated with that given absorber... A given absorber will 
        be cause by gas over many cells, and possibly spread over cells not in
        the line of sight... this then really just picks the highest density
        region along the line of sight of the identified absorber, but not
        the center of mass or even highest density point of the whole 
        neutral gas region....
    """
    params = {}
    
    # length of ray in Mpc
    start, end = ray.start_point, ray.end_point
    los        = end - start
    ray_sort = np.argsort(ray['t'])


    l = ((start[0] - end[0])**2 + (start[1]-end[1])**2 + (start[2]-end[2])**2)**0.5
#    l = l * pf['Mpc'] # length in Mpc
    l = l.convert_to_units('Mpc')   

    H_o = pf.hubble_constant * 100.0 * yt.units.km / yt.units.s / yt.units.Mpc
    vH   = (ray['t'][ray_sort] * l) * H_o
    vH   = vH.convert_to_units('km/s')


    if pec_vel:
        vlos = losVel(los, ray)
    else:
        vlos = np.zeros(np.size(vH))
#    vH   = (ray['t'] * l) * pf.hubble_constant * 100.0 # in km/s
    vtot = vH + vlos
    
    ray_N = ray[('gas','H_number_density')][ray_sort] * ray['dx'][ray_sort] # cm^-2
    ray_N = ray_N.convert_to_units('cm**-2')
  #  ray_vx = ray['x-velocity']
  #  ray_vy = ray['y-velocity']
  #  ray_vz = ray['z-velocity']
  #  ray_vmag = 
    # the physical position of the absorber is going to be the cell that is
    # closest to the fitted cz of the absorber....
    # it may be best, however, to pick the top 5 cells (arbitrary number), then
    # take whichever of those 5 has the highest HI column ..... 
    
    absorber_index = find_absorber_index(ray_N, abs_N, vtot, cz, method='vselect')
    

    
   # if np.size(absorber_index) >= 1:
    try:
        try:
            absorber_index = absorber_index[0][0]
        except:   
            absorber_index = absorber_index[0]
    except:
        absorber_index = absorber_index
    
#    absorber_index = np.where( ray_N == ray_N[diffN.argsort()[0]] ) 
        
#    absorber_index = np.where( ray_N == np.max(ray_N[nearest]) )

    absorber_vH    = vH[absorber_index] # found the physical (vH) of the abs
    absorber_dist  = absorber_vH / H_o
    absorber_dist = absorber_dist.convert_to_units('Mpc')
   

    print absorber_index, absorber_dist, l, absorber_vH, vtot[absorber_index], cz, abs_N
   
    # find the x y z of the absorber
    pos = los * (absorber_dist / l) + start
    
    
    vel = np.array([ ray['x-velocity'][ray_sort][absorber_index].convert_to_units('cm/s').value,
                     ray['y-velocity'][ray_sort][absorber_index].convert_to_units('cm/s').value,
                     ray['z-velocity'][ray_sort][absorber_index].convert_to_units('cm/s').value])
                  
    
    params["Velocity"] = vel
    
    additional_fields = ["Temperature", "Metallicity", "Density",
                         "HI_NumberDensity", "CellVolume", 'NumberDensity',
                         'Pressure']#, 'MachNumber']
    ytnames = {'Temperature': ('enzo','Temperature'),
                                 'Metallicity': ('gas', 'metallicity'),
                                 'Density'    : ('enzo', 'Density'),
                                 'HI_NumberDensity': ('gas', 'H_number_density'),
                                 'CellVolume': ('index','cell_volume'),
                                 'NumberDensity': ('gas','number_density'),
                                 'Pressure' : ('gas','pressure')}
#                                 'MachNumber', : ('gas', 'mach_number')}
                                 
    ytunits = {'Temperature': 'K', 'Metallicity': 'Zsun', 'Density': 'g/cm**3',
               'HI_NumberDensity' : 'cm**-3', 'CellVolume' :'cm**3',
               'NumberDensity': 'cm**-3', 'Pressure' : 'dyne/cm**2'}
               
    
    for field in additional_fields:
        params[field] = ray[ytnames[field]][ray_sort][absorber_index].convert_to_units(ytunits[field]).value
        
    params['MachNumber'] = ray[('gas','mach_number')][ray_sort][absorber_index].value
    
    #params["Temperature"] = ray["Temperature"][absorber_index][0]
   # params["Metallicity"] = ray["Metallicity"][absorber_index][0]
   # params["Density"] = ray["Density"][absorber_index][0]
   # params["HI_Density"] = ray["HI_Density"][absorber_index][0]
    
    params['vh'] = absorber_vH.value
    params['vpec'] = vlos[absorber_index].value
    
    params['los'] = los.value
    
    return pos, params

def find_absorber_index(ray_N, abs_N, ray_v, abs_v, method='old'):

    if method == 'old':

        diff = np.abs(vtot - cz)
        nearest = diff.argsort()[:20] # may want to play around with ##....
        
        diffN = np.abs(ray_N - abs_N)

        # this is pretty gross, should rewrite to something nicer
        
        absorber_index = np.where( diffN == diffN[nearest][diffN[nearest].argsort()[0]]  )


    elif method == 'new':
        # calculate error for both... play around with 'minimum' error
    #    ray_N, abs_N = np.log10(ray_N), np.log10(abs_N)

        err_N = np.abs( (ray_N - abs_N)/ (abs_N) )
        err_v = np.abs( (ray_v - abs_v)/ (abs_v) )
        
        # try multiply first
        total_err = err_N + err_v
        
        absorber_index = np.argmin(total_err)
        
        
    elif method == 'distance':
    
      #  ray_N, abs_N = np.log10(ray_N), np.log10(abs_N)
        distance = (((ray_N-abs_N)/abs_N)**2 + ((ray_v-abs_v)/abs_v)**2)**0.5
        
        absorber_index = np.argmin(distance)
    
    elif method == 'vselect':
        vcut = 75.0 # within +/- 50 km/s
        
        # select nearest xx in velocity space
        shifted_v = ray_v - abs_v # cz is now v=0
        
        outside_range = np.abs(shifted_v.value) >= vcut
    
        ray_N[outside_range] = 1.0E50 / (yt.units.cm**2) # shift ray_N to silly values
        
        abs_N = abs_N.convert_to_units('cm**-2')
        
        err_N = np.abs( (ray_N.value - abs_N.value) / abs_N.value )
                
        absorber_index = np.argmin(err_N)
    
    
    return absorber_index



def plot_gas_3D(cluster,outdir,data, draw_velocity=True, pec_vel=True):
    """
        Function takes in the relevant pf, output directory, QSO
        data, virial radius, and center of the cluster in order
        to plot each fitted absorber in 3D in the cluster (0,0,0) is cluster
        center. The points are colored by the fitted HI column density,
        and a line is drawn indicating the velocity of the absorbing gas.
    """
    pf = cluster.ds
    center = cluster.center
    R_vir   = cluster.R_vir

    num_abs = np.size(data['QSO_id'])
    
    # define list to contain absorber objects
    all_absorbers = [None] * num_abs
    
    draw_los_arrow = True
    
    xpos,ypos,zpos,all_N = [],[],[],[]
    
    if num_abs > 40:
        draw_los_arrow = False
    # limits on the color scale  
    vmin, vmax = 1.0e9, 1.0e19
   # vmin, vmax = np.log10(vmin), np.log10(vmax)
                  
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')                
                
    # for each fittted line, compute the absorber location
    start_pos = np.array([data['x0'],data['y0'],data['z0']]).transpose()
    end_pos   = np.array([data['x1'],data['y1'],data['z1']]).transpose()
    
    i = 0
    for i in np.arange(num_abs):
        start = np.array([data['x0'][i], data['y0'][i], data['z0'][i]])
        end   = np.array([data['x1'][i], data['y1'][i], data['z1'][i]])
    
        # absorber information:
        N, b, z = data['HI_N'][i], data['b'][i], data['z'][i]
        cz      = z * 2.998E5 # abs velocity in km/s
        
        QSO_id  = data['QSO_id']
        
        # line of sight vector
        los   = end - start
    
        # make the ray
        ray = pf.h.ray(start, end)            
    
        # find absorber position along the ray
        pos, params = find_absorber_position_velocity(pf, ray, cz*yt.units.km/yt.units.s,
                                                      N/(yt.units.cm**2), pec_vel=pec_vel)
        
        
        # recenter so that the center of the cluster is at [0,0,0]
        # assuming pos and center are both in code length
        pos = pos.value - center
        
        # convert to Mpc then to virial radii 
        #pos = pos * pf['cm'] #/ R_vir
        pos = pos * pf.domain_right_edge[0].convert_to_units('cm').value

        # convert velocity to axis units
        params["Velocity"] = params["Velocity"]# / pf['cm'] * pf['Mpc'] / R_vir
        params['HI_Column_Density'] = N
        params['b'] = b
        params['cz'] = cz
        params['QSO_id'] = QSO_id

        print "absorber obs %.3E -- %.3E"%(params['HI_Column_Density'], params['HI_NumberDensity']*params['CellVolume']**(1.0/3.0))

        # draw the arrows
        # t is set right now at a Myr / 100.0
        
       
        # Store the results in list of absorber objects:
        all_absorbers[i] = Absorber(pos, params)      
 
 
        
   

    return all_absorbers
    # draw the arrow



