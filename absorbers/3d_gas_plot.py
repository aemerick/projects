### The purpose of this script is to 

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import plotting.scatter3D as s3d # for arrow drawing
from yt.mods import *

## mencoder "mf://*.png" -mf w=800:h=600:fps=2:type=png -ovc copy -oac copy -o movie.avi
    
to_directory        = 'SP/32x32/z_3d/'
pfname = "RD0020/RD0020"
line_file_name = to_directory + "QSO_data_absorbers.out"

center     = np.array([0.49730353,  0.50146712,  0.5017492]) # in code units
R_vir  = 1.85 # Mpc
vel_arrows = False

DEFAULT_ARROW_TIME = 3.15E13 / 75.0

class GalCluster:
    def __init__(self, pf, center, R_vir):
        self.pf = pf
        self.center= center
        self.R_vir = R_vir # in Mpc!!



class Absorber:
    """ making a class for the absorbers to store data better for
        post processing after initial render """
    absorberCount = 0
        
    
    def __init__(self, pos, vel, N, plot = True):
        self.pos = pos          # 3d array giving position
        self.vel = vel          # 3d array giving velocity
        self.N   = N            # Fitted HI column
        self.plot = plot        # bool ... should this be plotted?
        Absorber.absorberCount +=1
        
    # Functions to set attributes if needed
    def setPos(self, pos):
        self.pos = pos
        
    def setVel(self, vel):
        self.vel = vel
    
    def setN(self, N):
        self.N = N
        
    def setPlot(self, plot):
        self.plot = plot

    # Calculates arrow end point for velocity    
    def arrowEndPoint(self, t = DEFAULT_ARROW_TIME ):
        """
        Calculates arrow end point
        """
        
        return t*self.vel + self.pos
    
        
class AbsorberPlot:
    """
    Class for the absorber plot. Initialize with a list of absorbers to plot
    """

    def __init__(self, absorbers, cluster):
        """
        Must have a list of absorber class objects and a galaxy cluster
        class object to initialize
        """
        self.absorbers = absorbers
        self.cluster = cluster
    

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

    def set_N_cut(self, Nmin, Nmax):
        self.Nmin = Nmin
        self.Nmax = Nmax

    def set_velocity_cut(self, vmin, vmax, absvel=False):
        a = self.absorbers[0]
        pf = self.cluster.pf
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
                    v = a.vel
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
                    v = a.vel
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
            
                if a.N > Nmax or a.N < Nmin:
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
            
                if a.T > Tmax or a.T < Tmin:
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
                if a.Z > Zmax or a.Z < Zmin:
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


    def plot_points(self):
        """
        Nmin, Nmax, vmin, vmax are the HI column and velocity
        cutoffs. Velocity cutoff is vel magnitude AND is in km/s.
        absvel=True set the velocity cutoffs as absolute velocity
        """
        cmin,cmax = 1.0E9,1.0E19
            
        xpos = []
        ypos = []
        zpos = []
        all_color   = []
        
        try:
            self.points.remove()
            print "Cleared previous ponits"
        except AttributeError:
            print "Making new plot"

        for a in self.absorbers:
 
            if a.plot:
                pos = a.pos

                xpos.append(pos[0])
                ypos.append(pos[1])
                zpos.append(pos[2])

                # this is bad... need to fix to generalize to take an arbitrary
                # color value decided by user at any point.
                if self.cparam == "HI_Column_Density":
                    all_color.append(a.N)
                elif self.cparam == "Metallicity":
                    all_color.append(a.Z)
                elif self.cparam === "Temperature":
                    all_color.append(a.T)

                                
        # draw the points... make the points mappable
        points = ax.scatter(xpos,ypos,zpos, c = np.log10(all_N),
                            s=50, cmap='spectral',
                            vmin=np.log10(cmin),vmax=np.log10(cmax))

        # set the points object
        self.points = points
        # redraw
        self.fig.canvas.draw()


    # draw a sphere of radius r in plot units and color        
    def drawSphere(self, r, color):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x=np.cos(u)*np.sin(v) 
        y=np.sin(u)*np.sin(v)
        z=np.cos(v)
        
        self.ax.plot_wireframe(r*x,r*y,r*z, color=color)
        
    def drawVelocityArrows(self, t=DEFAULT_ARROW_TIME):
        """
        Draw velocity arrows for all absorbers
        """
        # need to draw arrows all at once.... not one by one
        
        for a in self.absorbers:
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
    
    
  #  def rotate_plot(self, elevation, azim):
  #"""
 #   """
    
#        for az in azim:
            #self.ax.view_init(elev=elevation,azim=az)
            #self.fig.canvas.draw()
        
    


def losVel(los,ray):
    """
    Calculate the scalar line of sight velocity given the line of sight
    vector and the yt ray object. Does this for all cells along the ray
    """
    
    vx = ray['x-velocity']
    vy = ray['y-velocity']
    vz = ray['z-velocity']
    
    vlos = (vx*los[0] + vy*los[1] + vz*los[2])/(np.sum(los*los)**0.5)
    
    return vlos/(1000.0*100.0) # convert from cm/s to km/s

def find_absorber_position_velocity(pf, ray, cz):
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
    
    # length of ray in Mpc
    start, end = ray.start_point, ray.end_point
    los        = end - start

    l = ((start[0] - end[0])**2 + (start[1]-end[1])**2 + (start[2]-end[2])**2)**0.5
    l = l * pf['Mpc'] # length in Mpc
    
    vlos = losVel(los, ray)
    vH   = (ray['t'] * l) * pf.hubble_constant * 100.0 # in km/s
    vtot = vH + vlos
    
    ray_N = ray['HI_NumberDensity'] * ray['CellVolume']**(1.0/3.0) # cm^-2
    
    # the physical position of the absorber is going to be the cell that is
    # closest to the fitted cz of the absorber....
    # it may be best, however, to pick the top 5 cells (arbitrary number), then
    # take whichever of those 5 has the highest HI column ..... 
    diff = np.abs(vtot - cz)
    nearest = diff.argsort()[:8] # may want to play around with ##....
    
    absorber_index = np.where( ray_N == np.max(ray_N[nearest]) )

    absorber_vH    = vH[absorber_index][0] # found the physical (vH) of the abs
    absorber_dist  = absorber_vH / (pf.hubble_constant * 100.0) 
   
    
    # find the x y z of the absorber
    pos = los * (absorber_dist / l) + start
    
    vel = np.array([ ray['x-velocity'][absorber_index][0],
                     ray['y-velocity'][absorber_index][0],
                     ray['z-velocity'][absorber_index][0]])
                  
    
    return pos, vel


def plot_gas_3D(cluster,outdir,data, draw_velocity=True):
    """
        Function takes in the relevant pf, output directory, QSO
        data, virial radius, and center of the cluster in order
        to plot each fitted absorber in 3D in the cluster (0,0,0) is cluster
        center. The points are colored by the fitted HI column density,
        and a line is drawn indicating the velocity of the absorbing gas.
    """
    pf = cluster.pf
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
    
        # line of sight vector
        los   = end - start
    
        # make the ray
        ray = pf.h.ray(start, end)            
    
        # find absorber position along the ray
        pos, vel = find_absorber_position_velocity(pf, ray, cz)
        
        
        # recenter so that the center of the cluster is at [0,0,0]
        pos = pos - center
        
        # convert to Mpc then to virial radii 
        pos = pos * pf['Mpc'] / R_vir
    
        # add to lists
        xpos.append(pos[0])
        ypos.append(pos[1])
        zpos.append(pos[2])
        all_N.append(N)        

        # make an array to hold all the data point data and save it somehow
        # maybe build an absorber class with position, velocity, N, 
    

        # convert velocity to axis units
        vel = vel * 100.0 * 1000.0 / pf['cm'] * pf['Mpc'] / R_vir
    
        # draw the arrows
        # t is set right now at a Myr / 100.0
                    
        # Store the results in list of absorber objects:
        all_absorbers[i] = Absorber(pos, vel, N)      
        
        if draw_velocity:
            s3d.drawArrow(pos, ax, v=vel, t = DEFAULT_ARROW_TIME)
        
        # draw los vector arrow
        if draw_los_arrow:
            start = (start - center) * pf['Mpc'] / R_vir
            end   = (end   - center) * pf['Mpc'] / R_vir
            los = end - start
            
            los_arrow = np.array([start[0],start[1],start[2],end[0],end[1],end[2]])
            s3d.drawArrow(los_arrow,ax,color='r')
            
            
    
    points = ax.scatter(xpos,ypos,zpos, c=np.log10(all_N), s=50, cmap='spectral',
                           vmin=np.log10(vmin), vmax=np.log10(vmax))  
    
    

    # Draw a circle at 1 and 2 virial radii... red and blue respectively
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    r = 1.0
    x=np.cos(u)*np.sin(v) 
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    ax.plot_wireframe(r*x, r*y, r*z, color="r")
    
    r = 2.0
    ax.plot_wireframe(r*x, r*y, r*z, color="b")
    
    #lim = np.array([np.min([xlim,ylim,zlim]),np.max([xlim,ylim,zlim])])
    
    # set the x,y,z limits of the plot and add labels
        # center the plot in all 3 axes
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()
    lim = [-2.5,2.5]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_zlim(lim)


    
    # if lines of sight are not drawn for the absorbers,
    # then draw a line of sight arrow to indicate the observing direction
    # places this arrow at the edge of the plot
    if not draw_los_arrow:
        
        start = (start - center) * pf['Mpc'] / R_vir
        end   = (end   - center) * pf['Mpc'] / R_vir
        los = end - start
    
        for i in np.arange(np.size(los)):
            if los[i] == 0:
                start[i], end[i] = lim[0], lim[0]
    
        los_arrow = np.array([start[0],start[1],start[2],end[0],end[1],end[2]])
        s3d.drawArrow(los_arrow,ax,color='r')
    
    
    # Do some stuff to make the plot look nice
    # add axes labels
    ax.set_xlabel(r"R$_{vir}$")
    ax.set_ylabel(r"R$_{vir}$")
    ax.set_zlabel(r"R$_{vir}$")
    
    # adding color bar and color bar label
    cbar = fig.colorbar(points)
    cbar.set_label("log(N)")
    
    
    
    fig.savefig(outdir + "/3d_test.png")
    return fig, ax, points, all_absorbers
    # draw the arrow

# load the pf
pf = load(pfname)


# load the fitted lines data


cluster = GalCluster(pf,center,R_vir)

data = np.genfromtxt(line_file_name,names=True)

fig, ax, points, all_points = plot_gas_3D(cluster,to_directory,data,draw_velocity=vel_arrows)


aplot = AbsorberPlot(all_points,cluster)
aplot.set_plot(fig,ax,points)

#
    


