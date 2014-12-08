import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors,rc
import glob

from mpl_toolkits.axes_grid1 import make_axes_locatable

def heatmap2D(x, y, C, xlabel='', ylabel='', clabel='',\
                mesh_data = False, ax = None, vmin = None, vmax = None):
    """
    
    mesh_data : optional
        If data is not meshed, mesh it. Requires x and y to be single dim.
        arrays, and C to be single dim. array.
    
    """
#    if np.size(x) == np.size(C) :

    if ax is None:
        ax = plt.gca()
   
    if vmin is None:
        vmin = np.min(C)
    if vmax is None:
        vmax = np.max(C)


    if mesh_data == 0:
        N = (np.size(C))**0.5

        if not (np.size(x) == N + 1) or not (np.size(y) == N + 1) :
            print "ERROR. WHEN USING mesh_data = True , DIMENSIONS OF X" 
            print "AND Y MUST EACH BE A ONE DIMENSIONAL ARRAY WITH ONE MORE"
            print "ELEMENT THAN C"
        
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        
        x = np.append(x,x[-1] + dx)
        y = np.append(y,y[-1] + dy)
        
        x = x - 0.5*dx
        y = y - 0.5*dy
        
        xmesh, ymesh = np.meshgrid(x,y)    
        
        C = np.array(C).reshape((N,N))
    elif mesh_data == 1:
        N = (np.size(C))**0.5
        
        C = np.array(C).reshape((N,N))
        xmesh, ymesh = np.array(x).reshape((N+1,N+1)), np.array(y).reshape((N+1,N+1))
    else:
        xmesh , ymesh = x , y





    #fig = plt.figure(figsize=([6,6]))
    #ax  = fig.add_subplot(111)
    fig = plt.figure(figsize=([6,6]))
    ax = fig.add_subplot(111)
    ax.set_axis_bgcolor("#000000")
    
    plt1 = ax.pcolormesh(xmesh,ymesh, C, cmap='spectral', vmin = vmin, vmax = vmax)
    
    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(np.min(y),np.max(y))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right",size="2.5%",pad=0.001)
    plt.tight_layout()
    plt.colorbar(plt1,label=clabel,cax=cax)
    
    return plt1 

    

