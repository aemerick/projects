import numpy as np

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
