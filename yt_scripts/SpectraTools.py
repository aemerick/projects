import numpy as np
import astropy.constants as C  

def _generate_species():
    """
    Populates a dictionary containing the constants for species used to 
    calculate equivalent width.
    """
    speciesList = {'Lya':\
                        {'f': 0.4164,        # oscillator strength  
                         'glambda' : 7616.0, # cm/s
                         'wl0'     : 1215.7, # Angstroms
                        }\
                  }
                  
    return speciesList
    

def Wlambda(NArray, bArray, species=None,speciesDict=None):
    """
    Given paramaters, calculate the equivalent width of 
    a given line.
    
    
    species:
        Can specify a species to load pre-calculated constants (oscillator 
        strength, gamma, etc.). Currently supports: Lyman Alpha as "Lya"

    speciesDict:
        Can be used to specify constants for the species. Has the following
        values: . If both species and speciesDict are specified, speciesDict
        takes priority.
        
    
       
    
    
    """
    c = C.c.value * 100.0 # speed of light in cm/s
     
    if not speciesDict == None:
        #Do something here
        print "yay things"
    else: # else, check species name against known list
    
        speciesList = _generate_species() # make dict of species
        
        speciesDict = speciesList[species]
        
    Wl = np.zeros( np.size(NArray) )
    
    
    t0 = _tau_o(NArray,bArray,speciesDict)
    
    glambda = speciesDict['glambda']
    limit = 1.25393 # as defined in Draine 9.27
    for i in np.arange(np.size(Wl)):
        N   = NArray[i]  # cm^-2
        b   = bArray[i] * 1000.0 * 100.0 # km/s -> cm/s
        tau = t0[i]
        
        if tau <= limit:
            Wl[i] = (np.pi)**0.5 * b * tau / (c * (1.0 + tau/(2.0*(2.0**0.5))))
        else: 
            Wl[i] = (\
                    (2.0*b/c)**2 * np.log(tau/np.log(2.0)) +\
                    b*glambda*(tau-limit)/(c*c*(np.pi)**0.5)\
                    )**0.5
    

    return Wl*speciesDict['wl0']*1000.0 # returns wlambda
        
def _tau_o(NArray,bArray,speciesDict):
    """
    Calculate tau o as per Draine 9.8 - 9.10. Ignoring the correction for 
    stimulated emission, as valid only for non-radio transitions.
    """
    
    # convert things to the right units and pull from speciesDcit
    f = speciesDict['f']
    wl0 = speciesDict['wl0'] * 1.0E-8 # converts A -> cm
    bArray = bArray * 1000.0 * 100.0  # converts km/s -> cm/s
    
    const = 1.497E-2 # constant. cm^2/s

    to = const * NArray*f*wl0/bArray

    return to
    
    
def calcN(Wlambda,species,b=30.0):
    """
    Calculates the column density of a given line provided the 
    equivalent width and species. For saturated lines, the doppler broadening
    value is important. If no b values is provided, b is assumed 30 km/s
    
    Parameters
    ----------
    W  :
        Equivalent widths in mA
    species : string
        Species name. Currently supports "Lya"
    b : optinoal
        Doppler broadening values. Must be an array of length equal to that
        of the equivalent width array. If not, b is assumed to be 30 km/s.
        Default = 30 km/s
    
    """
    
    speciesList = _generate_species()
    speciesDict = speciesList[species]
    
    wo = speciesDict['wo'] # natural line wavelength in A
    f  = speciesDict['f']  # oscillator strength
    
    if np.size(Wlambda) != np.size(b):
        b = 30.0 # km/s
    
    
    def linear(Wl,b):
    
    
        return N
    
    def flat(Wl,b):
        
        return N
        
    def damped(Wl,b):
    
        return N
        
        
        
    cases = {'linear' : linear,
             'flat'   : flat,
             'damped' : damped}
             
    for i in np.arange(0,np.size(Wlambda)):
        Wl = Wlambda[i]
        W  = Wlambda / (1000.0*wo) # mA / 1000*A - > mA/mA
        
        if Wl :
            case = "linear"
        elif Wl :
            case = "flat"
        elif Wl : 
            case = "damped"
            
        N = cases[case](Wl,b)
        
    
    
    
    
    
    
    return N
