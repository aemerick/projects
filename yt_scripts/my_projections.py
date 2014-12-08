from yt.mods import *
import os.path
import numpy as np

def my_proj(pf,fields=None,weight_fields=None, sub_volume=None,axis="z",
            view_width = None,cmap_color="spectral", name='', save=True):
    """
    make some projections. Output to working directory in new folder named
    "PFNAME_proj"
    
    pf : 
        yt pf object
    fields : list, optional
        optional list of fields to project. Default is Density,
        Dark_Matter_Density, and Temperature
    weight_fields : list, optional
        optional list of fields to weight by. Default is Density
    sub_volume : optional
        optional sub volume of pf to use for the projection. Default is the
        entire box
    axis : optional
        optional axis to project along. Default is 'z'
    view_width: float, optional
        optional width of projection view. Default is full box
    cmap_color: str, optional
        optional cmap name to use. Default is spectral.. yt default is algea
    
    Returns
    -------
    projDict : dictionary
       A dictionary of dictionaries. The upper level dictionary is sorted by
       the weight field. This contains a dictionary of each projection made 
       associating field name with the yt projection object. This is made
       mainly in case some extra stuff needs to be done with the projectionss    
    """
    projDict = {}
    
    # use the default fields and weight fields 
    if fields==None:
        fields = ['Density', 'Dark_Matter_Density',
                  'Temperature']

    if weight_fields == None:
        weight_fields = ['Density']

    # default to using the entire volume        
    if sub_volume == None:
        center = [0.5,0.5,0.5]
        lc, rc = [0.0,0.0,0.0], [1.0,1.0,1.0]
        reg = pf.h.region(center,lc,rc)
    else:
        reg = sub_volume
        center = reg.center
        
    # default to full box width
    if view_width==None:
        view_width = pf['Mpc']
    
    # make the ouptut directory
    outDir = os.getcwd() + "/" +str(pf) +"_proj/"
    if not os.path.exists(outDir): os.makedirs(outDir)
   
   
    
    
    # loop over fields and weights
    for field_weight in weight_fields:
        temp_dict = {}

        for field in fields:
            proj = ProjectionPlot(pf, axis, field, weight_field=field_weight,
                                  width = (view_width, 'Mpc'),
                                  center = center, data_source=reg)
                                  
            proj.set_buff_size(2056)
            
            proj.set_cmap(field,cmap_color)
            
            if save: proj.save(os.getcwd() + "/"+str(pf)+"_proj/" + name)
            
            print "completed ", field, " weighted by ", field_weight
            
            temp_dict[field] = proj
            
        projDict[field_weight] = temp_dict
            
    # return 
    return projDict
    
    
def annotate_halo_mass(pf,proj,halos_list,axis,max_number=None,text_args=None,
                       virial_overdensity=100.0):
    """
        Need to pass center coordinates of image in full simulation units,
        along with the number of halos and the halos list. Pass the proj.
        Get the projection width to obtain images coordinates.
    """
     
    if axis == 'z':
        x1, x2 = 0, 1
    elif axis == 'x':
        x1, x2 = 1, 2
    elif axis == 'y':
        x1, x2 = 0, 2
    else:
        print "WARNING CHOSE AN INCORRECT PROJECTION AXIS"
        print "Please use 'x', 'y', or 'z' "
        return 0
  
    
    if max_number == None:
        max_number = np.size(halos_list)
        if max_number >= 20:
            print "Warning... annotating for ", max_number, " halos."
            print "Use max_number option to reduce this"    
    
    
    
    width  = proj.width[0]
    hwidth = 0.5 * width
    center = proj.center 
    
    if text_args == None:
        text_args = {'size': 'x-large'}
 
    
    for halo in halos_list[0 : max_number]:
        # get halo center of mass and cut extra coordinate
        COM = halo.center_of_mass()
        COM = np.array([ COM[x1], COM[x2] ])

        # pos vec and distance 
        vec  = np.array([center[0]-COM[0],center[1]-COM[1]])
        dist = (vec[0]**2 + vec[1]**2)**0.5
        
        # convert to image coordinates
        distIMG = dist / hwidth
        vecIMG  = vec / hwidth
        
        # find center of annotations relative to image center
        centText = np.array([0.5,0.5]) - vecIMG
        
        # if it is within the image, annotate!
        if all(centText) >= 0. and all(centText) <= 1.0:
            massStr = "%.4E" % (halo.virial_mass(virial_overdensity=virial_overdensity))
            radiusStr = "%.4F" % (halo.virial_radius(virial_overdensity=virial_overdensity)\
                                  *pf['Mpc'])
                   
            massStr = "M = " + massStr
            radiusStr = "R = " + radiusStr            
            proj.annotate_text(centText,massStr,text_args=text_args)
            centText = np.array([centText[0],centText[1]-0.05])
            proj.annotate_text(centText,radiusStr,text_args=text_args)
    
    
    return proj

