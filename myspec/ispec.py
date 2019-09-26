# Reads in a fits file to arrays and does a basic plot
# use reload(module) to refresh while editing
# ipython --pylab ORR from pylab import *
# plt.plot()
# plt.xlim() plt.ylim()


# CLEAN PLEASE:
import sys
import pyfits # to read in python
import matplotlib
import matplotlib.pyplot as plt # for plotting
import numpy as np
import os.path
from scipy.optimize import curve_fit
from pylab import *


## For latex in plotting
from matplotlib import rc
rc('text',usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
##

#fitsfile = str(sys.argv[1]) # Fits file taken as command argv

c = 0
wv_lya = 0
v_lya = 0
z_coma = 0 
z_virgo = 0
s_coma = 0
units = {}

def paramDefine():
    """
    Function to define some convenient numbers
    """
    global c, wv_lya, v_lya, z_coma, z_virgo, s_coma, units
    
    c       = 299792458.0*100.0 # cm/s
    wv_lya  = 1215.668          # Angstroms
    v_lya   = 1215.668*1.0E-12/c        # Hz
    z_coma  = 0.0231
    z_virgo = 0.0038
    
    s_coma  = 1000.0 # km/s

    units = { c : "cm/s", wv_lya : "Angstroms", v_lya : "Hz",\
              z_coma : "n/a", z_virgo : "n/a", s_coma : "km/s"}


#DECLARING GLOBAL VARIABLES
# Load some preset variables to make me have to type less... 
# not terribly versatile, but works...
#if 'INIT_CHCK' in globals():    # if the module is loaded once already, does not erase vars
#	print "Variables definedalready , not re-defined"
#else: 				# defines vars if they are not already there
	#fx    = [0]*100		# flux list
	#err   = [0]*100		# error list
	#w     = [0]*100		# wavelength list
	#fitsfile = "NONE"		# fits file name
	#INIT_CHCK = "initialized"	# dummy var for initialization (see above)

#----------------------------------------------------------------
#READS FITS FILE INTO DATA VARIABLES
def readFits(x):
#	global fx, err, w, fitsfile	
	fitsfile = x		
	fx    = np.array([])		# flux list
	err   = np.array([])	    # error list
	w     = np.array([])  	# wavelength list


	if os.path.isfile(x):	
		hdulist = pyfits.open(x)
		print "Fits file opened!"
        else:
		print "FITS FILE DOES NOT EXIST... CHECK FILENAME"
	
	#Import data to arrays defined at the global levels
	if len(hdulist)==2: # one of the coadd files probably
 		tbdata = hdulist[1].data
		cols = hdulist[1].columns

		if 'WAVELENGTH' in cols.names:
			w = tbdata.field('WAVELENGTH')
		elif 'WAVE' in cols.names:
			w = tbdata.field('WAVE')
		else:
			print "NO WAVELEGNTH DATA FOUND"

		if 'ERROR' in cols.names:
			err = tbdata.field('ERROR')
		else:
			print "NO ERROR DATA FOUND"
	
		if 'FLUX' in cols.names:
			fx = tbdata.field('FLUX')
		else:
			print "NO FLUX DATA FOUND"
	
	elif len(hdulist)==3: #one of the xpsec
		fx = hdulist[0].data # Flux
		err = hdulist[1].data # Error
		w = hdulist[2].data # Wavelength bins

	elif len(hdulist)>3:
		print "CHECK FITS FILE HEADERS... DIFFERENT FORMAT THAN NORMAL"


	print "Data loaded to 'f', 'err', and 'w' lists"
	print "Fits filename saved to global 'fitsfile' "
	return w, fx, err
	hdulist.close()
# END READ DATA FUNCTION
#-----------------------------------------------------------------

def shift_w(w, z):
    """
    Shift to different reference frame
    
    """
    
    return w/(1.0+z)

# Automatically plots and spits out .eps files
# Takes in data in the x and y arrays.
# over until the entire spectrum is output to file.
def plotSpec(x_orig,y_orig,xlow,xhigh,fitsname,*arguments,**keywords):  
    print "WTF"
    if 'persist' in keywords.keys():
        print 'persisting'
    else:
        clf()
    print "qwefqwefqwef"

    if len(x_orig) < len(x_orig[0]):
        x = x_orig[0]
    if len(y_orig) < len(y_orig[0]):
        y = y_orig[0]
    #x = x_orig
    #y = y_orig
	#fitsname = fitsfile.replace('.fits','')
#	if keyword_args.has_key('IDfile'):
#		idfile = keyword_args['IDfile']
#		with open(idfile,'r') as lines
#			for line in lines

    if 'auto' in keywords.keys():
        print "in auto"
        if keywords['auto'] == True:
            if keywords.keys():
                DX = keywords['dx']
            else:
                DX = 20
            print "here"
            xmin = x[0]
            xmax = x[-1]
            print xmin, xmax
            xlow  = xmin      # min of xrange
            xhigh = xmax + DX # max of xrange
            while xlow < xmax:
                plt.plot(x,y,color="black",linewidth=0.75)
                plt.xlabel(r'Wavelength ($\AA$)')
                plt.ylabel(r'Flux ($ergs s^{-1} cm^{-2} \AA^{-1}$)')
                plt.xlim(xlow,xhigh)
                axes = plt.axes()
                plt.savefig(fitsname + "_" + str(int(xlow)) + "-" + str(int(xhigh)) +'.png')
                xlow = xlow + (DX - 4)
                xhigh = xlow + DX
	
    else:
			
        plt.plot(x,y,color="black",linewidth=0.75,drawstyle="steps")
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'Flux (ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
        plt.xlim(xlow,xhigh)
		
        if xlow != x[0]:
            if xhigh != x[-1]:
                dwv    = x[1]-x[0]
                binl = int((xlow - x[0])/(dwv))
                binh = int((xhigh - x[0])/(dwv))

                yavg = sum(y[binl:binh]) / float(len(y[binl:binh]))
                ystd = np.std(y[binl:binh])
                ymax = yavg + 6*ystd
                ymin = yavg - 4*ystd	
                plt.ylim(ymin,ymax)
		
        increment = (xhigh-xlow)/5.0
        xticks = np.arange(xlow,xhigh+2,increment)
        plt.xticks(xticks)	
        axes = plt.axes()
		# insert stuff for file here

		#if xlow < 1220.1080 & xhigh > 1220.1080:
		#	plt.(1220.1080,axes[3]*0.75,color="red")

#END PLOT FUNCTION







#old
def plotSpectra(x,y,xmin,xmax):
	if xmin < 0:
		xmin = x[0]
        if xmax < 0:
		xmax = x[-1]

	# use math to make nice-ish y max and mins
	yavg = sum(y) / float(len(y))
	ystd = np.std(y)
	ymax = yavg + 8*ystd
	ymin = yavg - 4*ystd


	# loop here from i=0 until i*shift + lambda_o > lambda_f
	plt.plot(x,y,color="black",linewidth=0.25) 
	plt.xlabel('Wavelength (Angstrom)')
	plt.ylabel('Flux')
	plt.xlim(xmin,xmax)
	plt.ylim(ymin,ymax)
	#plt.axis([xmin,xmax,ymin,ymax])
	plt.show() 

###########



def zlya(obs):
	lya = 1215.668
	c   = 299792.458

	z = (obs - lya)/lya
	vel = z*c

	print "z = " + str(z)
	print "vel = " + str(vel)

def lya(z,disp,**keywords):
	lya = 1215.668 #angstroms	
	c   = 299792.458 #km/s
	
	if 'obj' in keywords.keys():
		print 'yay'
		if keywords['obj'] == 'Coma':
			z = 0.0231
			disp = 1000 #km/s
		elif keywords['obj'] == 'Virgo':
			z = 0.0038
			disp = 1000

	zh = (z*c + disp)/c
	zl = (z*c - disp)/c

	central = (lya*z) + lya
	low     = (lya*zl) + lya
	high    = (lya*zh) + lya

	print "z   = " + str(z)
	print "vel = " + str(c*z)
	print str(zl*c) + ' - ' + str(zh*c)
	print "wvobs = " + str(central)
	print str(low) + " - " + str(high)



# in tha workssss
def fitGaus(xData,yData,xmin,xmax):
	# fits a gaussian to the given dataset
	mean = sum(xData*yData)
	sigma = sum(yData*(xData-mean)**2)	

	dx = xData[1]-xData[0]
	binMin = int((xmin - xData[0])/dx)
	binMax = int((xmax - xData[0])/dx)
	popt, pcov = curve_fit(gaussian,xData[binMin:binMax],yData[binMin:binMax], p0=[1,mean,sigma])
	print popt
	return popt


def gaussian(x, a, b, c):
	return a * np.exp(-(x - b)**2/(c**2))


