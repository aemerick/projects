import sys
import pyfits # to read in python
import matplotlib
import matplotlib.pyplot as plt # for plotting
import numpy as np
import os.path
from scipy.optimize import curve_fit
from pylab import *

def plotSpec(x_orig,y_orig,xlow,xhigh,*arguments,**keywords):  
	if 'persist' in keywords.keys():
	   print 'persisting'
	else:
	    clf()

    
	if len(x_orig) < len(x_orig[0]):
		x = x_orig[0]
	if len(y_orig) < len(y_orig[0]):
		y = y_orig[0]


	if 'auto' in keywords.keys():
		if keywords['auto'] == 'True':
			if keywords.keys():
				DX = keywords['dx']
			else:
				DX = 20

			xmin = x[0]
			xmax = x[-1]

			xlow  = xmin      # min of xrange
			xhigh = xmax + DX # max of xrange
			while xlow < xmax:
				plt.plot(x,y,color="black",linewidth=0.75)
				plt.xlabel(r'Wavelength ($\AA$)')
				plt.ylabel(r'Flux ($ergs s^{-1} cm^{-2} \AA^{-1}$)')
				plt.xlim(xlow,xhigh)
				axes = plt.axes()
				plt.savefig(fitsname + "_" + str(int(xlow)) + "-" + str(int(xhigh)) +'.eps')
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





