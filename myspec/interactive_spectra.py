# Reads in a fits file to arrays and does a basic plot

import sys
import pyfits # to read in python
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt # for plotting
import numpy as np

#fitsfile = str(sys.argv[1]) # Fits file taken as command argv

#read a fits file into an array
def readFits(tbdata0,tbdata1,tbdata2,fitsfile):
	hdulist = pyfits.open(fitsfile)

	#Import data to arrays
	tbdata0 = hdulist[0].data # Flux?
	tbdata1 = hdulist[1].data # Error
	tbdata2 = hdulist[2].data # Wavelength bins


def plotSpectra(xmin,xmax,tbdata0,tbdata2):
	deltaLambda = abs(tbdata2[0] - tbdata2[1])
	Lambda_o = tbdata2[0]  #finds lowest wavelength bin
	Lambda_f = tbdata2[-1] #finds last wavelength bin

	yavg = sum(tbdata0) / float(len(tbdata0))
	ystd = np.std(tbdata0)
	ymax = yavg + 6*ystd
	ymin = yavg - 6*ystd

	# loop here from i=0 until i*shift + lambda_o > lambda_f
	plt.plot(tbdata2,tbdata0,color="black",linewidth=0.25) 
	plt.xlabel('Wavelength (Angstrom)')
	plt.ylabel('Flux')
	plt.axis(xmin,xmax,ymin,ymax)
	#plt.savefig('test' + str(i) + '.eps')
	plt.show() 
