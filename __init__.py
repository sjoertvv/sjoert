"""
Collection of tools that I'm gathering/writing while doing astronomy/astrophysics
research with Python.

2010-2012, Sjoert van Velzen
"""
import time 
time0 = time.time()

print ''
from astropy.io import fits as pyfits
from astropy.io import ascii as asciitable
import astropy.wcs as pywcs
print 'base',time.time()-time0

import string2radec

import stellar
print 'stellar',time.time()-time0
import io
print 'io',time.time()-time0
import rec
print 'rec',time.time()-time0
import dirs
print 'dirs',time.time()-time0
import html
print 'html',time.time()-time0
import simstat
print 'simstat',time.time()-time0
import plot
print 'plot',time.time()-time0
import misc
print 'misc',time.time()-time0
import catalogs #should become seperate rep
print 'catlogs',time.time()-time0
import latex
print 'latex',time.time()-time0
import sdss
print 'sdss',time.time()-time0
import fits
print 'fits',time.time()-time0
import cutout
print 'cutout',time.time()-time0


__all__ = ['io','rec','dirs','cutout','html', \
			'stellar','catalogs','simstat','plot',\
			'misc', 'latex', 'sdss', 'fits']
