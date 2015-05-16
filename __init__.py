"""
Collection of tools that I'm gathering/writing while doing astronomy/astrophysics
research with Python.

2010-2012, Sjoert van Velzen
"""
import time 
#time0 = time.time()

from astropy.io import fits as pyfits
#from astropy.io import ascii as asciitable # slow (0.3 s)
import astropy.wcs as pywcs

import stellar
import simstat
import io
import rec
import dirs
import html
import plot
import misc
import catalogs
import latex
import sdss
import fits
import cutout


__all__ = ['io','rec','dirs','cutout','html', \
			'stellar','catalogs','simstat','plot',\
			'misc', 'latex', 'sdss', 'fits']
