"""
Collection of tools that I'm gathering/writing while doing astronomy/astrophysics
research with Python.

2010-2018, Sjoert van Velzen
"""
from __future__ import print_function 

import time 

from astropy.io import fits as pyfits
from astropy.io import ascii as asciitable 
import astropy.wcs as pywcs

from . import stellar
from . import simstat
from . import io
from . import rec
from . import dirs
from . import html
from . import plot
from . import catalogs
from . import latex
from . import sdss
from . import fits
from . import cutout

#from . import misc

__all__ = ['io','rec','dirs','cutout','html', \
			'stellar','catalogs','simstat','plot',\
			'latex', 'sdss', 'fits']
