"""
Collection of tools that I'm gathering/writing while doing astronomy/astrophysics
research with Python.

2010-2020, Sjoert van Velzen
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
#from . import dirs # doesnt work on windows
from . import html
from . import plot
from . import online
from . import latex
#from . import sdss
from . import fits
from . import cutout
from . import gaia
from . import swift

#from . import misc

__all__ = ['io','rec','dirs','cutout','html', \
			'stellar','online','simstat','plot',\
			'latex', 'sdss', 'fits', 'gaia', 'swift']
