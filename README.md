Collection of tools that I'm gathering while doing astronomy/astrophysics
research with Python. For example,

 - Create, merge, append with numpy.recarrays (`rec.py`)
 - Make latex tables from (record) arrays (`latex.py`)
 - Convert between magnitude and luminosity (`stellar.py`)
 - I/O convenience for fits/ascii/json tables  (`io.py`)
 - Convert MJD to anything (`simtime.py`)
 - Binning, e.g. to make luminosity functions (`simstat.py`)


Requirements: 

 - numpy
 - astropy - http://www.astropy.org/

Optional:
 - matplotlib
 - scipy
 - APLpy  - http://aplpy.github.com/
 - k3match - http://pschella.github.com/k3match

To use, at your own risk, run `python setup.py install` (or place this folder in your $PYTHONPATH)
