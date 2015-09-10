'''
this module is merly a container with directories:
catalogs, metadata, or projects
'''

import os


# the catalogs
catdir = os.environ['HOME']+'/project/catalogs/'
metadir = os.environ['HOME']+'/project/meta/'
fltdir = os.environ['HOME']+'/project/meta/filters/'

nvdir = catdir+'NVSS/'
sudir = catdir+'SUMSS/'
xsdir = catdir+'2MASS/'
audir = catdir+'AUGER/'
fidir = catdir+'FIRST/'


# project specific folders
frankdir = os.environ['HOME']+'/project/frank/'
frdir = frankdir+'Frank/'

# dir for latex scripts
latexdir= os.environ['HOME']+'/latex/'

def mkdir(dir):
    if not(os.path.isdir(dir)):
           os.system('mkdir '+dir)
    return
