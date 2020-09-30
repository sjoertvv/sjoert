'''
this module is merely a container with directories:
catalogs, metadata, or projects
'''

import os


# the catalogs
home_folder = os.path.expanduser('~') 
catdir = home_folder+'/Documents/project/catalogs/'
metadir = home_folder+'/Documents/project/meta/'
fltdir = home_folder+'/Documents/project/meta/filters/'

nvdir = catdir+'NVSS/'
sudir = catdir+'SUMSS/'
xsdir = catdir+'2MASS/'
audir = catdir+'AUGER/'
fidir = catdir+'FIRST/'


# project specific folders
frankdir = home_folder+'/Documents/project/frank/'
frdir = frankdir+'Frank/'

# dir for latex scripts
latexdir= home_folder+'/Documents/latex/'

def mkdir(dir):
    if not(os.path.isdir(dir)):
           os.system('mkdir '+dir)
    return
