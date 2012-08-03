'''
function to interface with NED (and DSS)
you can get an SED, currently only via the source name

the html pages are downloaded to NEDdir, by default
this is sjoert.dirs.catdir+NED/

2011 - Sjoert
'''

import sjoert
import os
import re
import numpy as np
import string2radec
import time
import pyfits

import shlex, subprocess


# change this variable to change the default input to all NED functions
# or give NEDdir explicitly to each call
NEDdir = sjoert.dirs.catdir+'NED/'


#NEDurl = '"http://ned.ipac.caltech.edu/cgi-bin/nph-datasearch?search_type=Photo_id&objid=80631&objname=_SOURCE_&img_stamp=YES&hconst=73.0&omegam=0.27&omegav=0.73&corr_z=1&of=table"' # doesn't work because of objid key
NEDurl = '"http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=_SOURCE_&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=NO"'

def get_NED_name(name, NEDdir=NEDdir):
    '''
    download or read NED page, return NED name of the galaxy
    >> get_NED_name(name, NEDdir=NEDdir)
    '''
    local_name = name.strip().replace(' ', '_')
    fname  = NEDdir+local_name+'.html'
    if not(os.path.isfile(fname)):
        print fname
        print 'not on disk, sendig query to NED'
        download_NED(name, local_name, NEDdir=NEDdir)
        time.sleep(0.5)
        
    return read_NED_mainpage(local_name)

def read_NED_mainpage(page_name, verbose=False, NEDdir=NEDdir):
    '''
    return NED name of object (more info can be added later)
    '''
    page = open(NEDdir+page_name+'.html')
    
    lines = page.readlines()
    i = 0
    name = None
    while i<len(lines):
        if re.search('INDEX for', lines[i]):
            if verbose: print lines[i]
            name = lines[i].split('INDEX for')[1].split('</td')[0].strip()
        i+=1

    return name

def get_NED_SED(name, NEDdir=NEDdir):
    '''
    download or read NED SED page
    out is frequency (Hz) and flux (Jy)
    >> nu, flux = get_NED_SED('M87')
    '''
    
    local_name = name.strip().replace(' ', '_')
    SEDpage_name = local_name+'_SED'
    
    if not(os.path.isfile(NEDdir+SEDpage_name+'.html')):

        download_NED(name, local_name, SEDpage_name, NEDdir=NEDdir)
        
    nu, flux = read_NED_SEDpage(SEDpage_name, NEDdir=NEDdir)
    return nu, flux

def download_NED(name=None, local_name=None,
                 SEDpage_name=None, NEDdir=NEDdir, debug=False):
    '''
    download summary page and SED page from NED
    >> download_NED(name=NVSS 143714-324447, local_name=None)

    may want to add add ra=, dec= options later
    '''

    if not(local_name):
        local_name=name.replace(' ','_')
    #if not(SEDpage_name):
    #     local_name+'_SED'

    # get the main NED page
    html_name = name.strip().replace(' ','_')
    command_line = 'wget '+NEDurl.replace('_SOURCE_', html_name)+ \
          ' -O '+NEDdir+local_name+'.html'
    args = shlex.split(command_line)

    if debug:
        print 'pasing to command line  :\n', command_line,' \n\n'
        print  args
        print 'NEDdir:', NEDdir
        
    #proc = subprocess.Popen(args cwd=NEDdir, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    p = subprocess.call(command_line, shell=True)
    

    # read the url to the SED page
    page = open(NEDdir+local_name+'.html','r')
    lines = page.readlines()
    i = 0


    if not SEDpage_name: return

    while i<len(lines):

        if debug: print lines[i]

        # get the link to the SED page
        if (re.search('Photometric', lines[i])):

            link = '"http://ned.ipac.caltech.edu/'+ \
                   lines[i].split('REF=')[1].split('TARGET')[0].strip()+'"'

            com = ' wget '+link+' -O  '+NEDdir + SEDpage_name+'.html'
            print 'pasing to command line  :\n', com, ' \n\n'
            os.system(com)
            return

        i+=1
        
        
    print 'WARNING No link to SED page found on NED main page.'


def read_NED_SEDpage(SEDpage_name, NEDdir=NEDdir):
    '''
    read the download SED page of NED
    return nu(Hz), flux(Jy)
    '''
    page = open(NEDdir+SEDpage_name+'.html')

    nu = []
    flux = []

    lines = page.readlines()
    i = 0

    while i<len(lines):
        if (re.search('A NAME="No', lines[i])):
            i+=3
            lflux = lines[i]
            #print ' line: ', lflux
            # check if not an upper limit
            if (re.search('=', lflux)):                   
                lflux = lflux.split('=')[1].split('W')[0].strip()
                flux.append(float(lflux) *1e23/1e4*1e7)
                #print 'flux (SI, Jy):',float(lflux), flux
                i+=3
                lnu = lines[i]
                #print ' line:', lnu
                lnu = lnu.split('=')[1].split('H')[0].strip()
                nu.append(float(lnu))
                #print 'nu (Hz)', nu
            
        i+=1
        
    return np.array(nu), np.array(flux)


def get_dss(ra, dec, size=1/60., name=None, sdir=None, color='red',
            debug=False, noclobber=False):
    '''
     wrapper to get the DSS image, requires that dss2 application is installed
    (see http://archive.eso.org/cms/tools-documentation/the-eso-st-ecf-digitized-sky-survey-application)
    >> fitsim = get_dss(ra, dec, size=1/60., name=None, sdir=None, color='red')
   
    input:
     ra, dec, size in degree
     sdir: location to save image, default is sjoert.catdir/DSS/fits/
     name: image name default is name dss+pl_Id (pl_Id=dss plate naming convention)
     color: blue, red, IR
    '''

    change_name = True
    if name is None:
        name = 'dss'
        change_name = False
    name = name
    
    if sdir is None:
        sdir = sjoert.dirs.catdir+'DSS/fits/'
        indir = sjoert.dirs.catdir+'DSS/infiles/' #keep log somewhere else
    else:
        indir=sdir

    new_name = sdir+name+'.fits'
    if noclobber and os.path.isfile(new_name):
        print 'not rerunning, reading:', new_name
        hdu = pyfits.open(new_name)
        return hdu


    infile_name = indir+name+'.in'
    logfile_name = indir+name+'.log'
    errfile_name = indir+name+'.err'

    infile = file(infile_name,'w')

    sradec = string2radec.radec2string(ra, dec)
    sra = sradec[0:11]
    sdec = sradec[11:]

    long_sra = sra[0:2]+' '+sra[3:5]+' '+sra[6:]
    long_sdec = sdec[0:3]+' '+sdec[4:6]+' '+sdec[7:]
    long_ssize = str(size*60)+' '+str(size*60)

    line =name+'  '+long_sra+' '+long_sdec+ ' '+long_ssize        

    infile.writelines(line+ ' \n')
    infile.close()
    
    command_line = 'cd '+sdir+' ; dss2 '+color+' -i '+infile_name +' > '+logfile_name #+' 2> '+ errfile_name

    command_line = 'dss2 '+color+' -i '+infile_name #+' > '+logfile_name #+' 2> '+ errfile_name

    if debug:
        print sradec
        print line
        print command_line

    #p = subprocess.Popen(command_line, shell=True)
    #t = os.popen(command_line)
    #p = subprocess.call(command_line, shell=True)

    os.system('touch '+logfile_name)
    logfile = open(logfile_name, 'r')

    args = shlex.split(command_line)
    if debug: print  args
    #p3 = subprocess.call('dss2 red -i /Users/sjoertvanvelzen/Documents/project/catalogs/DSS/fits/sjoert.in', shell=True)
    proc = subprocess.Popen(args, cwd=sdir, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    sleep_time = (size*60)
    if sleep_time <10: sleep_time = 10
    print 'sleeping for (s):', sleep_time
    time.sleep(sleep_time )
    # at this point we either have an image, or we give up
    proc.terminate()
    
    if proc.wait() < 0:
        
        print 'problem with dss call'
        print 'id:', proc.pid
        #proc.kill()
        #os.kill(proc.pid,  signal.SIGHUP) # doesn't work ?
        return None
    #time.sleep(1)
    #proc.kill()


    #stdout_value = proc.communicate()[0]
    stderr_value = proc.stderr.read()
    print 'stderr:', stderr_value
    stdout_value = proc.stdout.read()


    if debug:
        print '\n stdout: ', stdout_value

    lines = (stdout_value).split('\n')

    if not(len(lines)):
        print 'empty log file', logfile_name
        return None

    i=0
    while not re.search("Field", lines[i]) and (i < (len(lines)-1)):
        i+=1

    if i == (len(lines)-3):
        print 'no info from log file:',logfile_name
        return None

    info_line = lines[i+2].split()
    if len(info_line) < 7:
        print 'no info on plate id:', info_line
        return None
    
    Pl_Id = lines[i+2].split()[8].swapcase()
    if debug: print 'Pl_Id:', Pl_Id
    force_name = sdir+name[0:20]+Pl_Id+'.fits' # max 20 characters...
    if debug: print 'dss image name:', force_name
    
    if change_name:
        os.system('mv '+force_name+' '+new_name)
        if debug: print 'new name:', new_name

        try:
            hdu = pyfits.open(new_name)
            print 'created:', new_name
        except IOError:
            print 'file not created:', new_name
            return None

    else: 
        try:
            hdu = pyfits.open(force_name)
            print 'created:', force_name
        except IOError:
            print 'file not created:', force_name
            return None
        
    return hdu
