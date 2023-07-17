'''
function to interface with PS1, NED, SDSS, DSS
you can get an SED, currently only via the source name

the html pages are downloaded to NEDdir, by default
this is sjoert.dirs.catdir+NED/

2011 - Sjoert started
2020 - adding neowise   
'''
from __future__ import print_function 

import os, sys, time, re
import numpy as np
import shlex, subprocess
import requests
import json

import sjoert.dirs as dirs
import sjoert.rec as rec
from sjoert.io import readascii, pyfits, json2rec
from sjoert.stellar import iau_name
import sjoert.simstat

import astropy.time
import astropy.io.ascii
from astropy import units as u
from astropy.coordinates import SkyCoord

import matplotlib
from matplotlib import pyplot as plt

selcols_str = 'selcols=ra,dec,sigra,sigdec,mjd,w1mpro,w1sigmpro,w1flux,w1sigflux,w1sky,w2flux,w2sigflux,w2mpro,w2sigmpro,w2sky,cc_flags,nb,saa_sep,qual_frame,w1mpro_allwise,w1sigmpro_allwise,w2mpro_allwise,w2sigmpro_allwise,w3mpro_allwise,w3sigmpro_allwise,w4mpro_allwise,w4sigmpro_allwise'

url_to_fill =\
    'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=neowiser_p1bs_psd&objstr={0}&radius=2.5&spatial=Cone&radunits=arcsec&'+selcols_str

def get_WISE(ra, dec, name='', t0=58119, wait=False, verbose=True, redo=False):
    
    '''
    query IPAC to get allWISE (and neoWISE)
    >> catalog, class_str = get_WISE(ra, dec, name='', t0=58119,  wait=False, verbose=True, redo=False)

    input:
        - name:  if given, we save the ipac table to this path 
        - t0=:   normalization for plotting
        - redo:  forve new download even if catalog name is found

    '''    


    if name and os.path.isfile(name+'-neoWISE.ipac'):
        
        if verbose:
            print ('reading:', name+'-neoWISE.ipac')
        
        astro_tab = astropy.io.ascii.read(name+'-neoWISE.ipac')

        if not('w2flux' in astro_tab.dtype.names):
            print ('old format, need to download again')
            redo = True
    else:
        redo = True

    if redo or not(name):
        obs_str = '{0:0.6f}+{1:0.6f}'.format(ra, dec)

        url = url_to_fill.format(obs_str)
        if verbose:
            print ('getting neoWISE data from:')
            print(url)

        try:
            r = requests.get(url, timeout=120)

            content=str(r.content)

            urlindex=content.find('/workspace/')
            urlindex2=content.find('.tbl')

            table=requests.get('https://irsa.ipac.caltech.edu'+content[urlindex:urlindex2]+'.tbl')

        except Exception as e:
            print ('plot_neoWISE: no connection to IRSA; url was:\n', url)
            print (e)
            return [], None

        astro_tab = astropy.io.ascii.read(table.text)

        # save the data
        if name:
            with open(name+'-neoWISE.ipac', 'w') as the_file:
                the_file.write(table.text)

    
    dd = np.array(astro_tab)

    if len(dd)==0:
        return ["neoWISE: no match"], None


    # apply some quality cuts
    iqual = (dd['saa_sep']>5) * (dd['qual_frame']>0) * (dd['cc_flags']=='0000')
    if verbose:
        print ('# neoWISE observations', len(dd))
        print ('# data point removed by quality cuts', sum(iqual==False))
    
    if sum(iqual)<3:
        if verbose:
            print ('no enough good data:', astro_tab[('qual_frame', 'saa_sep', 'cc_flags','w1mpro', 'w1sigmpro')])
        return ["neoWISE: not enough detections (N={0})".format(sum(iqual))], None

    dd = dd[iqual]

    xx = dd['mjd']-t0


    # attempt to make smart bins that are catch the entire 
    roundtime = np.unique(np.round(xx, decimals=-2)) # bins in 100 days
    neartime = [ xx[np.argmin(np.abs(xx-rt))] for rt in roundtime ] # catch nearest data
    bins = [] # find the edges of the bins
    for nt in neartime:
        ii = abs(xx-nt)<30
        if sum(ii): 
            nb = [min(xx[ii])-0.1, max(xx[ii])+0.1] # some padding
            # avoid duplicates
            if not nb in bins: 
                bins.append(nb)
    
    bins = np.array(bins).flatten()

    # compute: W1
    ii1 = (dd['w1sigmpro']>0) 

    print ('# W1 (significant) detections: ({0}) {1}'.format(sum(ii1), sum(dd['w1sigflux']>0)))
    xbin1, ybin1 = sjoert.simstat.binthem(xx[ii1], dd['w1mpro'][ii1],dd['w1sigmpro'][ii1], use_wmean=True,std=True, sqrtN=True, bins=bins, silent=True) 
    inz1 = ybin1[0,:]>0 

    if sum(ii1):
        w1_med = np.median(dd['w1mpro'][ii1])
        w1_N = sum(ybin1[3,:])
        #w1_chi2 = sum((dd['w1mpro'][ii1]-w1_med)**2 / dd['w1sigmpro'][ii1]**2) / sum(ii1)  # use all data with reported uncertainty
        inz1_chi=ybin1[3,:]>2 # need enough detection to compute scatter
        w1_chi2 = sum((ybin1[0,inz1_chi]-w1_med)**2 / ybin1[1,inz1_chi]**2) / sum(inz1_chi)          # use uncertainty estimated from scatter in each bin
    else:
        w1_N, w1_med, w1_chi2=0, np.nan, np.nan


    # plot and compute: W2
    
    ii2_mag = (dd['w2sigmpro']>0)   # significant detections
    ii2 = dd['w2sigflux']>0         # all detections
    print ('# W2 (significant) detections: ({0}) {1}'.format(sum(ii2_mag), sum(ii2)))
    if sum(ii2) and sum(ii2_mag):

        zp2 = np.median(dd['w2mpro'][ii2_mag]+2.5*np.log10(dd['w2flux'][ii2_mag]))
        if verbose:
            print ('mean zeropoint W2:', zp2)

        xbin2, ybin2 = sjoert.simstat.binthem(xx[ii2], dd['w2flux'][ii2],dd['w2sigflux'][ii2], use_wmean=True,std=True, sqrtN=True, bins=bins, silent=True)
        
        inz2 = ybin2[0,:]>0 # note, this removed epoch where the mean flux is negative 

        # convert mean flux to mag
        w2binflux = ybin2[0,:].copy()
        ybin2[0,:]= -2.5*np.log10(w2binflux)+zp2 
        ybin2[1,:] = 2.5/np.log(10) * ybin2[1,:]/w2binflux

        w2_med = -2.5*np.log10(sjoert.simstat.wmean(dd['w2flux'][ii2],dd['w2sigflux'][ii2]))+zp2 # convert to mag
        w2_N = sum(ybin2[3,:])
        
        #w2_chi2 = sum((dd['w2mpro'][ii2]-w2_med)**2 / dd['w2sigmpro'][ii2]**2) / sum(ii2) # use all data with reported uncertainty
        inz2_chi=ybin2[3,:]>2 # need enough detection to compute scatter
        w2_chi2 = sum((ybin2[0,inz2_chi]-w2_med)**2 / ybin2[1,inz2_chi]**2) / sum(inz2_chi)         # use uncertainty estimated from scatter in each bin
    else:
        w2_N, w2_med, w2_chi2  =0, np.nan, np.nan

    # compute W1-W2 in each bin
    w1minw2 = []
    if (sum(ii1)>10) and (sum(ii2)>10) and (sum(ii2_mag)>5):        
        for x in xbin2[inz2]:
            # check that we have data in the same bins
            i1_near = abs(xbin1[inz1]-x)<20
            i2_near = abs(xbin2[inz2]-x)<20
            n1_near = sum(ybin1[3,:][inz1][i1_near])
            n2_near = sum(ybin2[3,:][inz2][i2_near])
            if (n1_near>2) and (n2_near>2):
                w1_near = np.mean(ybin1[0,:][inz1][i1_near])
                w2_near = np.mean(ybin2[0,:][inz2][i2_near])
                w1minw2.append(w1_near-w2_near)
                if verbose:
                    print ('N(W1)={0:3.0f} N(W2)={1:3.0f} <W1>={2:0.2f} <W2>={3:0.3f} <W1>-<W2>={4:5.3f}'.format(n1_near,n2_near,w1_near,w2_near, w1minw2[-1]))

    # do classification
    wise_class = None

    # and make a nice string that we can return
    if w2_N==0 and w1_N==0:
        wise_info = "neoWISE: no match"
    else:
        wise_info = "neoWISE: N(w1)={0:0.0f}; N(w2)={1:0.0f}".format(w1_N, w2_N)

    if w1_N>0:
        wise_info += '; <w1>={0:0.2f}'.format(w1_med)
    if w2_N>0:
        wise_info += '; <w2>={0:0.2f}'.format(w2_med)

    if not np.isnan(np.median(w1minw2)):
        wise_info += '; <w1-w2>={0:0.2f}'.format(np.median(w1minw2))

    if not np.isnan(w1_chi2):
        wise_info += '; chi2(w1)={0:0.1f}'.format(w1_chi2)

    if not np.isnan(w2_chi2):
        wise_info += '; chi2(w2)={0:0.1f}'.format(w2_chi2)
 
    # we want at least 10 W2 detections over 3 epochs
    if w2_N>10:
        if sum(ybin2[3,inz2])>10 and (len(xbin2[inz2])>3) :
            
            # Stern12
            # if (w2_med<15) and (np.median(w1minw2)>0.8): 
            #   wise_class = 'AGN'
            
            # Assaf+13
            if np.median(w1minw2)>0.662*np.exp(0.232*(np.clip(w2_med-13.97,0,1e99))**2):
                wise_class = 'AGN'
                wise_info+='; AGN classification based on color (Assef+12)'

        # significant variability
        if w1_chi2>10:
            wise_class = 'AGN'
            wise_info+='; significant variability'
        # evidence for variability
        elif w1_chi2>5:
            wise_class = 'AGN?'
            wise_info+='; some evidence for variability'

    info_list = [wise_info]

    # also add allWISE info
    if max(dd['w2mpro_allwise']) and max(dd['w1mpro_allwise']):
        w1w2_allwise = max(dd['w1mpro_allwise'])-max(dd['w2mpro_allwise'])
        w1w2_allwise_err = np.sqrt(max(dd['w1sigmpro_allwise'])**2+max(dd['w2sigmpro_allwise'])**2)
        all_info = 'allWISE: w1-w2={0:0.2f}+/-{1:0.2f}'.format(w1w2_allwise,w1w2_allwise_err)

        if np.median(w1w2_allwise)>0.662*np.exp(0.232*(np.clip(dd[0]['w2mpro_allwise']-13.97,0,1e99))**2):
            wise_class = 'AGN'
            all_info+='; AGN classification based on color (Assef+12)'

        info_list += [all_info]
    


    if verbose:
        print ('allWISE: median(w1)-median(w2) = {0:0.3f}'.format(w1_med-w2_med))
        [print(x) for x in info_list]
        print ('neoWISE classification:', wise_class)


    # make a plot
    if not os.path.isfile(name+'-neoWISE.pdf') or redo or wait:

        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if w1_N:
            line = ax.errorbar(xx[ii1], dd['w1mpro'][ii1],dd['w1sigmpro'][ii1], fmt='s', alpha=0.5,label='W1')
            ax.plot(xbin1[inz1], ybin1[0,inz1],'d--', color=line[0].get_color())
            ax.fill_between(xbin1[inz1], ybin1[0,inz1]+ybin1[1,inz1], ybin1[0,inz1]-ybin1[1,inz1], color=line[0].get_color(), alpha=0.5)
            # if dd[0]['w1mpro_allwise']:
            #   ax.plot(xx, np.repeat(dd[0]['w1mpro_allwise'], len(xx)), ':', lw=2, color=line[0].get_color(), alpha=0.6)
        
            if max(dd['w1mpro_allwise']):
                xw = (min(xbin1[inz1])-10, max(xbin1[inz1]))
                plt.plot(xw, np.repeat(max(dd['w1mpro_allwise']), 2), 
                    ':',lw=1,alpha=0.8, color=line[0].get_color())

        if w2_N:
            line = ax.errorbar(xx[ii2], dd['w2mpro'][ii2],dd['w2sigmpro'][ii2], fmt='o', alpha=0.5,label='W2')
            ax.plot(xbin2[inz2], ybin2[0,inz2],'d--', color=line[0].get_color())
            ax.fill_between(xbin2[inz2], ybin2[0,inz2]+ybin2[1,inz2], ybin2[0,inz2]-ybin2[1,inz2], color=line[0].get_color(), alpha=0.5)
            # if dd[0]['w2mpro_allwise']:
            #   ax.plot(xx, np.repeat(dd[0]['w2mpro_allwise'], len(xx)), ':', lw=2, color=line[0].get_color(), alpha=0.6)

            if max(dd['w2mpro_allwise']):
                xw = (min(xbin2[inz2])-10, max(xbin1[inz2]))
                plt.plot(xw, np.repeat(max(dd['w2mpro_allwise']), 2),
                 ':', lw=1, alpha=0.8,color=line[0].get_color())



        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel('MJD - {0:0.0f}'.format(t0))
        ax.set_ylabel('mag (Vega)')
        ax.set_title("N(w1)={0:0.0f}; N(w2)={1:0.0f}; chi2(w1)={2:0.1f}; chi2(w2)={3:0.1f} <w1-w2>={4:0.2f}".format(
                                w1_N, w2_N, w1_chi2, w2_chi2, np.median(w1minw2)))
        ax.legend()

        if name:
            ax.figure.savefig(name+'-neoWISE.pdf') 
        if wait:
            ax.figure.show()
            key = input()

    return info_list, wise_class



def get_PS(ra, dec, name='', radius=1,
            t0=58119, wait=False, verbose=False, redo=False, 
            skip_detections=False, fluxtype='kron'):
    '''
    download and safe PS1 DR2 stack (and detections) 
    default search radius is 1" 

    >>> cats_dict, info_dict = get_PS(350.95256, -1.13618, radius=1, name='TDE2') 

    output:
    - cats_dict: is a dict that contain the catalogs for 'stack' and 'detection'
    - info_dict: summary info on the detections table, contains the chi2, median, etc. of the different filters 

    the summary info of the detections (and the plot) only makes sense if you have one source in the detections table

    input:
    - name: path for saving data and plot detections
    - radius: match radius in arsec, default is 1"
    - skip_detections=False: dont safe and plot the detection table
    - t0: mjd, plot time relative to this; default is the year 2018 (58119)
    - redo=False: force download even if we have the data on disk
    - fluxtype=kron: can also be set to psf
    - verbose=False: talk to me 
    - wait=False: pause for interactive plots

    '''

    flt_dict = {1:'g', 2:'r', 3:'i', 4:'z', 5:'y'}

    url_2fill = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/' + \
                '{0}.csv?flatten_response=false&raw=false&sort_by=distance&{1}&radius={2:0.6f}'

    fluxtype = fluxtype.lower() # lets be helpful

    # deal with inconsitent columnnames between PS1 cvs files
    if fluxtype=='psf':
        stack_fluxtype = 'PSF'
    else:
        stack_fluxtype = 'Kron'


    # read/save different PS1 catalog
    astro_tab = {}

    if skip_detections:
        catnames = ['stack']
    else:
        catnames = ['detection', 'stack']

    for catname in catnames:

        fname = name+'-PS_'+catname+'.csv'
        if name and os.path.isfile(fname):
            
            if verbose:
                print ('reading:', fname)
            
            astro_tab[catname] = astropy.io.ascii.read(fname, delimiter=',', format='fast_csv', fast_reader=True) # this is kinda slow?

        else:
            redo = True

        if redo or not(name):
            obs_str = 'ra={0:0.6f}&dec={1:0.6f}'.format(ra, dec)

            url = url_2fill.format(catname, obs_str, radius/3600)
            if verbose:
                print ('getting PS data from:')
                print(url)

            try:
                r = requests.get(url, timeout=120)

                table=r.content.decode()

                # urlindex=content.find('/workspace/')
                # urlindex2=content.find('.tbl')

                # table=requests.get('https://irsa.ipac.caltech.edu'+content[urlindex:urlindex2]+'.tbl')

            except Exception as e:
                print ('PS: no connection to IRSA; url was:\n', url)
                print (e)
                return None, None

            if table:
                astro_tab[catname] = astropy.io.ascii.read(table)
            else:
                
                if verbose:
                    print ('no match')
                return None, None

            # save the data
            if name:
                with open(fname, 'w') as the_file:
                    the_file.write(table)

    
    if skip_detections:
        return astro_tab, {}

    # short-hand
    dd = astro_tab['detection']

    if len(dd)==0:
        return astro_tab, {}


    # apply some quality cuts
    iqual = (dd['infoFlag2']>-1) # this can be made better...
    iqual = (dd[fluxtype+'Flux']/dd[fluxtype+'FluxErr']>1) 
    if fluxtype=='kron':
        iqual *= (dd[fluxtype+'Rad']>0) # this can be made better...

    if verbose:
        print ('# of observations', len(dd))
        print ('# data point removed by quality cuts', sum(iqual==False))
    
    if sum(iqual)<1:
        if verbose:
            print ('no enough good data in detections:', astro_tab)
        return astro_tab, None

    dd = dd[iqual]

    xx = dd['obsTime']-t0

    # compute for each band
    N = {} # number of detection
    mjd = {}
    mag = {}
    mag_err = {}
    med = {} # median
    smag = {} # from stack
    chi2 = {} # chi2/dof

    for i in [1,2,3,4,5]:
        
        flt = flt_dict[i]
        idx = (dd['filterID']==i) * (dd[fluxtype+'Flux']/dd[fluxtype+'FluxErr']>1)
        N[flt] = sum(idx)

        if len(astro_tab['stack'][flt+stack_fluxtype+'Mag'])>1:
            print ('warning, mulitple matches in Stack... band={0}'.format(flt))
            print ('distance (arcsec) is ',np.array(astro_tab['stack']['distance'])*3600)
            print ('mag is ', np.array(astro_tab['stack'][flt+stack_fluxtype+'Mag']))
            smag[flt] = float(astro_tab['stack'][flt+stack_fluxtype+'Mag'][0])
            #key = input()
        else:
            smag[flt] = float(astro_tab['stack'][flt+stack_fluxtype+'Mag'])

        if verbose:
            print (flt, '# detections: {0}'.format(sum(idx)))


        if sum(idx):

            mjd[flt] = dd[idx]['obsTime']
            mag[flt] = -2.5*np.log10(dd[idx][fluxtype+'Flux']) + 8.9 # Jy to AB MAG
            mag_err[flt] = 2.5/np.log(10) * dd[idx][fluxtype+'FluxErr']/dd[idx][fluxtype+'Flux']

            med[flt] = np.median(mag[flt])
            chi2[flt] = sum( (mag[flt]-med[flt])**2 / mag_err[flt]**2) / (sum(idx)-1)  # use all data with reported uncertainty
            #inz1_chi=ybin1[3,:]>2 # need enough detection to compute scatter
            #chi2[flt] = sum((ybin1[0,inz1_chi]-w1_med)**2 / ybin1[1,inz1_chi]**2) / sum(inz1_chi)          # use uncertainty estimated from scatter in each bin
    

        else:
            N[flt], med[flt], chi2[flt]=0, np.nan, np.nan


    # and make a nice string that we can return
    if N['g']==0 and N['r']==0:
        ps_info = "PS: no match"
    else:
        ps_info = "PS: N(g)={0:0.0f}; N(r)={1:0.0f}".format(N['g'], N['r'])

    if N['g']>0:
        ps_info += '; <g>={0:0.2f}'.format(med['g'])
    if N['r']>0:
        ps_info += '; <r>={0:0.2f}'.format(med['r'])

    if not np.isnan(chi2['g']):
        ps_info += '; chi2(g)={0:0.1f}'.format(chi2['g'])

    if not np.isnan(chi2['r']):
        ps_info += '; chi2(r)={0:0.1f}'.format(chi2['r'])


    # make a plot 
    if not os.path.isfile(name+'-PS_detections.pdf') or redo or wait:

        import matplotlib
        from matplotlib import pyplot as plt

        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for i in [1,2,3,4,5]:
            
            flt = flt_dict[i]
            #idx = (dd['filterID']==i) 

            if N[flt]:
                line = ax.errorbar(mjd[flt]-t0,mag[flt],mag_err[flt], fmt='s', alpha=0.9,label=flt)
                ax.plot( (min(mjd[flt]-t0), max(mjd[flt]-t0)),(med[flt], med[flt]), '--', color=line[0].get_color())
                
                if smag[flt]>0:
                    ax.plot( (min(mjd[flt]-t0), max(mjd[flt]-t0)),(smag[flt],smag[flt]), ':', color=line[0].get_color())
            else:
                ax.plot( (min(dd['obsTime']-t0), max(dd['obsTime']-t0)), (smag[flt],smag[flt]), ':', label=flt)
                #ax.fill_between(xbin1[inz1], ybin1[0,inz1]+ybin1[1,inz1], ybin1[0,inz1]-ybin1[1,inz1], color=line[0].get_color(), alpha=0.5)


        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel('MJD - {0:0.0f}'.format(t0))
        ax.set_ylabel(fluxtype+' mag')
        ax.set_title("chi2(g)={0:0.1f}; chi2(r)={1:0.1f}; chi2(i)={2:0.1f} ".format(chi2['g'], chi2['r'], chi2['i']))
        ax.legend()

        if name:
            ax.figure.savefig(name+'-PS_detections.pdf')
        if wait:
            ax.figure.show()
            key = input()

    return astro_tab, {'median':med, 'N':N, 'chi2':chi2, 'stack':smag}



def get_SDSS_simple(ra, dec, rad=1/60., dir='./', name='', silent=False):
    '''
    >> data = get_SDSS_simple(ra, dec, rad=1, dir='./' name='mydata')

    call the DR12 skyserver with single [ra, dec] and download cvs file
    input:
     ra, dec (deg)
    optional input:
      radius=1/60. (deg).
      dir='./' directory for fits file
      name name for file, if not give we use IAU name of input coords
      silent=False shut it.

    note, make 1e5 lines are returned, only basic imaging data (magnitudes)

    '''

    SDSSurl_radec = '''"http://skyserver.sdss.org/dr12/en/tools/search/x_results.aspx?searchtool=Radial&uband=&gband=&rband=&iband=&zband=&jband=&hband=&kband=&TaskName=Skyserver.Search.Radial&whichphotometry=optical&coordtype=equatorial&ra=__RA__&dec=__DEC__&radius=__RAD__&min_u=0&max_u=23&min_g=0&max_g=23&min_r=0&max_r=23&min_i=0&max_i=23&min_z=0&max_z=23&min_j=0&max_j=23&min_h=0&max_h=23&min_k=0&max_k=23&format=csv&TableName=&limit=1000"'''
    url = SDSSurl_radec.replace('__RA__',str(float(ra))).replace('__DEC__', str(float(dec))).replace('__RAD__',str(rad*60))
    
    if not(name):
        name = 'SDSS-'+iau_name(ra, dec)+'.csv'

    com = '''curl '''+url+''' > ''' +dir+name

    if not silent: print(com)

    os.system(com)
    
    # the (stupid?) format with two tables
    f = open(dir+name,'r')
    lines = ''
    for l in f.readlines():
        #print(l)
        if l.count('#Table2')==1: 
            break
        lines+=l

    lines= lines.split('\n')
    print(lines[1])
    print(lines[2:])
    data = readascii(lines=lines[2:], delimiter=',', names=lines[1].split(','))
    
    return data

# ----
# some stuff for NED, used for "ragolu" paper (these function are very old by now and may not work with the latest NED)


# change this variable to change the default input to all NED functions
# or give NEDdir explicitly to each call
NEDdir = dirs.catdir+'NED/'

NEDurl = '"http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=_SOURCE_&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=NO"'
NEDurl_radec = '''http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=_RA_d&lat=_DEC_d&radius=_RAD_&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=NO'''



def get_NED_name(name=None,ra=None, dec=None, rad=.1/60., NEDdir=NEDdir, redo=False):
    '''
    download or read NED page, return NED name of the galaxy
    >> get_NED_name(name, NEDdir=NEDdir)
    >> get_NED_name(ra=130,dec=-12, rad=0.1/60, NEDdir=NEDdir)
    '''
    
    if name is not None:
        local_name = name.strip().replace(' ', '_')
    else:
        local_name = iau_name(ra, dec)
    
    fname  = NEDdir+local_name+'.html'
    
    if not(os.path.isfile(fname)):
        print(fname)
        print('not on disk, sendig query to NED')
    
    if not(os.path.isfile(fname)) or redo:
        download_NED(name=name, ra=ra,dec=dec, rad=rad, local_name=local_name, NEDdir=NEDdir)

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
            if verbose: print(lines[i])
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

def download_NED(name=None, ra=None, dec=None, rad=None, local_name=None,
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
    if name is not None: 
        html_name = name.strip().replace(' ','_')
        command_line = "wget '"+NEDurl.replace('_SOURCE_', html_name)+ \
          "' -O "+NEDdir+local_name+'.html'
    else:
        if rad is None: 
            rad = .1/60. 

        command_line = "wget '"+NEDurl_radec.replace('_RA_', str(ra)[0:9]).replace('_DEC_', str(dec)[0:9]).replace('_RAD_', str(rad*60))+ \
          "' -O "+NEDdir+local_name+'.html'

    print(command_line)
    args = shlex.split(command_line)

    if debug:
        print('pasing to command line  :\n', command_line,' \n\n')
        print( args)
        print('NEDdir:', NEDdir)
        
    #proc = subprocess.Popen(args cwd=NEDdir, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    p = subprocess.call(command_line, shell=True)
    

    # read the url to the SED page
    page = open(NEDdir+local_name+'.html','r')
    lines = page.readlines()
    i = 0


    if not SEDpage_name: return

    while i<len(lines):

        if debug: print(lines[i])

        # get the link to the SED page
        if (re.search('Photometric', lines[i])):

            link = '"http://ned.ipac.caltech.edu/'+ \
                   lines[i].split('REF=')[1].split('TARGET')[0].strip()+'"'

            com = ' wget '+link+' -O  '+NEDdir + SEDpage_name+'.html'
            print ('pasing to command line  :\n', com, ' \n\n')
            os.system(com)
            return

        i+=1
        
        
    print('WARNING No link to SED page found on NED main page.')


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
            # check if not an upper limit
            if (re.search('=', lflux)):                   
                lflux = lflux.split('=')[1].split('W')[0].strip()
                flux.append(float(lflux) *1e23/1e4*1e7)
                i+=3
                lnu = lines[i]
                lnu = lnu.split('=')[1].split('H')[0].strip()
                nu.append(float(lnu))
            
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
        sdir = dirs.catdir+'DSS/fits/'
        indir = dirs.catdir+'DSS/infiles/' #keep log somewhere else
    else:
        indir=sdir

    new_name = sdir+name+'.fits'
    if noclobber and os.path.isfile(new_name):
        print('not rerunning, reading:', new_name)
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
        print(sradec)
        print(line)
        print(command_line)

    #p = subprocess.Popen(command_line, shell=True)
    #t = os.popen(command_line)
    #p = subprocess.call(command_line, shell=True)

    os.system('touch '+logfile_name)
    logfile = open(logfile_name, 'r')

    args = shlex.split(command_line)
    if debug: print(args)
    proc = subprocess.Popen(args, cwd=sdir, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    sleep_time = (size*60)
    if sleep_time <10: sleep_time = 10
    print ('sleeping for (s):', sleep_time)
    time.sleep(sleep_time )
    # at this point we either have an image, or we give up
    proc.terminate()
    
    if proc.wait() < 0:
      
        print ('problem with dss call')
        print ('id:', proc.pid)
        #proc.kill()
        #os.kill(proc.pid,  signal.SIGHUP) # doesn't work ?
        return None
    #time.sleep(1)
    #proc.kill()


    #stdout_value = proc.communicate()[0]
    stderr_value = proc.stderr.read()
    print ('stderr:', stderr_value)
    stdout_value = proc.stdout.read()


    if debug:
        print ('\n stdout: ', stdout_value)

    lines = (stdout_value).split('\n')

    if not(len(lines)):
        print ('empty log file', logfile_name)
        return None

    i=0
    while not re.search("Field", lines[i]) and (i < (len(lines)-1)):
        i+=1

    if i == (len(lines)-3):
        print ('no info from log file:',logfile_name)
        return None

    info_line = lines[i+2].split()
    if len(info_line) < 7:
        print ('no info on plate id:', info_line)
        return None
    
    Pl_Id = lines[i+2].split()[8].swapcase()
    if debug: print ('Pl_Id:', Pl_Id)
    force_name = sdir+name[0:20]+Pl_Id+'.fits' # max 20 characters...
    if debug: print ('dss image name:', force_name)
    
    if change_name:
        os.system('mv '+force_name+' '+new_name)
        if debug: print ('new name:', new_name)

        try:
            hdu = pyfits.open(new_name)
            print ('created:', new_name)
        except IOError:
            print ('file not created:', new_name)
            return None

    else: 
        try:
            hdu = pyfits.open(force_name)
            print ('created:', force_name)
        except IOError:
            print ('file not created:', force_name)
            return None
        
    return hdu
