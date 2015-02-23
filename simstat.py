'''
simple/basic statistical tools mostly about counting things (ie, Poisson)
2010 - Sjoert van Velzen
'''
import numpy as np
from numpy.random import rand
from matplotlib import pyplot as plt
import scipy.stats.distributions as pydist
import time

def Gauss(x, mu=0, sigma=1, fwhm=None):
    '''
    Gaussian probablity distribution
    P = Gauss(x, mu=0, sigma=1)
    '''
    if fwhm is not None: 
        sigma = 0.5 * fwhm / np.sqrt(2*np.log(2))
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5 * (x-mu)**2 / (sigma**2) )

def wmean(d, ivar):
    '''
    >> mean = wmean(d, ivar)
    inverse-variance weighted mean
    '''
    d = np.array(d)
    ivar = np.array(ivar)
    return np.sum(d * ivar) / np.sum(ivar)

def mc_dice(array, nsucces=1):
    '''
    >> P = mc_dice(array, nsucces)
    find the probability to *at least* nsucces from array for probabilities to succes
    '''
    ntry = 5000*(max(len(array),2))
    #print 'drawing', len(array)*ntry, 'randon numbers'
    uni = rand(ntry*len(array))

    outcome = np.zeros(ntry)
    l=0
    
    for i in range(ntry):

        ns = 0
        
        for j in range(len(array)):
            if uni[l]<array[j]:
                ns+=1
            l+=1

        if ns>=nsucces:
            outcome[i] = 1
        
    return np.sum(outcome) / float(ntry)
        
    
def poisson(actual, mean):
    '''
    >> p = poisson(actual, mean)
    give Poisson probablity to find <actual> (int) if true value is <mean>
    uses pydist.poisson.pmf
    '''
    return pydist.poisson.pmf(actual,mean)

def poisson_limits(n_i, conf, silent=False):
    '''
    NAME:
    IM_POISSON_LIMITS()
  
    PURPOSE:
    function im_poisson_limits, n_i, conf
    For a given number of events, calculate the Poisson upper and
    lower limits for a given 1-sided confidence interval. 
 
    INPUTS: 
    n_i  - number of (counting) events or samples (must be a scalar) 
    conf - desired confidence interval for example, 0.8413, 0.9772,
    and 0.9987 correspond to the usual 1-sided intervals for 1-, 2-, 
    and 3-sigma Gaussian probabilities (must be a scalar)
 
    OUTPUTS: 
    limit - two-element array corresponding to the lower- and
    upper- confidence limit
    
    OPTIONAL OUTPUTS:
 
    COMMENTS:
    Formulae taken from Gehrels (1986).
 
    MODIFICATION HISTORY:
    J. Moustakas, 2008 Oct 17, NYU - written based on a Perl script
    kindly provided by G. Rudnick 
    S. van Velzen, 2012 Jan 18, Radboud - converted to python
    '''
    if not(np.isscalar(n_i)) or not(np.isscalar(conf)):
        print 'N_I and CONF must be scalars'
        return -1.0
    
    if (n_i < 0.0) :
        print 'N_I must be positive!'
        return -1.0

    if n_i == 0:
        if not(silent): print 'poisson_limits:: n_i==0, returning [0,2.4]'
        return np.array([0.,2.4])
    
    
# define a grid of confidence intervals
    conf_grid = np.array([0.8413, 0.90, 0.95, 0.975, 0.9772, 0.990, 0.995, 0.9987, 0.999, 0.9995])
    if (conf < min(conf_grid)) or (conf > max(conf_grid)) :
        print 'CONF must be in the interval [',min(conf_grid),'-',max(conf_grid),']'
        return -1.0
    

    # 'S' parameter corresponding to CONF_GRID
    SS = np.array([1.0, 1.282, 1.645, 1.960, 2.00, 2.326, 2.576, 3.000, 3.090, 3.291])

    # additional parameters needed to compute the lower limit
    beta = np.array([0.0, 0.010, 0.031, 0.058, 0.062, 0.103, 0.141, 0.222, 0.241, 0.287])
    gamma = np.array([0.0, -4.0, -2.50, -2.22, -2.19, -2.07, -2.0, -1.88, -1.85, -1.80])

    # eqs. (10) and (14) from Geherls (1986)    
    limit_up_grid = n_i + SS*np.sqrt(n_i+1.0) + (SS**2.0+2.0)/3.0
    if (n_i < 1.0):
        limit_dn_grid = 0.0
    else:
        limit_dn_grid = n_i * (1.0 - 1.0/(9.0*n_i) - SS/(3.0*np.sqrt(n_i)) + beta*n_i**(gamma))**3.0

    limit_up = np.interp(conf, conf_grid, limit_up_grid)
    limit_dn = np.interp(conf, conf_grid, limit_dn_grid)
    
    return np.array([limit_dn,limit_up])


def his2scat(arr, bins=10, range=None, conv=0.8413, logbin=False, return_zero=False, silent=False):
    '''
    make scatter points form array (eg, for luminosity function)
    >> xx, yy, yy_err = his2scat(arr, bins=10, range=None, conv=0.8413)

    output is: the mid of the bins, the normalized value, the uncertainty for conv interval
 
    normalization is such that integration over the bins yields the length of the array
    use logbin flag to if input in log10, but normalization needs to be linear
    (ie, input is logS, but output dN/dS) 
    
    note, negative bins are ignored
    todo: return upper limits for conv
    '''

#  make python histrogram
    isfin = np.where(np.isfinite(arr)==True)[0]
    if len(isfin)==0:
        print 'all bins are inf/nan?'
        return None
    hh = np.histogram(arr[isfin], bins=bins, range=range, normed=False)

    if not(silent):
        print 'number in bins', hh[0]
    ipos= np.where(hh[0]>0)[0]

    if (len(ipos) != len(hh[0])) and (return_zero==False):
        if not(silent):
            print 'warning we have empty bins! skipping these'
    else:           
        ipos = np.arange(len(hh[0]))
    xx = np.zeros(len(ipos))
    yy = np.zeros(len(ipos))
    yy_err = np.zeros((2,len(ipos)))



    for i, ip in enumerate(ipos):

        xx[i] = (hh[1][ip]+hh[1][ip+1])/2.

        bin_width =  hh[1][ip+1]-hh[1][ip]
        if logbin:
            bin_width = 10**(hh[1][ip+1])-10**(hh[1][ip])
            
        yy[i] = hh[0][ip] / bin_width
        plim = poisson_limits( hh[0][ip], conv)
        yy_err[:,i] =  abs(plim - hh[0][ip]) / bin_width
        
    return xx, yy, yy_err


def binthem(x, y, bins=10, range=[], use_mean=False,use_sum=False,
            cl=0.9, poisson=False, std=False, sqrtN=True, 
            silent=False):
    '''
    >> xmid, ymid = binthem(x, y)

    bin parameter y using bins in x-direction

    output:
     xmid  (N) the mean value of x for each bin
     ymid  (3,N) the median/mean/total value, uncertainty/dispersion 

    the uncertainty is computed from the values in bin,
    or using std, poisson stat if flags are set 

    input:
     - x,y  equal length arrays
     - bins  number of xbins or array with bins (default length is 10)
     - range=[xmin, xmax] range for the bins (default is [min(x),max(x)])
     - cl=0.9  confidence level for computing the uncertainty
     - poisson=False use Poisson stat to find uncertainty on number in bin
     - std=False  use standard deviation / sqrt(N) to compute uncertainty (ie, symtric errorbars)
     - sqrtN=True  set to False to use std only
     - use_mean=False  use mean (median is default) 
     - use_sum=False sum bin content
     - silent  shutup
    '''

    if np.isscalar(bins):
        if len(range)==2:
            x_bins = np.linspace(range[0], range[1], bins)
        else:
            x_bins = np.linspace(np.min(x), np.max(x), bins) #removed in 2015: max()*1.01
    else:
        x_bins = bins

    xmid = np.zeros(len(x_bins)-1)
    ymid = np.zeros((3,len(x_bins)-1))

    for i in np.arange(len(xmid)):

        ibin = np.where((x>=x_bins[i]) & (x<x_bins[i+1]))[0]

        xmid[i] = np.mean(x[ibin])

        if len(ibin)>0:
    
            y_arr = y[ibin]

            if use_mean:
                ymid[0,i] = np.mean(y_arr)
            elif use_sum:
                ymid[0,i] = np.sum(y_arr)
            else:
                ymid[0,i] = np.median(y_arr)

            if len(ibin)>1:
                if std:
                    ymid[[1,2],i] = np.std(y_arr)
                    if sqrtN:
                        ymid[[1,2],i]= ymid[[1,2],i]/ np.sqrt(len(ibin))
                elif poisson:
                    ymid[[1,2],i] = np.abs(poisson_limits(np.sum(y[ibin]), cl)/len(ibin)- ymid[0,i])
                else:
                    y_arrs = np.sort(y_arr)
                    ii = np.arange(len(ibin))
                    ymid[1,i] = np.abs(ymid[0,i]-np.interp( (cl/2.-0.5)*len(ibin), ii, y_arrs))
                    ymid[2,i] = np.abs(ymid[0,i]-np.interp( (cl/2.+0.5)*len(ibin), ii, y_arrs))


        if not silent:
            print '{0:0.2f} - {1:0.2f}  {2:0.0f}  [{3:0.2f}  {4:0.2f}  {5:0.2f}]'.format(x_bins[i],x_bins[i+1], len(ibin), ymid[0,i], ymid[1,i], ymid[2,i])

    return xmid, ymid


def cdf_match(x,y, new_len=None):
    '''
    >>> y_sampled, idx = cdf_match(x,y, oversample=False)
    
    Retrun subset of y, plus idx to this subset. 
    The subset is constructed to have the same 
    cummulative distribution function as x.
    Note that for default operation, we do not
    oversample. So len(y)>>len(x) is assumed. 
    
    optional input: 
    - new_len=N, make length of y equal to N 
    - oversample=False, allow multiple entries of y in y_sampled. 
    ''' 

    if new_len is None:
        new_len = len(x) 

    # drawn new y 
    x_cdf = np.sort(x), np.linspace(0, 1, len(x))
    y_pref =  np.interp(np.random.rand(new_len), x_cdf[1], x_cdf[0])

    # get index to true y
    #y_idx = np.round(np.interp(y_pref, np.sort(y), np.arange(len(y))))
    #y_idx = np.array(y_idx, dtype=np.int)
    #y_out = y[y_idx]
    
    y_out = np.zeros(new_len)
    y_idx = np.zeros(new_len, dtype=np.int)
    y_crop = y
    
    for i in range(new_len): 
        idx = np.int(np.floor(np.interp(y_pref[i], np.sort(y_crop), np.arange(len(y_crop)))))
        y_out[i] = y_crop[np.argsort(y_crop)[idx]]
        if not(oversample):
            y_crop = np.delete(y_crop, idx)

    for i in range(new_len): 
        i_ori  = np.where(y==y_out[i])[0]
        y_idx[i] = i_ori[np.int(np.random.rand()*len(i_ori))] # account for non-unique input or output

        
    return  y_out, y_idx
    
                     
                              
    
    
    


