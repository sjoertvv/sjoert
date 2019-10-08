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
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp( -0.5 * (x-mu)**2 / (sigma**2) )

def sammy(y,x, n=int(1e6)):
    '''
    returns x samples for a distribution of y(x) 
    x has to be monotonic increasing/decreasing, but doesnt have to be equally spaced
    '''
    y_int = np.zeros(len(x))
    x_int = np.zeros(len(x))
    x_int[0]=x[0]
    for i in range(1,len(y)):
        y_int[i] = np.trapz(y[0:i], x[0:i])
        x_int[i] = x[i]-(x[i]-x[i-1])*0.5
    y_int/=y_int[-1]
    
    return np.interp(np.random.rand(int(n)), y_int, x_int)

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
        print('N_I and CONF must be scalars')
        return None
    
    if (n_i < 0.0) :
        print('N_I must be positive')
        return None

    if n_i == 0:
        if not(silent): print('poisson_limits:: n_i==0, returning [0,2.4]')
        return np.array([0.,2.4])
    
    
# define a grid of confidence intervals
    conf_grid = np.array([0.8413, 0.90, 0.95, 0.975, 0.9772, 0.990, 0.995, 0.9987, 0.999, 0.9995])
    if (conf < min(conf_grid)) or (conf > max(conf_grid)) :
        print ('CONF must be in the interval [',min(conf_grid),'-',max(conf_grid),']')
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
        print('all bins are inf/nan?')
        return None
    hh = np.histogram(arr[isfin], bins=bins, range=range, normed=False)

    if not(silent):
        print('number in bins', hh[0])
    ipos= np.where(hh[0]>0)[0]

    if (len(ipos) != len(hh[0])) and (return_zero==False):
        if not(silent):
            print('warning we have empty bins! skipping these')
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


def binthem(x, y, yerr=None, bins=10, range=[], 
            use_mean=False,use_sum=False, use_wmean=False,
            cl=0.9, poisson=False, std=False, sqrtN=True, 
            smart_expand=True,
            silent=False):
    '''
    >> xmid, ymid = binthem(x, y)

    bin parameter y using bins in x-direction

    output:  
     - xmid  (N) the mean value of x for each bin  
     - ymid  (4,N) containing the [median/mean/total value, -/+ of uncertainty/dispersion, number in bin]

    input:
     - x,y  equal length arrays
     - bins  number of xbins or array with bins (default is 10)
    optional input:
     - range=[xmin, xmax] range for the bins (default is [min(x),max(x)])
     - cl=0.9  confidence level for computing the uncertainty
     - poisson=False use Poisson stat to find uncertainty on number in bin
     - std=False  use standard deviation / sqrt(N) to compute uncertainty (ie, symtric errorbars)
     - sqrtN=True  set to False to use std only
     - use_mean=False  use mean (median is default) 
     - use_sum=False sum bin content
     - use_wmean=False  compute weighted mean, requires yerr= input
     - smart_expand=False, if set, we try to increase bins to catch clusters of data (warning: experimentenal!)
     - silent  shutup
    '''

    if use_wmean:
        if yerr is None: 
            print('ERROR, please give yerr= input when wmean=True')
            return

    if np.isscalar(bins):
        if len(range)==2:
            x_bins = np.linspace(range[0], range[1], bins)
        else:
            x_bins = np.linspace(np.min(x), np.max(x), bins) 
    else:
        x_bins = bins.copy()

    if smart_expand:
        
        for i in np.arange(1,len(x_bins)):
            dbin = x_bins[i]-x_bins[i-1]
            iright = x>x_bins[i]            
            if sum(iright):
                # points closer to the current bins than next bins
                if min(x[iright]-x_bins[i])<dbin/2.:
                    if not silent:
                        print ('adjusting this bin', x_bins[i], 'by ', dbin/2.)
                        print ('new points added:', x[iright][x[iright]-x_bins[i]<dbin/2.])
                    x_bins[i] = x_bins[i] + dbin/2.

    xmid = np.zeros(len(x_bins)-1)
    ymid = np.zeros((5,len(x_bins)-1))

    for i in np.arange(len(xmid)):

        ibin = np.where((x>=x_bins[i]) & (x<x_bins[i+1]))[0]
        if i == (len(xmid)-1): # close last bin
            ibin = np.where((x>=x_bins[i]) & (x<=x_bins[i+1]))[0]

        if len(ibin)>0:

            xmid[i] = np.mean(x[ibin])
            ymid[4,i] = np.std(x[ibin])
    
            y_arr = y[ibin]
            ymid[3,i] = len(ibin)

            if use_mean:
                ymid[0,i] = np.mean(y_arr)
            elif use_sum:
                ymid[0,i] = np.sum(y_arr)
            elif use_wmean:
                y_arr_err = yerr[ibin]
                ymid[0,i] = wmean(y_arr, 1/y_arr_err*2)
            else:
                ymid[0,i] = np.median(y_arr)

            if std:
                ymid[[1,2],i] = np.std(y_arr)
                if sqrtN:
                    ymid[[1,2],i] = ymid[[1,2],i]/ np.sqrt(len(ibin))
            elif poisson and (len(ibin)>1):
                ymid[[1,2],i] = np.abs(poisson_limits(np.sum(y[ibin]), cl)/len(ibin)- ymid[0,i])
            elif use_wmean:
                ymid[[1,2],i] = 1/np.sqrt(sum(1/y_arr_err**2))
            else:
                y_arrs = np.sort(y_arr)
                ii = np.arange(len(ibin))
                ymid[1,i] = np.abs(ymid[0,i]-np.interp( (cl/2.-0.5)*len(ibin), ii, y_arrs))
                ymid[2,i] = np.abs(ymid[0,i]-np.interp( (cl/2.+0.5)*len(ibin), ii, y_arrs))
        else:
            xmid[i] = (x_bins[i] +x_bins[i+1])/2.

        if not silent:
            print('{0:0.2f} - {1:0.2f} ({2:0.2f})  {3:0.0f}  [{4:0.2f}  {5:0.2f}  {6:0.2f}]'.format(x_bins[i],x_bins[i+1], np.std(x[ibin]), ymid[3,i], ymid[0,i], ymid[1,i], ymid[2,i]))

    if sum(ymid[3,:]) != len(x):
        print ('binthem: WARNING: more points in bins ({0}) compared to lenght of input ({1}), please check your bins'.format(sum(ymid[3,:]), len(x)))
        key = input()

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
    - new_len=N        make length of y equal to N 
    - oversample=False if True, allow multiple entries of y in y_sampled. 
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
    
                     
def mad(arr):
    '''
        Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    '''
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))
                              
    
def unbinned_gauss(arr, sigma_cut=10, verbose=False):
    '''
    >> mean, sigma = unbinned_gauss([1,2,3])

    Fit normal distribution to array of data, with some robustness cooked in.

    As name suggest, no binning, we maximize the product of likelihoods
    of each element of the array (no weighting used). We guess sigma from 
    the MAD and mean from the median.

    optional input:
    - sigma_cut=10 clip the probability that this sigma (ie, less weight for outliers)

    '''

    mean_guess = np.median(arr)
    sigma_guess = mad(arr)*1.4826
    
    n = 100 # feels nice
    mean_arr = np.linspace(mean_guess-sigma_guess/1.5, mean_guess+sigma_guess/1.5, n) 
    sigma_arr = np.linspace(sigma_guess/1.5, sigma_guess*1.5, n)

    
    likel = np.zeros((n,n))
    for j,s in enumerate(sigma_arr):
        this_max = np.repeat(np.log(Gauss(sigma_cut*s, mu=0, sigma=s)), len(arr))        
        for i,m in enumerate(mean_arr):            
            likel[i,j] = np.sum(np.max([this_max, np.log(Gauss(arr, mu=m, sigma=s))],axis=0))
            #likel[i,j] = np.sum(np.log(Gauss(arr, mu=m, sigma=s)))
            #print(i,j,m,s,likel[i,j],sum(np.log(Gauss(arr, mu=m, sigma=s))<this_max[0]))

    mxl = np.where(likel == np.max(likel))
    
    if verbose: 
        print(mxl[0],mxl[1])
    
    if len(mxl[0])>1:
        print('unbinned_gaus, WARNING: more than one max likelihood, please check input')
        print(mxl)

    elif (mxl[0]==0) or (mxl[1]==0) or \
        (mxl[0]==0) or (mxl[1]==0):
        print('unbinned_gaus, WARNING: we hit an edge of sigma or mean array, please check input')
        print(mxl)

    return mean_arr[mxl[0]][0], sigma_arr[mxl[1]][0]

def test_unbinned(): 
    for x in np.linspace(0.8,3, 10):
        tt = np.append(np.random.normal(size=10**x), np.random.normal(-1, 5, size=10**x/5))
        m, s = unbinned_gauss(tt, sigma_cut=3)

        print('N true gauss',int(10**x), '  N false gauss', int(10**x/8))
        print('mean, median, likel  :',np.mean(tt), np.median(tt), m)
        print('std, mad*1.483, likel:', np.std(tt), mad(tt)*1.4826, s)
        print('')



