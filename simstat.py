'''
simple/basic statistical tools mostly about counting things (ie, Poisson)
2010 - Sjoert van Velzen
'''
import numpy as np
from numpy.random import rand
from matplotlib import pyplot as plt
import scipy.stats.distributions as pydist
import time

def Gauss(x, mu=0, sigma=1):
    '''
    Gaussian probablity distribution
    P = Gauss(x, mu=0, sigma=1)
    '''
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

def poisson_limits(n_i, conf):
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
        print 'poisson_limits:: n_i==0, returning [0,2.4]'
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


def his2scat(arr, bins=10, range=None, conv=0.8413, logbin=False, silent=False):
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

    xx = np.zeros(len(ipos))
    yy = np.zeros(len(ipos))
    yy_err = np.zeros((2,len(ipos)))
    if len(ipos) != len(hh[0]):
        if not(silent):
            print 'warning we have empty bins! skipping these'

    for i, ip in enumerate(ipos):

        xx[i] = (hh[1][ip]+hh[1][ip+1])/2.

        bin_width =  hh[1][ip+1]-hh[1][ip]
        if logbin:
            bin_width = 10**(hh[1][ip+1])-10**(hh[1][ip])
            
        yy[i] = hh[0][ip] / bin_width
        plim = poisson_limits( hh[0][ip], conv)
        yy_err[:,i] =  abs(plim - hh[0][ip]) / bin_width
        
    return xx, yy, yy_err


def binthem(x, y, bins=10,
            cl=0.9, poisson=False, std=False, sqrtN=True, use_mean=False,
            silent=False):
    '''
    bin x,y using bins in x-direction
    >> xmid, ymid = binthem(x, y)

    output:
     xmid  (N) the mean value of x for each bin
     ymid  (3,N) the median/mean value, uncertainty/dispersion 

    the uncertainty is computed from the values in bin,
    or using std, poisson stat if flags are set 

    input:
     - x,y  equal length arrays
     - bins  number of xbins or array with bins (default is 10)
     - cl=0.9  confidence level for computing the uncertainty
     - poisson=False use Poisson stat to find uncertainty on number in bin
     - std=False  use standard deviation / sqrt(N) to compute uncertainty (ie, symtric errorbars)
     - sqrtN=True  set to False to use std only
     - use_mean=False  use mean (median is default) 
     - silent  shutup
    '''

    if np.isscalar(bins):
        x_bins = np.linspace(np.min(x), np.max(x)*1.01, bins)
    else:
        x_bins = bins

    xmid = np.zeros(len(x_bins)-1)
    ymid = np.zeros((3,len(x_bins)-1))

    for i in range(len(xmid)):

        ibin = np.where((x>=x_bins[i]) & (x<x_bins[i+1]))[0]

        xmid[i] = np.mean(x[ibin])

        if len(ibin)>0:
    
            y_arr = y[ibin]

            if use_mean:
                ymid[0,i] = np.mean(y_arr)
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
                    ymid[1,i] = np.abs(ymid[0,i]-np.interp( (cl/2.+0.5)*len(ibin), ii, y_arrs))
                    ymid[2,i] = np.abs(ymid[0,i]-np.interp( (cl/2.-0.5)*len(ibin), ii, y_arrs))

        if not silent:
            print xmid[i], len(ibin), ymid[:,i]

    return xmid, ymid

def poisson_limits_frac(numer=2, denom=4, cl=0.9,contained=True,
                        nmc=1000, verbose=False):
    '''
    >> low, up = poisson_limits_frac(numer=2, denom=4, cl=0.9, verbose=True)
    return CL interval from Possion probablity for observed fraction 
    if both numerator denominator are drawn from Poisson distribution
    option input:
    - contained=True, the numerator is *contained* in the fraction (ie, ratio<=1 is enforced)
    - verbose=False, make plots
    - nmc=1000 number of Monte Carlo samples

    (never used this function for any publication, it may actualy be wrong)
    '''

    nn = 1000
    time0= time.time()
    possible_numer = np.linspace(max(numer-8*np.sqrt(numer),numer/100.),
                                 numer+8*np.sqrt(numer),nn) # possible underlying true value
    possible_denom = np.linspace(max(denom-8*np.sqrt(denom),denom/100.),
                                 denom+8*np.sqrt(denom),nn) # possible underlying true value

    cvis_numer = 1-pydist.poisson.cdf(numer, possible_numer)
    cvis_denom = 1-pydist.poisson.cdf(denom, possible_denom)

    if verbose:
        print 'time for direct prob', time.time()-time0
        plt.clf()
        vis = pydist.poisson.pmf(numer, possible_numer)
        vis2 = pydist.poisson.pmf(denom, possible_denom)
        plt.plot(possible_numer, vis, 'o-')
        plt.plot(possible_denom, vis2, 'x-')
        plt.plot(possible_numer, cvis_numer, 'x-')

        key = raw_input()

    # draw possible true means from poisson distribution
    time0= time.time()
    uni = np.random.rand(2,nmc)
    ratio = np.zeros(nmc)
    temp = np.zeros(nmc)

    for i in range(nmc):
        pn = np.interp(uni[0,i], cvis_numer, possible_numer)
        pd = np.interp(uni[1,i], cvis_denom, possible_denom)
        if verbose: print uni[0,i], uni[1,i],  pn, pd
        ratio[i] = pn/pd
        temp[i] = pn
    
    if contained: ratio = ratio[np.where(ratio<=1)]

    # make cumulative distribution
    ratio.sort()
    cum = np.arange(len(ratio))/float(len(ratio))
    temp = np.sort(temp)
    cum_temp = np.arange(len(temp))/float(len(temp))

    if verbose:
        hh = plt.hist(temp, bins=200, normed=True, alpha=0.5)
        xx=np.zeros(len(hh[0]))
        yy = np.zeros(len(hh[0]))
        for i in np.arange(len(hh[0])-1)+1:
            ww = hh[1][i+1]-hh[1][i]
            xx[i] = hh[1][i]+ww*0.5
            yy[i] = yy[i-1]+hh[0][i]*ww
        plt.plot(xx,yy)
        key = raw_input(' next, hist')

        plt.hist(temp, bins=100, normed=True, cumulative=True, alpha=0.5)
        print 'time for mc of ratio (s)', time.time()-time0
        print 'mean, median ratio', np.mean(ratio), np.median(ratio)
        print 'p=0.5', np.interp(0.5, cum, ratio)
        print 'p=0.05',np.interp(0.05, cum, ratio) 
        print 'p=0.95', np.interp(0.95, cum, ratio)
        print '90% CL only on numinator', poisson_limits(numer, 0.9)/denom

        print ''
        print 'p=0.05, numer',np.interp(0.05, cum_temp, temp), temp[0.05*nmc]
        print 'p=0.05, numer after his', np.interp(0.05,yy,xx)
        print 'p=0.975, numer', np.interp(0.975, cum_temp, temp),temp[0.975*nmc]
        print 'p=0.01, numer',np.interp(0.01, cum_temp, temp)
        print 'p=0.99, numer', np.interp(0.99, cum_temp, temp)
        print '90% CL  numinator', poisson_limits(numer, 0.9)
        print '95% CL  numinator', poisson_limits(numer, 0.95)
        print '99% CL  numinator', poisson_limits(numer, 0.99)


        key = raw_input(' next, ratio hist')
        plt.hist(ratio, bins=30, alpha=0.5, normed=True)
        plt.plot(ratio,cum)
        key = raw_input()

    
    return np.interp(0.5-cl/2., cum, ratio), np.interp(.5+cl/2., cum, ratio)
