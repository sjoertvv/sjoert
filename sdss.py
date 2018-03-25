'''
countainer with SDSS related functions
'''

import numpy as np

bvalues = np.array([1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10])
bvalue_82 = np.array([1.0e-11, 0.43e-11, 0.81e-11, 1.4e-11, 3.7e-11])

def lup_to_maggie(lup, lup_sig=[], do_82=False):

    '''
    from https://github.com/ubutsu/kcorrect
    - if lup_sig is set, we return the sigma (not ivar) 
    - if do_82, using Stripe 82 softening parameters
    '''

    b = bvalues
    if do_82: 
        b = bvalue_82
    maggie = 2. * b * np.sinh(-np.log(b) - 0.4 * np.log(10.) * lup)
    if len(lup_sig) != 0:
        
        print ('uncertainties given,; omputing the ivar of maggies')
        if lup.shape != lup_sig.shape:
            raise ValueError('Array size different for lup and lup_sig.')
        maggie_sigma = (2. * b * np.cosh(-np.log(b) - 0.4 * np.log(10.) * lup)
                        * 0.4 * np.log(10.) * lup_sig)

        return maggie, maggie_sigma
    return maggie



sdss_l = np.array([3543, 4770, 6231, 7625, 9134]) # SDSS filter in A

def sdss_nu(k):
    return 3e8 / (sdss_l[k]*1e-10)

def sdss_lambda(k):
    return sdss_l[k]

def get_nu(band):
    nu = None
    if band == 'FUV': nu = 3e8 / (1528 *1e-10)
    if band == 'NUV': nu = 3e8/ (2271 *1e-10)
    if band == 'u': nu = sdss_nu(0)
    if band == 'g': nu = sdss_nu(1)
    if band == 'r': nu = sdss_nu(2)
    if band == 'i': nu = sdss_nu(3)
    if band == 'z': nu = sdss_nu(4)
    if not(nu):
        print ('please use band=[FUV, ..., z]')
    return nu