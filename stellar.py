'''
various function for astronomy. 
 - wrappers around cosmolopy (lumdis(z), dismod(z) etc.)
 - very simle friends-of-friends (one_group)
 - luminosity functions (schechter)
 - angular distance, and more
2010 May - started (Sjoert van Velzen, Radboud, NL)
2012 Aug - updated docstrings, import checks (SVV)
'''

try:
    from k3match import celestial as spherematch
except ImportError:
    print 'k3match import failed'
    print 'the fof function (one_group) wont work'
    print '(you can obtain k3match from pschella.github.com/k3match/)'

try:
    import cosmolopy as cospy
    sjoert_cosmo = cospy.parameters.WMAP7_BAO_H0_mean(flat=True, extras=False)
except ImportError:
    print 'sjoert.stellar ImportError: cosmolopy import failed, function using distances will fail'
    print '(get it from roban.github.com/CosmoloPy/)'


# import some useful stuff from pyspherematch's util.starutil_numpy
from starutil_numpy import radectolb, lbtoradec, ra2hmsstring, dec2dmsstring, hmsstring2ra, dmsstring2dec
import starutil_numpy as sutil


from numpy import floor, double
import numpy as np

# some constants
parsec = 3.08568e18 # cm
parsec_in_cm = parsec  
year_in_sec = 31557600.0

Mpc2cm = parsec_in_cm *1e6
cm2Mpc = 1/(parsec_in_cm *1e6)
deg2rad = np.pi/180.
rad2deg = 180. / np.pi
rad2ac = rad2deg*3600
sqdegonsky = 129600. /np.pi


def ahav(a):
    '''
    inverse Haversine
    '''
    return np.arccos(1 - 2 * a)
def hav(theta):
    '''
    Haversine
    '''
    return (1 - np.cos(theta))/2.

def ang_sep(ra1, dec1, ra2, dec2):
    '''
    >> dist = ang_sep_xyz(ra1, dec1, ra2, dec2)

    input/output in degree
    input can be:
     - two equal length arrays, we compute distance between elements
     - one array, one scalar, we compute distance between array and scalar

     uses Havesine functions
    '''
    ra1 = np.asarray(ra1*deg2rad, dtype=np.float64)
    ra2 = np.asarray(ra2*deg2rad, dtype=np.float64)
    dec1 = np.asarray(dec1*deg2rad, dtype=np.float64)
    dec2 = np.asarray(dec2*deg2rad, dtype=np.float64)

    radsep = ahav(hav(dec2-dec1) + np.cos(dec1) * np.cos(dec2) * hav(ra2 - ra1))
    
    return radsep / deg2rad

def ang_sep_xyz(ra1, dec1, ra2, dec2):
    '''
    >> dist = ang_sep_xyz(ra1, dec1, ra2, dec2)

    input/output in degree
    input can be:
     - two equal length arrays, we compute distance between elements
     - one array, one scalar, we compute distance between array and scalar

    method uses starutil of astrometry.net
    use starutil_numpy.degrees_between(), to get distance between all
    combinations of ra1, ra2 (ie, a matrix). 
    '''

    # (avoid rounding problems for inits)
    ra1 = np.asarray(ra1, dtype=np.float64)
    ra2 = np.asarray(ra2, dtype=np.float64)
    dec1 = np.asarray(dec1, dtype=np.float64)
    dec2 = np.asarray(dec2, dtype=np.float64)
    
    xyz1 = sutil.radectoxyz(ra1, dec1)
    xyz2 = sutil.radectoxyz(ra2, dec2)

    d2 = np.sum((xyz1-xyz2)**2, axis=1)
    dist_rad = np.arccos(1. - 0.5*d2)

    # convert single element array to float 
    if len(dist_rad) == 1:
        dist_rad = dist_rad[0]
        
    return dist_rad * rad2deg


def one_match(ra0, dec0, ra, dec, rad):
    '''
    matches between one coordindate and an array
    >> m1, d = one_match(ra0, dec0, ra, dec, rad)
    '''
    dist = ang_sep(ra0, dec0, ra, dec)
    index = np.where(dist < rad)[0]
    return index, dist[index]


def one_group(ra0, dec0, ra, dec, link_length, group_in=None, 
              verbose=False, apl=False, cat=False):
    '''
    the most simple friends-of-friends I could think of
    returns indices of (ra,dec) that belong to the group
   >> m1 =  one_group(ra0, dec0, ra, dec, link_length)
   >> m1 =  one_group(None, None, ra, dec, link_length, group_in=[0,3])

    it can return only one group
    input:
     - link_lenght: in deg
     - ra0, dec0: start of matching, to be picked near the center of the group
     - group_in: start with these indices to (ra, dec). If set, ra0, dec0 are ignored 
     - apl=APLpy.FITFigure, used to test the preformance
     - cat=np.recarray, containing ra dec columns, only used if apl keyword is True
     - verbose=False
    '''

    if apl:
        verbose=True
        apl.show_markers(ra0, dec0,
                         edgecolor='red', marker='x', alpha=.8, s=200)            

    # first simply match all object to the center of the group
    if group_in is None:
        m1,m2,d12 =spherematch(ra,dec, ra0,dec0, link_length)
        ss1, d12 = one_match(ra0, dec0, ra, dec, link_length)
        
        if len(m1) < 1:
            if verbose:print ' < 1 elements in group, not matching:', m1
            return m1

        prev_members = m1
        current_members = m1
        if verbose: print 'setting prev_members to first match with ra0,dec0:', prev_members


    else:

        m1, m2, d12 = spherematch(ra, dec, group_in[0], group_in[1], 0.1/3600.)
        prev_members = m1
        current_members = m1
        if verbose: print 'setting prev_members to first match with group_in:', prev_members

        
    if verbose: print ' one_group: #potential ra:', len(ra)
    left = np.arange(len(ra))
    left = np.setdiff1d(left, prev_members)
    nleft =   len(left)
    if verbose:
        print ' one_group: #left to match:', nleft

    if nleft == 0:
        print ' first group already contains everyting'
        return m1
    
    if apl:
        try:
            dum = len(cat)
            apl.show_markers(np.array(cat[prev_members]['ra']),
                         np.array(cat[prev_members]['dec']),
                         edgecolor='cyan', marker='d', alpha=.7, s=350)
            print ' overplotted prev_members in cyan diamonds:', prev_members
            key = raw_input()
        except TypeError:
            print '\n no cat=catalog with ra, dec columns given to overplot coordinates'

    # loop until now new groups are found, or all is in the group
    while len(left):

        if verbose:
            print ' current members', current_members
            print ' # left ',len(left)

        # find new matches
        new_members = np.array([],dtype='<i4')
        #for mm in prev_members:
        m1,m2,d12 =spherematch(ra[left],dec[left],
                               ra[prev_members],dec[prev_members], link_length)
        if len(m1):
            um1 = np.unique(m1)
            new_members  = np.append(new_members,left[um1])

        if verbose:
            print ' new members from match', new_members

        # add new matches, or stop
        if len(new_members) > 0:
            current_members = np.append(current_members, new_members)
            prev_members = new_members
            left = np.setdiff1d(left, prev_members) # remove elements that are in group

            if apl:
                try:
                    dum = len(cat)
                    apl.show_markers(np.array(cat[new_members]['ra']),
                                     np.array(cat[new_members]['dec']),
                                     edgecolor='green', alpha=.7, s=300)
                    print ' overplotted new_members with green circles:', new_members

                except TypeError:
                    print '\n no cat=catalog given to overplot coordinates'
            if verbose: key = raw_input()
            
        else:
            break


    return current_members



def iau_name(ra,dec,prefix='',precision=1, verbose=False):
    '''
    make IAU name from coordinates
    copied from hogg_iau_name (by D.W. Hogg)
    
    >> iau_name(350.95257, -1.1361928, prefix='TDE ',precision=1)
    >> 'TDE 232348.6-010810.3'
    '''
    
    rah= int(floor(ra/15.0))
    ram= floor(((ra/15.0)-double(rah))*60.0)
    ras= (((ra/15.0)-double(rah))*60.0-double(ram))*60.0
    rasformat= '%'+str(precision+3)+ '.'+str(precision)+'f'
    if (precision == 0):
        ras = np.round(ras)
        rasformat= '%2d'

    if verbose: print 'ras, rasformat:', ras, rasformat 

    desgn= '+'
    if dec <0: desgn= '-'

    adec= abs(dec)
    ded= floor(adec)
    dem= floor((adec-double(ded))*60.0)
    des= ((adec-double(ded))*60.0-double(dem))*60.0
    desformat= '%'+str(precision+3)+ '.'+str(precision)+'f'
    if (precision == 0):
        des = np.round(des)
        desformat= '%2d'

    if verbose: print 'des, desformat', des, desformat

    adstr= '%2.2d'%rah + '%2.2d'%ram + rasformat %ras \
           +desgn +'%2.2d'%ded +'%2.2d'%dem +desformat %des

    adstr=adstr.replace(' ','0')
    adstr=adstr.replace(' ','0')

    return prefix+adstr


def lumdis(z, h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    return luminosity distance in cm
    use WMAP+SN+BOA comology
    with omega_m_0=0.3, omega_l_0=0.3, h=0.72
    can be overwritten using keywords
    >> ld_cm = lumdis(z, h=.72, omega_m_0=.3, omega_l_0=.7)
    '''
    cosmo = sjoert_cosmo
    cosmo['h'] = h
    cosmo['omega_M_0'] = omega_m_0
    cosmo['omega_lambda_0'] = omega_l_0
    
    return cospy.distance.luminosity_distance(z, **cosmo)*1e6*parsec_in_cm

def comdis(z, h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    return comosving distance in cm
    use WMAP+SN+BOA comology
    with omega_m_0=0.3, omega_l_0=0.3, h=0.72
    can be overwritten using keywords
    >> cd_cm =  comdis(z, h=.72, omega_m_0=.3, omega_l_0=.7)
    '''
    cosmo = sjoert_cosmo
    cosmo['h'] = h
    cosmo['omega_M_0'] = omega_m_0
    cosmo['omega_lambda_0'] = omega_l_0
    
    return cospy.distance.comoving_distance(z, **cosmo)*1e6*parsec_in_cm

def angdis(z, h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    return comosving distance in cm
    use WMAP+SN+BOA comology
    with omega_m_0=0.3, omega_l_0=0.3, h=0.72
    can be overwritten using keywords
    >> ad_cm =  angdis(z, h=.72, omega_m_0=.3, omega_l_0=.7)
    '''
    cosmo = sjoert_cosmo
    cosmo['h'] = h
    cosmo['omega_M_0'] = omega_m_0
    cosmo['omega_lambda_0'] = omega_l_0
    
    return cospy.distance.angular_diameter_distance(z,
                                                    **cosmo)*1e6*parsec_in_cm

def lum2flux(L, z=None, cm=None, nu=None, band=None,
                     h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    erg/s to Jansky
    >> flux = lum2flux(L, z, nu=1.4e9) # in Jy
    input: 
     - L: luminsity in erg/s
     - z:  reshift
     - nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
     - cm is set, redshift is ignored
    '''
    if not(nu) and not(band):
        print 'please give nu= (in Hz) or band=[FUV, NUV, u,g,r,i,z]'
        return None

    if band is not(None):  nu = _get_nu(nu, band)
    if cm is None: cm = lumdis(z, h=h, omega_m_0=omega_m_0, omega_l_0=omega_l_0)

    return L / (nu * 4*np.pi * cm**2) *1e23 
          
def flux2lum(S, z=None, cm=None, nu=None, band=None,
                     h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    Jansky to erg/s
    >> nuLnu = flux2lum(S, z, nu=None, band=None)

    input: 
     - S: flux in Jansky
     - z:  reshift
     - nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
     - cm is set, redshift is ignored
    can overwrite default cosmology using keywords (see lumdis docstring)
    '''
    if not(nu) and not(band):
        print 'please give nu= (in Hz) or band=[FUV, NUV, u,g,r,i,z]'
        return None

    if band is not(None):  nu = _get_nu(nu, band)
    if cm is None: cm = lumdis(z, h=h, omega_m_0=omega_m_0, omega_l_0=omega_l_0)

    return 4*np.pi * cm**2 * S*1e-23 * nu


def dismod(z, h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    distance modulus: 5*np.log10(lumdis/10.)
    >> dm = dismod(z, h=.72, omega_m_0=.3, omega_l_0=.7)
    '''
    d = lumdis(z, h=h, omega_m_0=omega_m_0, omega_l_0=omega_l_0)/parsec
    return 5*np.log10(d/10.)


abs_const = 4*np.pi* 3631e-23 *  (10*parsec_in_cm)**2

def abs2lum(M, nu=None, band=None):
    '''
    absolute AB manitude to nu*L_nu
    >> lum = abs2lum(-12, nu=1e15)

    input nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
    '''

    if nu is None and band is None:
        print 'please give nu= (in Hz) or band=[FUV, NUV, u,g,r,i,z]'
        return None

    if nu is None:
        nu = _get_nu(band)

    return 10**(-0.4*M) * abs_const * nu


def lum2abs(L, nu=None, band=None):
    '''
    nu Lnu (erg/s) to AB absolute mag
    '''
    if not(nu) and not(band):
        print 'please give nu= (in Hz) or band=[FUV, NUV, u,g,r,i,z]'
        return None
    if not(nu):
        nu = _get_nu(band)    

    return -2.5*np.log10(L/(abs_const*nu))


def flux2mag(flux):
    '''
    >> AB_mag = flux2mag(flux) # flux in Jy
    '''
    return -2.5*np.log10(flux*1e-23) - 48.6

def mag2flux(mag):
    '''
    AB magnitude to Jansky
    '''
    return 10**(-0.4*(mag + 48.6)) *1e23

def mag2lum(M, z, nu=None, band=None,
                      h=.72, omega_m_0=.3, omega_l_0=.7):
    '''
    lum  = mag2lum(22.5, 0.1, band='r')
    convert AB magnitude to luminosity in erg/s
    input nu in Hz, choose from band=[FUV, NUV, u,g,r,i,z]
    can overwrite default cosmology using keywords (see lumdis docstring)
    '''

    return abs2lum(M - dismod(z, omega_m_0=omega_m_0, omega_l_0=omega_l_0), nu, band)

def tundo(Mr):
    '''
    return black hole mass based on SDSS r-band absolute magnitude
    from Tundo et al. (2007)  - 2007ApJ...663...53T
    >> Mbh = tundo(-22)
    '''
    slope = -1.31 / 2.5
    log_mass = (Mr+22)*slope +8.69 
    return 10**log_mass


# --- 
# some filters and functions
sdss_l = np.array([3543, 4770, 6231, 7625, 9134]) # SDSS filter in A

def sdss_nu(k):
    return 3e8 / (sdss_l[k]*1e-10)

def sdss_lambda(k):
    return sdss_l[k]

def _get_nu(band):
    nu = None
    if band == 'FUV': nu = 3e8 / (1528 *1e-10)
    if band == 'NUV': nu = 3e8/ (2271 *1e-10)
    if band == 'u': nu = sdss_nu(0)
    if band == 'g': nu = sdss_nu(1)
    if band == 'r': nu = sdss_nu(2)
    if band == 'i': nu = sdss_nu(3)
    if band == 'z': nu = sdss_nu(4)
    if not(nu):
        print 'please use band=[FUV, ..., z]'
    return nu


def gamma2beta(gamma):
    return (1-gamma**(-2))**(0.5)

def Doppler(gamma=5, i_obs=0.1):
    '''
    Doppler fractor: 1/(gamma * (1-beta*cos(i_obs)))
    >> delta = Doppler(gamma=5, i_obs=0.1)
    '''
    
    beta = gamma2beta(gamma)
    return 1/(gamma * (1-beta * np.cos(i_obs)) )

def Rs(Mbh):
    '''
    return Schwarzschild Radius in cm
    input: Mbh in Msun
    note, only accurate to one decimal :)
    '''
    return 2*Mbh*2e30*6.7e-11/9e16 *1e2 # cm


def schechter(M, h=0.72, paper='Blanton01'):
    '''
    the Schelter function for the 
    possible papers: r-band: Blanton01  (default)
                               K-band: Simth09, Bell03, Cole01, Kochanek01, Jone06

    >> rho  = schechter([-22, -23], h=0.72, paper='Blanton01')
    '''
    if paper== 'Blanton01':

        # wrong r-band?
        #M_s = -20.44 + 5*np.log10(h)
        #alpha = -1.05
        #psi_s = 1.49e-2 *h**3

        # r-band
        M_s = -20.83 + 5*np.log10(h)
        alpha = -1.20
        psi_s = 1.46e-2 *h**3

    if paper== 'Blanton01_gband':
        M_s = -20.04 + 5*np.log10(h)
        alpha = -1.26
        psi_s = 2.06e-2 *h**3

    if paper== 'Gaia':
        M_s = -20.04 + 5*np.log10(h)
        alpha = -1.0
        psi_s = 4e-3 *h**3

    if paper == 'Smith09': # K-band
        M_s = -23.19 + 5*np.log10(h)
        alpha = -0.81
        psi_s = 0.0166 *h**3

    if paper == 'Bell03': # K-band, valid for bright-end only...
        M_s = - 23.33 + 5*np.log10(h)
        alpha = (-0.88)
        psi_s = 0.0149 *h**3

    if paper == 'Cole01': # 2001MNRAS.326..255C
        M_s = - 23.44 + 5*np.log10(h)
        alpha = -0.96
        psi_s = 0.0108 *h**3

    if paper == 'Kochanek01': #2001ApJ...560..566K 
        M_s = - 23.39 + 5*np.log10(h)
        alpha = -1.09
        psi_s = 0.0116 *h**3

    if paper == 'Jones06': #2006MNRAS.369...25J
        M_s = - 23.83 + 5*np.log10(h)
        alpha = -1.16
        psi_s = 10**(-2.126) *h**3

    
    return 0.4*np.log(10) * psi_s * 10.0**(0.4*(alpha+1)*(M_s-M) ) * \
        np.exp( -10.0**( 0.4*(M_s-M) ) )
    
