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


import astropy.cosmology as cospy

# import some useful stuff from pyspherematch's util.starutil_numpy
from starutil_numpy import radectolb, lbtoradec, ra2hmsstring, dec2dmsstring, hmsstring2ra, dmsstring2dec


import numpy as np
from numpy import floor, double


# Some Physical constants (cgs)
h = 6.6261e-27      # Planck's constant (erg s)
c = 2.9979e+10      # Speed of light (cm/s)
k = 1.3807e-16      # Boltzman's (erg/K)
e = 1.6022e-12      # Electronvolt (erg)
G = 6.674e-8        # Gravitational (kg m^3 / kg^2/s^2)
sigma_SB =5.6704e-5 # Stefan-Boltzmann constant (erg/s/cm^2/K^4)

# Some Astro constants
parsec = 3.08568e18 # cm
parsec_in_cm = parsec  
Msun = 1.989e33 # gram
Rsun = 6.96e10   #cm
 
# Conversions
Mpc2cm = parsec_in_cm *1e6
cm2Mpc = 1/(parsec_in_cm *1e6)
deg2rad = np.pi/180.
rad2deg = 180. / np.pi
rad2ac = rad2deg*3600
sqdegonsky = 129600. /np.pi
year_in_sec = 31557600.


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


def lumdis(z, h=.72, Om0=.3):
    '''
    return luminosity distance in cm
    use FlatLambdaCDM
    with Om0=0.3, h=0.72
    can be overwritten using keywords
    >> ld_cm = lumdis(z, h=.72, Om0=.3)
    '''
    return cospy.FlatLambdaCDM(H0=h*100, Om0=Om0).luminosity_distance(z).value *1e6*parsec_in_cm

def comdis(z, h=.72, Om0=.3):
    '''
    return comosving distance in cm
    use FlatLambdaCDM
    with Om0=0.3, h=0.72
    can be overwritten using keywords
    >> cd_cm =  comdis(z, h=.72, Om0=.3)
    '''
    return cospy.FlatLambdaCDM(H0=h*100, Om0=Om0).comoving_distance(z).value *1e6*parsec_in_cm

def angdis(z, h=.72, Om0=.3):
    '''
    return comosving distance in cm
    use FlatLambdaCDM
    with Om0=0.3, h=0.72
    can be overwritten using keywords
    >> ad_cm =  angdis(z, h=.72, Om0=.3)
    '''
    
    return cospy.FlatLambdaCDM(H0=h*100, Om0=Om0).angular_diameter_distance(z).value *1e6*parsec_in_cm

def pc2as(z, h=.72, Om0=.3):
    '''
    convert parsec to arcsec
    >> as = pc2mas(0.1) * 10
    '''
    return 1/(angdis(z, h, Om0)/parsec_in_cm)/np.pi*180*3600



def lum2flux(L, z=None, cm=None, nu=None, band=None,
                     h=.72, Om0=.3):
    '''
    erg/s to Jansky
    >> flux = lum2flux(L, z, nu=1.4e9) # in Jy
    input: 
     - L: luminsity in erg/s
     - z:  reshift
     - nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
     - cm if set, redshift is ignored

    note, no K-correction
    '''
    if not(nu) and not(band):
        print 'please give nu= (in Hz) or band=[FUV, NUV, u,g,r,i,z]'
        return None

    if band is not(None):  nu = _get_nu(nu, band)
    if cm is None: cm = lumdis(z, h=h, Om0=Om0)

    return L / (nu * 4*np.pi * cm**2) *1e23 

def lum2mag(L, z=None, cm=None, nu=None, band=None,
                     h=.72, Om0=.3):
    '''
    erg/s to AB mag
    >> flux = lum2flux(L, z, nu=1.4e9) # in Jy
    input: 
     - L: luminsity in erg/s
     - z:  reshift
     - nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
     - cm if set, redshift is ignored

    note, no K-correction 
    '''

    return flux2mag(lum2flux(L, z=z, cm=cm, nu=nu, band=band,
                     h=h, Om0=Om0))


def flux2nuFnu(S, nu):
    '''
    nuFnu = flux2nuFnu(S, nu)
    Jansky to erg/s/cm^2
    '''
    return S*1e-23*nu
          
def mag2nuFnu(mag, nu):
    '''
    nuFnu = mag2nuFnu(ABmag, nu)
    AB mag to erg/s/cm^2
    '''
    return mag2flux(mag)*1e-23*nu

def flux2lum(S, z=None, cm=None, nu=None, band=None,
                     h=.72, Om0=.3):
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
    if cm is None: cm = lumdis(z, h=h, Om0=Om0)

    return 4*np.pi * cm**2 * S*1e-23 * nu


def dismod(z, h=.72, Om0=.3):
    '''
    distance modulus: 5*np.log10(lumdis/10.)
    >> dm = dismod(z, h=.72, Om0=.3)
    '''
    d = lumdis(z, h=h, Om0=Om0)/parsec
    return 5*np.log10(d/10.)


abs_const = 4*np.pi* 3631e-23 *  (10*parsec_in_cm)**2

def abs2mag(M, z=None):
    '''
    absolute AB manitude to nu*L_nu
    >> mag = abs2flux(-12, z, nu=1e15)

    input nu in Hz, or choose from band=[FUV, NUV, u,g,r,i,z]
    '''
    lum=abs2lum(M, nu=1.)

    return lum2mag(lum, z=z, nu=1.)


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
                      h=.72, Om0=.3):
    '''
    lum  = mag2lum(22.5, 0.1, band='r')
    convert AB magnitude to luminosity in erg/s
    input nu in Hz, choose from band=[FUV, NUV, u,g,r,i,z]
    can overwrite default cosmology using keywords (see lumdis docstring)
    '''

    return abs2lum(M - dismod(z, Om0=Om0), nu, band)


def mag2diff(mag, mag_err, baseline, baseline_err):
    '''
    >> diff_mag, diff_mag_err = mag2diff(mag, mag_err, baseline, baseline_err)

    compute the diffence between two magnitudes, plus the uncertainty 
    '''

    flux_baseline = 10**(-0.4*(baseline))
    flux = 10**(-0.4*mag)
    flux_diff = flux - flux_baseline
    flux_diff_err = np.sqrt((flux*mag_err)**2 + (flux_baseline*baseline_err)**2)
    mag_diff = -2.5*np.log10(flux_diff)
    mag_diff_err = flux_diff_err / flux_diff
    return mag_diff, mag_diff_err


def tundo(Mr):
    '''
    return black hole mass based on SDSS r-band absolute magnitude
    from Tundo et al. (2007ApJ...663...53T)

    note, this relation gives ~2 higher masses compared to McConnel & Ma

    >> Mbh = tundo(-22)
    '''
    slope = -1.31 / 2.5
    log_mass = (Mr+22)*slope +8.69 
    return 10**log_mass

def MBH_scaling(bulge_mass, scat=0, 
        do_McConnellMa=True,
        do_HarningRix=False, 
        do_KormendyHo=False):
    '''
    Container with scaling relations between bulge mass (dynamical or stellar)
    and black hole mass. Input and output in solar mass (linear scale).

    do_McConnellMa (default):
        McConnell & Ma (2013), based on stellar mass of the bulge

    do_HarningRix:
      Harning and Rix (2008) relation with bulge dynamical mass.
    
    do_KormendyHo:
        Kormendy & Ho (2014), Eq. 10, based on fixed K-band M/L


    optional input:
     scat=0.2 to add 0.2 dex scatter (default is zero)
     
    '''

    if scat>0:
        MBH_scat = np.random.normal(0, scat, size=len(bulge_mass))
    else:   
        MBH_scat = 0.
    
    if do_KormendyHo:
        norm = 8.69 # +/- 0.05
        slope = 1.16 # +/- 0.08
    elif do_HarningRix:
        norm = 8.20 # +/- 0.1
        slope = 1.12 # +/- 0.06
    elif do_McConnellMa:
        norm = 8.46 
        slope = 1.05

    return 10**(norm + MBH_scat + slope* np.log10(bulge_mass/1e11))


# --- 
# some filters and functions (moved to sdss.py)
from sdss import sdss_nu, sdss_lambda
from sdss import get_nu as _get_nu


def gamma2beta(gamma):
    return (1-gamma**(-2))**(0.5)

def Doppler(gamma=5, i_obs=0.1):
    '''
    Doppler fractor: 1/(gamma * (1-beta*cos(i_obs)))
    >> delta = Doppler(gamma=5, i_obs=0.1)
    '''
    
    beta = gamma2beta(gamma)
    return 1/(gamma * (1-beta * np.cos(i_obs)) )


def ev2K(ev):
    '''
    convert black body temperature from eV to K
    '''
    return ev*e/k
def ev2Hz(ev):
    '''
    convert frequency in eV to Hz
    '''
    return ev*e/h

def Planck(nu=None, T=None):
    '''
    >> I = Planck(nu=1e14, T=1e4)

    return black body intensity (power per unit area per solid angle per frequency)
    '''
    if (nu is None) or (T is None):
        print 'ERROR, please give input: '
        print 'Planck(nu=nu, T=T)' 
        return np.nan
    
    return 2*h/c**2 * nu**3 / (np.exp(h*nu/(k*T))-1) 
        
#return (2*h*c**2 / wave**5) / (np.exp(h*c/(wave*k*T))-1) 

def dPlanckdT(nu=None, T=None):
    if (nu is None) or (T is None):
        print 'ERROR, please give input: '
        print 'Planck(nu=nu, T=T)'
        return np.nan

    exp_part = np.exp(h*nu/(k*T))
    
    return 2*h/c**2 * h*nu/k*exp_part / (T**2 * (exp_part-1)**2)

def dPlanckdnu(nu=None, T=None):
    if (nu is None) or (T is None):
        print 'ERROR, please give input: '
        print 'Planck(nu=nu, T=T)'
        return np.nan
    hkT=h/(k*T)
    exp_part = np.exp(nu*hkT)
    
    return 2*h/c**2 * (3*nu**2/exp_part  + nu**3/(hkT*(exp_part-1)**2))



def Rs(Mbh):
    '''
    return Schwarzschild Radius in cm
    input: Mbh in Msun
    note, only accurate to one decimal :)
    '''
    return 2*Mbh*2e30*6.7e-11/9e16 *1e2 # cm

def Msigma(sigma):
    '''
    Kormendy and Ho (2014) M-sigma relation 
    give sigma in km/s, return Mbh in Msun
    '''
    return 10**8.5*(sigma/200.)**4.4 # in Msun

def BH_influence(sigma):
    '''
    give sigma in km/s, return sphere of influence in cm
    '''
    return G * Msun * Msigma(sigma) / (sigma*1e5)**2 # this is propto sigma^2.4

def schechter(M, h=0.72, paper='Blanton03'):
    '''
    the Schelter function for the 
    possible papers: r-band: Blanton01  (default)
                               K-band: Simth09, Bell03, Cole01, Kochanek01, Jone06

    >> rho  = schechter([-22, -23], h=0.72, paper='Blanton01')
    '''
    # r-band Blanton01 (z=0)
    if paper== 'Blanton01':

        M_s = -20.83 + 5*np.log10(h)
        alpha = -1.20
        psi_s = 1.46e-2 *h**3

    # r-band shifted to z=0.1
    if paper== 'Blanton03':

        M_s = -20.44 + 5*np.log10(h)
        alpha = -1.05
        psi_s = 1.49e-2 *h**3

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
  

def deVau(r, r_mid=1):
    '''
    De Vaucouleurs
    >> I = DeVau(r, r50)
    '''

    return np.exp(-7.669*(r/r_mid)**(1/4.)-1)


'''
next three functions from starutil_numpy of Dustin Lang
'''
import datetime
from datetime import datetime as dt
import time
def jdtomjd(jd):
        return jd - 2400000.5
def mjdtojd(mjd):
        return mjd + 2400000.5
def mjdtodate(mjd):
        jd = mjdtojd(mjd)
        return jdtodate(jd)
def jdtodate(jd):
        unixtime = (jd - 2440587.5) * 86400. # in seconds
        return datetime.datetime.utcfromtimestamp(unixtime)
def mjdtoyear(mjd):
    if np.isscalar(mjd):
        return datetoyear(mjdtodate(mjd))
    else:
        out = np.zeros(len(mjd))
        for i in range(len(mjd)):
            out[i] = datetoyear(mjdtodate(mjd[i]))
        return out
def timedeltatodays(dt):
        return dt.days + (dt.seconds + dt.microseconds/1e6)/86400.
def datetomjd(d):
        d0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
        dt = d - d0
        # dt is a timedelta object.
        return timedeltatodays(dt)


def datetoyear(date):
    '''
    import is datetime object
    from http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    '''

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

