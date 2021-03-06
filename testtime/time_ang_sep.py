'''
test different angular seperation functions
'''

import time
import numpy as np
import sjoert.stellar
from matplotlib import pyplot as plt

# we compare to pyspherematch's util.starutil_numpy
import starutil_numpy as sutil

ra, dec = 12, 0
time0 = time.time()
for i in range(1):
    l, b = sutil.radectolb(ra, dec)
print ('starutil', time.time() - time0)
time0 = time.time()
for i in range(1):
    l, b = sjoert.stellar.radectolb(ra, dec)
print ('astropy ', time.time() - time0)


ra, dec = np.random.rand(100000)*300, 90-np.random.rand(100000)*90
time0 = time.time()
l, b = sutil.radectolb(ra, dec)
print ('starutil', time.time()- time0)
time0 = time.time()
l1, b1 = sjoert.stellar.radectolb(ra, dec)
print ('astropy ', time.time() -time0)
print ('')
print ('max difference in l/b', max(abs(l-l1)), max(abs(b-b1)))
key = input('done.')


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


#Haversine functions
deg2rad = np.double(np.pi/180.)
    
def ahav(a):
    return np.arccos(1 - 2 * a)
def hav(theta):
    return (1 - np.cos(theta))/2.

def asep(ra1, dec1, ra2, dec2):
    ra1 = np.asarray(ra1*deg2rad, dtype=np.float64)
    ra2 = np.asarray(ra2*deg2rad, dtype=np.float64)
    dec1 = np.asarray(dec1*deg2rad, dtype=np.float64)
    dec2 = np.asarray(dec2*deg2rad, dtype=np.float64)

    radsep = ahav(hav(dec2-dec1) + np.cos(dec1) * np.cos(dec2) * hav(ra2 - ra1))
    
    return radsep / deg2rad



nn = 200000 
ra1 = np.random.rand(nn)*0.01-0.0001
ra2 =  np.random.rand(nn)*0.01-0.0001
dec1 =  np.random.rand(nn)*0.01-0.0001
dec2 =  np.random.rand(nn)*0.01-0.0001

time0 = time.time()
su_dist = sjoert.stellar.sutil.degrees_between(ra1, dec1, ra2[0], dec2[0])
print ('time for starutil (ra1, dec1, ra2[0], dec2[0])', time.time() - time0)

time0 = time.time()
su_dist2 = sjoert.stellar.sutil.degrees_between(ra1[0], dec1[0], ra2, dec2)
print ('time for starutil (ra1[0], dec1[0], ra2, dec2)', time.time() - time0)

time0 = time.time()
sj_dist = ang_sep_xyz(ra1, dec1, ra2[0], dec2[0])
print ('time for ang_sep_xyz', time.time() - time0)

time0 = time.time()
ha_dist = sjoert.stellar.ang_sep(ra1, dec1, ra2[0], dec2[0])
print ('time for Havesine', time.time() - time0)

time0 = time.time()
sq_dist = np.sqrt( (ra1-ra2[0])**2 + (dec1-dec2[0])**2)
print ('time for simple Euclid', time.time() - time0)

delta = su_dist-sj_dist
#plt.clf()
print ('mean, max difference between sjoert.angsep and starutil', np.mean(delta), np.max(np.abs(delta)))

delta = su_dist-ha_dist
print ('mean, max difference between starutil and Havesin', np.mean(delta), np.max(np.abs(delta)))
#plt.hist(delta, bins=100)

ii = np.argsort(np.random.rand(1000)*nn)
plt.clf()
xx = np.linspace(min(ha_dist), max(ha_dist))
plt.scatter(ha_dist[ii],  sj_dist[ii], color='r', alpha=0.5)
plt.plot(xx,xx)
plt.xscale('log')
plt.yscale('log')
key = input()

plt.clf()
plt.scatter(ha_dist[ii], ha_dist[ii] - sj_dist[ii], color='b', alpha=0.5)
plt.xlabel('Havesin distance (degree)')
plt.ylabel('Delta (degree)')
#plt.scatter(su_dist[ii], su_dist[ii] - sj_dist[ii], color='r', alpha=0.5)

#plt.scatter(sj_dist[ii], sj_dist[ii] - sq_dist[ii], color='b', alpha=0.5)

#plt.scatter(sj_dist[ii], sj_dist[ii] - sq_dist[ii], color='r', alpha=0.5)
plt.xscale('log')
#plt.show()
key = raw_input()


