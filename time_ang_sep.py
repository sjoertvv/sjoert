import time
import numpy as np
import sjoert.stellar
from matplotlib import pyplot as plt

nn = 100000 
ra1 = np.random.rand(nn)-0.5
ra2 =  np.random.rand(nn)-0.5
dec1 =  np.random.rand(nn)-0.5
dec2 =  np.random.rand(nn)-0.5

time0 = time.time()
su_dist = sjoert.stellar.sutil.degrees_between(ra1, dec1, ra2[0], dec2[0])
print 'time for starutil (ra1, dec1, ra2[0], dec2[0])', time.time() - time0

time0 = time.time()
su_dist2 = sjoert.stellar.sutil.degrees_between(ra1[0], dec1[0], ra2, dec2)
print 'time for starutil (ra1[0], dec1[0], ra2, dec2)', time.time() - time0


time0 = time.time()
sj_dist = sjoert.stellar.ang_sep(ra1, dec1, ra2[0], dec2[0])
print 'time for ang_sep', time.time() - time0


time0 = time.time()
sq_dist = np.sqrt( (ra1-ra2[0])**2 + (dec1-dec2[0])**2)
print 'time for simple Euclid', time.time() - time0

delta = su_dist-sj_dist
#plt.clf()
#plt.hist(delta, bins=100)
print 'mean, max difference between mine and starutil', np.mean(delta), np.max(np.abs(delta))

