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
print 'time for starutil', time.time() - time0


time0 = time.time()
sj_dist = sjoert.stellar.ang_sep(ra1, dec1, ra2[0], dec2[0])
print 'time for ang_sep', time.time() - time0

delta = su_dist-sj_dist
plt.clf()
plt.hist(delta, bins=100)
print 'mean, max difference', mean(delta), max(np.abs(delta))

