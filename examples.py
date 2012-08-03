'''
not a class but simply a collection of how to's
'''
###!/usr/bin/env python

##########################
# get input from command line #
##########################
import sys, getopt, os

try:  
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "v:b:m:c:hd", ['nologger'])
except getopt.GetoptError:
    print 'unkown parse parameter' 
    #usage()
    sys.exit(2) 

odict = dict(opts)

################
# two window plot #
################
"""
To create plots that share a common axes (visually) you can set the
hspace bewtween the subplots close to zero (do not use zero itself).
Normally you'll want to turn off the tick labels on all but one of the
axes.

In this example the plots share a common xaxis but you can follow the
same logic to supply a common y axis.
"""
#from pylab import *
import matplotlib.pyplot as plt
import numpy as np

t = np. arange(0.0, 2.0, 0.01)
s1 = np.sin(2*np.pi*t)
s2 = np.exp(-t)
s3 = s1*s2
plt.close()

# axes rect in relative 0,1 coords left, bottom, width, height.  Turn
# off xtick labels on all but the lower plot

f = plt.figure()
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(211)

norm = 1e-17
ax1.plot(t, s1, label='L1')
plt.minorticks_on()
plt.legend()
plt.xlabel('X label')


ax2 = plt.subplot(212, sharex=ax1)
plt.plot(t,s2)
plt.xlabel('X label')
plt.yticks(np.arange(0, 1, .2))

#plt.ylabel('Normalized Flux')

plt.minorticks_on()
xticklabels = ax1.get_xticklabels()#+ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

plt.show()




