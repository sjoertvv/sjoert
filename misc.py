'''
Some functions from the scipy Cookbook
'''

from scipy import *
import numpy as np

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape)/asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))

def bino_rebin(t, a, binwidth):
    '''
    a,t = bino_rebin(t, a, binwidth)
    rebin array a[t] using binwidth of t
    size has to be 2**n
    '''
    nn = len(a)

    if nn != len(t):
        raise ValueError('t and a need to have same dimension')
    if np.log2(nn) / np.floor(np.log2(nn))  !=1:
        print 'len(a):', nn
        raise ValueError('input array needs to be binomial (2**n)')

    dt = t[2] - t[1]
    nbins_float = binwidth / dt
    nbins = 2**np.round(np.log2(nbins_float))

    if nbins > nn/2: nbins = nn
    if nbins < 1: nbins = 1

    new_nn = nn/nbins
    
    return rebin(t, new_nn), rebin(a, new_nn)


