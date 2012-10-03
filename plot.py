'''
setup matplotlib for IDL-style high quality plots
'''


import re
import os

from matplotlib import rc
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import NullFormatter


def init(fig_width=7,fig_height=7,axes_labelsize=21,
         text_fontsize=20,legend_fontsize=15,
         xtick_labelsize=16,ytick_labelsize=16,
         xtick_minor_size=2,ytick_minor_size=2,
         xtick_major_size=4,ytick_major_size=4,
         subplot_bottom=0.16, subplot_top=.90,
         subplot_left=0.15,subplot_right=0.95,
         golden=False):
    """
    NAME:
       init
    PURPOSE:
       setup a figure for plotting
    INPUT:
       golden=False - use golden ratio (fig_heigh=fig_width/1.618)
       fig_width=7 - width in inches
       fig_height=7 - height in inches
       axes_labelsize=21 - size of the axis-labels
       text_fontsize=20 - font-size of the text (if any)
       legend_fontsize=15 - font-size of the legend (if any)
       xtick_labelsize=16 - size of the x-axis labels
       ytick_labelsize=16 - size of the y-axis labels
       xtick_minor_size - size of the minor x-ticks
       ytick_minor_size - size of the minor y-ticks
       subplot_bottom=0.16 - off set from bottom of x-axis
       subplot_top=0.9 - off set from bottom of x-axis 

    OUTPUT:
       (none)
    HISTORY:
       2009-12-23 - Written - Bovy (NYU)
       2010-01-23 - New defauts more options, backend:ps - Sjoert
    """

    if golden:
        fig_height=fig_width/1.6180
    
    #    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('font',**{'family':'serif','serif':['Times']})
    rc('text', usetex=True)
    plt.minorticks_on()
    
    plt.close()
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': axes_labelsize,
              'text.fontsize': text_fontsize,
              'legend.fontsize': legend_fontsize,
              'xtick.labelsize':xtick_labelsize,
              'ytick.labelsize':ytick_labelsize,
              'figure.figsize': fig_size,
              'xtick.major.size' : xtick_major_size,
              'ytick.major.size' : ytick_major_size,
              'xtick.minor.size' : xtick_minor_size,
              'ytick.minor.size' : ytick_minor_size,
              'figure.subplot.bottom':subplot_bottom,
              'figure.subplot.right':subplot_right,
              'figure.subplot.left':subplot_left,
              'figure.subplot.top':subplot_top}

    plt.rcParams.update(params)
    rc('text.latex', preamble=r'\usepackage{amsmath}')

def print_end(filename,nosticks=False,**kwargs):
    """
    NAME:
       print_end
    PURPOSE:
       saves the current figure(s) to filename
    INPUT:
       filename - filename for plot (with extension)
    OPTIONAL INPUTS:
       format - file-format
    OUTPUT:
       (none)
    HISTORY:
       2009-12-23 - Written - Bovy (NYU)
       2012-01-30 use plt.minticks_on() instead of _add_ticks() - SVV (Radboud)
    """

    #if not(nosticks): _add_ticks()
    plt.minorticks_on()
    if kwargs.has_key('format'):
        plt.savefig(filename,format=kwags['format'])

    else:
        plt.savefig(filename)

             

def _add_ticks():
    """
    NAME:
       _add_ticks
    PURPOSE:
       add minor axis ticks to a plot
    INPUT:
       (none; works on the current axes)
    OUTPUT:
       (none; works on the current axes)
    HISTORY:
       2009-12-23 - Written - Bovy (NYU)
    """
    ax=plt.gca()
    xstep= ax.xaxis.get_majorticklocs()
    xstep= xstep[1]-xstep[0]
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(xstep/5.))
    ystep= ax.yaxis.get_majorticklocs()
    ystep= ystep[1]-ystep[0]
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(ystep/5.))
