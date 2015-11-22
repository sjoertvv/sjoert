'''
setup matplotlib for high-quality plots
'''
import re
import os
import numpy as np
     
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
         golden=False, usetex=True):
    """
    NAME:
       init
    PURPOSE:
       setup a figure for plotting
       use golden ratio by default
    INPUT:
       golden=False - serif on the axis (IDL-like)
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
       2010-01-23 - New defauts more options - Sjoert
    """
    fig_height=fig_width/1.6180
    if golden:    
        rc('font',**{'family':'serif','serif':['Times']})    
    else: 
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica','Bitstream Vera Sans', 'sans-serif']})

    rc('text', usetex=usetex)
    
    plt.close()
    fig_size =  [fig_width,fig_height]
    params = { 'backend':'ps',
       'axes.labelsize': axes_labelsize,
       'font.size': text_fontsize,
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
    
    if not golden:
      rc('text.latex', preamble='\usepackage{sfmath} \usepackage{amsmath}') 
    else:
      rc('text.latex', preamble='\usepackage{amsmath}')
    plt.rcParams.update(params)


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


import scipy as sc
from scipy import special

def density_contours(x,y, levels=[0.95],colors=['k'], bins=[20,20],
    xrange=None, yrange=None, cmap=None):
  '''
  overplot density contours
  optional input:
   - xrange, yrange
   - bins=[20,02]
   - levels=[0.95], levels are probability masses contained within the contour
   - colors=['k']
   - cmap=None : if len(levels) != len(colors), we use the cubehelix_r colormap

  sniplets taken from Jo Bovy's bovy_plot.py
  '''

  if xrange is None:
    xrange = [min(x),max(x)]
  if yrange is None:
    yrange = [min(y),max(y)]
  
  if type(colors) is str:
    colors = [colors]
  if type(levels) is float:
    levels = [levels]
  
  if (len(colors) != len(levels)) and (cmap is None):
    colors=None
    cmap = 'cubehelix_r'
  
  data= sc.array([x,y]).T
  hist, edges= sc.histogramdd(data,bins=bins,range=[xrange,yrange])
  X = hist.T

  extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]]
  aspect= (xrange[1]-xrange[0])/(yrange[1]-yrange[0])

  #  make umulative 
  sortindx= sc.argsort(X.flatten())[::-1]
  cumul= sc.cumsum(sc.sort(X.flatten())[::-1])/sc.sum(X.flatten())
  cntrThis= sc.zeros(sc.prod(X.shape))
  cntrThis[sortindx]= cumul
  cntrThis= sc.reshape(cntrThis,X.shape)

  # plot
  plt.contour(cntrThis, levels, colors=colors, cmap=cmap, extent=extent, aspect=aspect)



def bovy_dens2d(x,y, *args,**kwargs):
  '''
  wrapper around bovy_dens2d
  hist_2d = bovy_dens2d(x,y, *args,**kwargs)
  '''

  try:
    import bovy_plot
  except ImportError:
    print 'Import bovyplot failed'

  if kwargs.has_key('xrange'):
      xrange=kwargs['xrange']
      kwargs.pop('xrange')
  else:
      xrange=[x.min(),x.max()]
  if kwargs.has_key('yrange'):
      yrange=kwargs['yrange']
      kwargs.pop('yrange')
  else:
      yrange=[y.min(),y.max()]
  ndata= len(x)
  if kwargs.has_key('bins'):
      bins= kwargs['bins']
      kwargs.pop('bins')
  else:
      bins= round(0.3*sc.sqrt(ndata))
  if kwargs.has_key('aspect'):
      aspect= kwargs['aspect']
      kwargs.pop('aspect')
  else:
      aspect= (xrange[1]-xrange[0])/(yrange[1]-yrange[0])
  if kwargs.has_key('weights'):
      weights= kwargs['weights']
      kwargs.pop('weights')
  else:
      weights= None
  if kwargs.has_key('levels'):
      levels= kwargs['levels']
      kwargs.pop('levels')
  else:
      levels= special.erf(0.5*sc.arange(1,4))

  hh_2d, edges= sc.histogramdd(sc.array([x, y]).T,
                                          bins=bins, range=[xrange ,yrange])

  bovy_plot.bovy_dens2d(hh_2d.T,
          contours=True,levels=levels,cntrmass=True,
          cmap='gist_yarg',origin='lower',
          xrange=xrange, yrange=yrange, aspect=aspect,
          interpolation='nearest',  retCumImage=True, **kwargs)

  return hh_2d
    
def rtext(line,x,y,s, golden=False, **kwargs):
  '''
  add text to line
  rtext(line,x,y,s, **kwargs)

  input line argument is ax.lines[0] (and ax = plt.gca())
  
  taken from http://stackoverflow.com/questions/17252790/matplotlib-adding-text-with-more-than-one-line-adding-text-that-can-follow-the
  '''
  from scipy.optimize import curve_fit
  xdata,ydata = line.get_data()
  dist = np.sqrt((x-xdata)**2 + (y-ydata)**2)
  dmin = dist.min()
  TOL_to_avoid_rotation = 0.3
  if dmin > TOL_to_avoid_rotation:
      r = 0.
  else:
      index = dist.argmin()
      xs = xdata[ [index-2,index-1,index,index+1,index+2] ]
      ys = ydata[ [index-2,index-1,index,index+1,index+2] ]
      def f(x,a0,a1,a2,a3):
          return a0 + a1*x + a2*x**2 + a3*x**3
      popt, pcov = curve_fit(f, xs, ys, p0=(1,1,1,1))
      a0,a1,a2,a3 = popt
      ax = plt.gca()
      derivative = (a1 + 2*a2*x + 3*a3*x**2)
      derivative /= ax.get_data_ratio()
      r = np.arctan( derivative )

  rot = np.rad2deg(r)
  if golden: 
      rot /= 1.6180
      
  return plt.text(x, y, s, rotation=rot, **kwargs)
