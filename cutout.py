'''
get subimage from fits image
wcs info in header is updated
'''

import numpy as np
from matplotlib import pyplot as plt
import os,sys,time
import pyfits, pywcs
from sjoert.stellar import iau_name

# input
def cutout(filename='', hdu=None, center=[None, None],im_size=0.1,
           writepdf=False, writefits=False):
    '''
    >> sub_hdu = cutout(filename='file.fits', center=[12,-1.], im_size=0.1,
                                        writepdf=False, writefits=False)
    >> sub_hdu = cutout(hdu=hdu, center=[12,-1.], im_size=0.1,
                                        writepdf=False, writefits=False
    returns:
      pyfits.core.PrimaryHDU 
    input:
     - hdu (pyfits.core.HDUList or pyfits.core.PrimaryHDU ) or filename of fits image
     - center=[ra,dec] (deg)
     - im_size=0.1 (deg)
    optional input:
     - writepdf='name.pdf' or True (write filename_sub.pdf)
     - writefits='name.fits' or True (fileanme+sub.fits

    the creation of new WCS header doesnt seem to work for all cases...
    '''

    if (filename=='') and not(hdu):
        print 'please give filename or hdu'
        return
    if not(hdu) and not(os.path.isfile(filename)):
        print 'file :'+ filename+ '\n not found' 
        return
    
    if not np.isscalar(center[0]):
        print 'please give center=[ra, dec]'
        return

    center = np.array(center)

    # read or index to PrimaryHDU 
    if not(hdu):
        hdu = pyfits.open(filename)[0]
    elif type(hdu) == pyfits.core.HDUList:
        hdu = hdu[0]
    elif type(hdu) != pyfits.core.PrimaryHDU:
        print 'type(hdu:)', type(hdu)
        print 'ERROR: type(hdu) needs to be <pyfits.core.HDUList> or pyfits.core.PrimaryHDU'
        return

    #  get WCS, get image
    wcs = pywcs.WCS(header=hdu.header,  naxis=2) #fobj=hdulist,
    im = hdu.data
    

    print 'image shape', im.shape
    if (len(im.shape) == 3):
        print 'using first layer of image'
        im = im[0,:,:]
        print 'new image shape', im.shape
    if (len(im.shape) == 4):
        print 'using first layer of image'
        im = im[0,0,:,:]
        print 'new image shape', im.shape

    # conver center en overplot coordinates
    pcenter = (wcs.wcs_sky2pix([center], 1))[0]

    # slice the image
    print 'center in pixel coordinates', pcenter

    # use true sky location (creates non-square images)
    #delta_ax1 = (np.abs(wcs.wcs_sky2pix([ [center[0]-im_size, center[1]] ], 1) - pcenter)/2.)[0][0]
    #delta_ax2 = (np.abs(wcs.wcs_sky2pix([ [center[0], center[1]-im_size] ], 1) - pcenter)/2.)[0][1]
    #delta = [delta_ax1, delta_ax2]
    #delta=  np.abs(wcs.wcs_sky2pix( [center-im_size], 1)[0] - pcenter)/2.

    # use pixel scale to define boundaries
    delta = np.abs(0.5*im_size / np.array([hdu.header['CDELT1'], hdu.header['CDELT2']]))

    ax1l = pcenter[0]-delta[0]
    ax1u=pcenter[0]+delta[0]
    ax2l = pcenter[1]-delta[1]
    ax2u=pcenter[1]+delta[1]

    # check limits
    if (ax1l < 0): ax1l = 0
    if (ax1u > im.shape[1]): ax1u = im.shape[1]
    if (ax2l < 0): ax2l = 0
    if (ax2u > im.shape[0]): ax2u = im.shape[0]

    print 'ax2l,ax2u, ax1l,ax1u', ax2l,ax2u,ax1l,ax1u

    subim = im[ax2l:ax2u,ax1l:ax1u ]

    print 'subim shape:',subim.shape
    if min(subim.shape) == 0:
        print 'ERROR. something is wrong with the cutout (center+/-im_size out of bounds?) '
        return None

    # set the base of the output file 
    if (writepdf==True) or (writefits==True):
        srem = ['.FITS', '.FIT', '.fit', '.fits']
        outname = filename
        for srm in srem: outname=outname.split(srm)[0]
        if filename =='':
            print 'warning, not filename given pdf or fits file, using IAU name'
            outname = iau_name(center[0], center[1])
        print 'base for output file:', outname
            

    if writepdf:
        # determine color range, we use arcsinh scaling
        scaled_subim = np.arcsinh(subim)
        vals = scaled_subim.reshape(subim.shape[0]*subim.shape[1])
        #nbins=50
        #h_vals = np.histogram(np.log10(vals), new=False, bins=nbins)
        svals = np.sort(vals)
        slen = len(svals)
        vmin = svals[0.01*slen]
        vmax = svals[0.99*slen] #np.power(10, h_vals[1][nbins-10])
        plt.imshow(scaled_subim, vmin=vmin, vmax=vmax,
                   cmap='hot', interpolation='bilinear')

        # set filename (if not given)
        if not(type(writepdf) is str):
            outfile = outname+'_sub.pdf'
        else: outfile = writepdf
        print 'writing', outfile
        plt.savefig(outfile, format='pdf')

    #make new WCS and save fits
    sub_hdr = hdu.header.copy() #sub_hdu.header

    ## this wrong: slightly rotates the image
    #sub_hdr['CRVAL1'] = center[0]
    #sub_hdr['CRVAL2'] = center[1]
    #sub_hdr['CRPIX1'] = pcenter[0]-ax1l
    #sub_hdr['CRPIX2'] = pcenter[1]-ax2l

    ## keep CRVAL1 the same, but swift its pixel coodindates
    sub_hdr['CRPIX1'] = np.floor(hdu.header['CRPIX1'] - ax1l ) +1
    sub_hdr['CRPIX2'] = np.floor(hdu.header['CRPIX2'] - ax2l ) +1
    print 'orgininal CRVAL (', hdu.header['CRVAL1'], hdu.header['CRVAL2'],\
                              ') in sub image pixel coordindates:', sub_hdr['CRPIX1'] , sub_hdr['CRPIX2'] 

    # note, after this call sub_hdu.header is differnt the input than sub_hdr
    # (eg, NAXIS, NAXIS1, are updated) 
    sub_hdu = pyfits.PrimaryHDU(data=subim, header=sub_hdr)

    # write files
    if writefits:
        if not(type(writefits) is str):
            outfile = outname+'_sub.fits'
        else: outfile = writefits
        print 'writing:', outfile
        sub_hdu.writeto(outfile,clobber=True)


    return sub_hdu

