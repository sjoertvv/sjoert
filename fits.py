'''
small collection of functions for pyFITS tables
(see sjoert.io for I/O of FITS files)

'''
from astropy.io import fits as pyfits

# this will be needed in the future, when STScI stops supporting the stand-alone version
#from astropy.io import fits as pyfits 

def append_field(orig_table, name, data=[], format='D', out_name=''):
    '''
    append a field (a column) to a FITS table 

    new_table = append_field(orig_table, 'name', data=data)

    optional input 
     - data=[]: rows of the new column
     - format='D', default format is Doubble
     - out_name='', if given we write new FITS table to disk

    notes:
     - we assume data is in hdu[1]
     - requires pyfits 3.3 or later
    '''

    #orig_table = orig_hdu[1].data
    orig_cols = orig_table.columns
    
    if len(data)==len(orig_table):
        new_data = data
    else:
        new_data = np.zeros(len(orig_table))

    new_cols = pyfits.ColDefs([pyfits.Column(name=name, format=format, array=new_data)])
    new_hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)
 
    if out_name: 
        new_hdu.writeto(out_name, clobber=True)
    
    return new_hdu.data
