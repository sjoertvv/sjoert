import pyfits

def append_field(orig_hdu, name, data=[], format='D', out_name=''):
    '''
    append a field to a FITS table 

    new_hdu = FITS_append_field(orig_hdu, 'name', data=data)

    optional input 
     - data=[]: rows of the new column
     - format='D', default format is Doubble
     - out_name='', if given we write new FITS table to disk

    note, we assume data is in hdu[1]
    '''
    orig_table = orig_hdu[1].data
    orig_cols = orig_table.columns
    
    if len(data)==len(orig_table):
        new_data = data
    else:
        new_data = np.zeros(len(orig_table))

    new_cols = pyfits.ColDefs([pyfits.Column(name=name, format=format, array=new_data)])
    new_hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)
 
    if out_name: 
        new_hdu.writeto(out_name, clobber=True)
    
    return new_hdu
