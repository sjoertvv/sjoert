'''
reading and writing, mostly fits/ascii tables and record arrays
2010 June - started (Sjoert van Velzen)
2012 Aug - updated docstrings (SVV)
'''

import pyfits
import numpy as np
import pickle

def readascii(filename='', names='', comment='#',
              delimiter='', write_pickle=False, silent=False):
    '''
    read asciitable, return record array
    >> rec = readascii(filename='', names='',formats='', comment='#',
                                  delimiter='', write_pickle=False, silent=False):
                                  
    if header contains a line that equals the number of columns
    of the first row of data we use this to label the columns

    optinal keywords:
     names: names of cols (optional)
     write_pickple: name of the pickle, if True we write filename+.pickple
    '''
    
    try:         
        f = open(filename, 'r')
    except IOError:
        print 'file not found'
        raise(IOError)
    
    # read through header, find number of columns
    com_lines = []
    line = f.readline()
    com_lines.append(line)
    while line[0] == comment: 
        line = f.readline() #read untill we reach data
        com_lines.append(line)
    sline = line.strip(comment).split()
    if (delimiter !=''):
        sline = line.strip(comment).split(delimiter)
    ncols = len(sline)
    f.close() 
    if not silent: print 'first line of data:', str(sline)
    if not silent: print 'number of colums:', ncols

    #Try to get col names from header 
    if (names == '') and (len(com_lines)>1):
        com_lines.reverse()
        for cl in com_lines[1:]:
            scl = cl.strip(comment).split()
            if delimiter:
                scl_del = cl.strip(comment).split(delimiter)
            else:
                scl_del = scl
            if  (len(scl) == ncols): # & (cl.split()[0][0] == comment) ):
                names = scl
                if not silent: print 'using col names from header:', names
                break
            if (len(scl_del) == ncols):
                names = scl_del
                if not silent: print 'using col names from header:', names
                break
            
    if len(names) != ncols:
        #print 'making columns with default names'
        names= [] 
        for i in range(ncols):
            names.append('field-'+np.str(i))

    # Read file, using loadtxt or recfromtxt
    dd = np.recfromtxt(filename, delimiter=delimiter, names=names)
    # notes:
    # recform finds formats automatically, 
    # while genfromtxt and loadtxt do not

    # Make a pickle
    if write_pickle:
        if not(type(write_pickle) is str): fout = filename+'.pickle'
        else: fout = write_pickle
        if not silent: print 'writing:', fout
        f = open(fout, 'w')
        pickle.dump(dd,f)
        f.close()
    return dd

     
def fits2rec(filename='', silent=False, ihdu=1, vizmod=False, verbose=False):
    '''
    read from disk and convert fits table to python record array
    >> rec = fits2rec(filename='', silent=False,
                                ihdu=1, vizmod=False, verbose=False)

    use ihdu keyword if data is not in HDU[1]
    vizmod keyword is used to rename '_RAj2000' to 'ra'
    fails if fits table contains vector data 
    '''
    
    if not(filename) :
        print 'give flinename= for output' 
        return
    
    hdulist = pyfits.open(filename)
    tbdata = hdulist[ihdu].data
    cols = tbdata.dtype
    hdulist.close()

    
    # give a better name to Vizier output
    rec_names = np.array(cols.names)
    if vizmod:
        for i,cc in enumerate(cols.names):
            if cc == '_RAJ2000': rec_names[i] = 'ra'
            if cc == '_DEJ2000': rec_names[i] = 'dec'

    pre_dtype = []
    for i, coln in enumerate(cols.names):
        col_data  = np.asarray(tbdata.field(coln))
        pre_dtype += [(rec_names[i], col_data.dtype)]

    dd = np.empty(len(tbdata.field(0)), dtype=np.dtype(pre_dtype))
    
    for i, coln in enumerate(cols.names): 
        if not(silent): print coln        
        coldata = tbdata.field(coln).copy()
        if verbose:
            print coldata
            print dd[coln].dtype
        dd[rec_names[i]] = coldata

    return dd


def readfits(filename='', ihdu=1):
    '''
    >> tbdata = readfits(filename='', ihdu=1)
    '''

    if not(filename) :
        print 'give flinename= for output' 
        return
    
    hdulist = pyfits.open(filename)
    tbdata = hdulist[ihdu].data
    cols = hdulist[ihdu].columns  
    hdulist.close()
    
    for colname in cols.names: print colname
    return tbdata

def rec2fits(rec=None, filename=''):
    '''
    write records array to fits table
    >> rec2fits(rec, filename)
    '''
    if not(filename) :
        print 'give flinename= for output' 
        return
    
    clist = []
    for dt in rec.dtype.descr:
        name = dt[0]
        form = dt[1]
        print name, form,
        fits_form = 'D'
        if len(form.split('i')) >1:
            fits_form = 'K'
        if len(form.split('S'))>1:
            fits_form='A'+form.split('S')[1]
        print fits_form
        cc = pyfits.Column(name=name, format=fits_form, array=np.array(rec[name]))
        clist.append(cc)

    cols = pyfits.ColDefs(clist) # make class that contains all colums
    tbhdu=pyfits.new_table(cols)
    print 'writing:', filename
    tbhdu.writeto(filename, clobber=1)
    return

    
def writerec(rec=None, filename='', delimiter='\t', wtype='w'):
    '''
    write all collumns of a records array to a ascii file
    >> writerec(rec=rec, filename=filename, delimiter='\t', wtype='w')

    the column names are printed to the header
    (see rec2fits for writing to fits table)
    '''

    f = open(filename, wtype)
    print 'writing', filename
    
    keys =rec.dtype.names # column names

    # make the delimiter vector (last has to be empty)
    delimiter = np.repeat(delimiter, len(keys))
    delimiter[-1] = ''
    
    # write header, unless we append
    if wtype != 'a':
        f.writelines('# ')
        for i, k in enumerate(keys):
            f.writelines(str(k) + delimiter[i])
        f.writelines(' \n')
            
    #write body
    for r in rec:
        for i, k in enumerate(keys):
            f.writelines(str(r[k]) + delimiter[i])
        f.writelines('\n')
    f.close()

def rec2ascii(rec=None, filename='', delimiter='\t', wtype='w'):
    '''
    write all collumns of a records array to a ascii table
    >> rec2ascii(rec=rec, filename=filename, delimiter='\t', wtype='w')
    
    use wtype='a' to append
    (see rec2fits for writing to fits table)
    '''
    writerec(rec, filename, delimiter='\t', wtype='w')
    return

def writecols(cols=[],filename='',  names=[], delimiter='\t', wtype='w'):
    '''
    write al list of columns to an ascii file
    >> writecols(cols=[col1, col2],filename=filename,
                         names=['ra', 'flux'], delimiter='\t', wtype='w')
    '''

    # the lengths
    ncols = len(cols)
    nrows = len(cols[0])
    for i, c in enumerate(cols[-1:]):
        if (len(c) != nrows):
            print 'column ',i, ' not of same length as first column. ', len(c), nrows
            return

    f = open(filename, wtype)
    print 'writing:\n', filename

    # make the header
    if len(names):
        if len(names) != ncols:
            print 'names not equal to number of collumns:', ncols, names
            return
        f.writelines('#')
        for n in names:
            f.writelines(' '+n)
        f.writelines('\n')
    
    
    # make the delimiter vector (last has to be empty)
    delimiter = np.repeat(delimiter, ncols)
    delimiter[-1] = ''

    #write body
    for i in range(nrows):
        for j, c in enumerate(cols[0:]):
            f.writelines(str(c[i]) + delimiter[j])
        f.writelines('\n')
    f.close()
    
