'''
reading and writing, mostly fits/ascii tables and record arrays
2010 June - started (Sjoert van Velzen)
2012 Aug - updated docstrings (SVV)
'''

import pyfits
import numpy as np
import pickle

def readpickle(filename):
   '''
   rec = readpickle(filename)
   '''
   f = open(filename,  'r')
   rec = pickle.load(f)
   f.close()
   return rec

def readascii(filename='', names='', comment='#',
              delimiter='', write_pickle=False, silent=False, verbose=False):
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
    while (line[0] == comment) | (len(line) <= 1):
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
            if delimiter:
                scl_del = cl.strip(comment).strip('\n').split(delimiter)
            else:
                scl_del = cl.strip(comment).strip('\n').split()
            
            if  (len(scl_del) == ncols):
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
    >> rec = fits2rec(filename='', ihdu=1, vizmod=False,
                                silent=False, verbose=False)

    use ihdu keyword if data is not in HDU[1]
    vizmod keyword is used to rename '_RAj2000' to 'ra'
    
    (Since PyFITS version 3, HDU table data inherits from numpy.ndarray,
    so a conversion is no longer required)
    '''
    
    if not(filename) :
        print 'give flinename= for output' 
        return
    
    hdulist = pyfits.open(filename)
    tbdata = hdulist[ihdu].data
    if verbose: print 'FITS dtype:', tbdata.dtype
    hdulist.close()

    
    new_dtype = tbdata.dtype # requires Pyfits version >3
    
    # give a better name to Vizier output
    if vizmod:
        new_dtype = []
        for i, name in enumerate(tbdata.dtype.names):
            coln = name
            # get dtype from column data (FITS dtype can be weird)
            # at this stage we actually read the column 
            dt = tbdata[coln][0:1].dtype.descr[0][1] 
            if coln  == '_RAJ2000': coln= 'ra'
            if coln == '_DEJ2000': coln = 'dec'

            if not(silent): print coln        
            # Vizier output always has len(dt)==2
            new_dtype.append((coln, dt)) 


    newrec = np.empty(len(tbdata.field(0)), dtype=new_dtype)

    
    for i, coln in enumerate(newrec.dtype.names): 
        coldata = tbdata.field(tbdata.dtype.names[i]).copy()
        if verbose:
            print 'data:\t', coldata
            print 'dtype:\t', newrec[coln].dtype
            print 'name in rec:\t',coln
        newrec[coln] = coldata

    return newrec


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

    if pyfits.__version__.split('.')[0]>=3:

       tbhdu=pyfits.new_table(rec)

    # if you have an old Pyfits,
    # we have to things the messy way       
    else:
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
    
