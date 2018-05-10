'''
reading and writing, mostly fits/ascii tables and record arrays

2010 June - started (Sjoert van Velzen)
2012 Aug - updated docstrings (SVV)
2017 Nov - fun with json 
2018 Apr - updates from astropy 3.0
'''

from astropy.io import fits as pyfits
import numpy as np
import pickle
import json
from six import string_types
import sys

def json2rec(jin, silent=False, verbose=False):
    '''
    take a list of dictionaries from json.read() or a filename, 
    read all the columns and make one big python rec array
    '''
    if isinstance(jin, string_types):    
        jdict_list = json.loads(open(jin).read())
    else:
        jdict_list = jin

    if type(jdict_list) != list:
        print('cant work with this, file is not a list; consider going a level deeper')
        return None

    # get all the keys and check if they are string or float, or another dict
    dt_dict = {}
    for jd in jdict_list:
        for ku in jd.keys():
            k = str(ku)
            
            # check the type            
            this_type = type(jd[k])
            
            # check if key contain str 
            if isinstance(jd[k], string_types):
                # check of the str uncodes a number
                try:
                    yup = float(jd[k])        
                    dt = 'f8'                                            
                except ValueError:
                    dt = 'U40'#+str(len(jd[k]))
            elif this_type==bool:
                dt = 'i8'
            # make an exception for entries with two values
            elif this_type ==list:
                if len(jd[k])==2:
                    dt = 'f8'
            else:
                dt = '_nofloat'
            if verbose:
                print('key           :', k)
                print('data          :', jd[k])
                print('current type  :',type(jd[k]))
                print('use type      :', dt)
            
            # add to our dict of dtypes
            if not k in dt_dict:
                dt_dict[k] = dt
            else:
                
                # if we already have this key, make sure the string type is long enough    
                if dt_dict[k][0]=='S':
                    if len(jd[k])>float(dt_dict[k][1:]):
                        dt_dict[k] = 'S'+str(len(jd[k]))            
                
                # or if the new dtype is different from the old one, we can't use this key for our table
                # (can this even happen in proper json tables?)
                elif (dt != dt_dict[k]) and (dt_dict[k]!='_skip'):
                    if verbose or not silent:
                        print('warning, key with different dtype, skipping this one:', k, jd[k], type(jd[k]))
                        print('current type', dt_dict[k])
                        print('new proposed type', dt)
                    dt_dict[k] = '_skip'

    if verbose or not(silent):
        print('# of keys', len(dt_dict.keys()))
    if verbose:
        print(dt_dict)
   
    # make the dtype for the rec array
    dt_list = []
    final_keys = []
    for k in dt_dict.keys():
        if dt_dict[k][0]!='_':
            dt_list.append((k,dt_dict[k]))
            final_keys.append(k)

    if verbose:
        print('final dtype:\n',dt_list)

    recarr = np.zeros(len(jdict_list), dtype=dt_list)    

    # loop over all json entries, read
    for i, jd in enumerate(jdict_list):
        for k in final_keys:
            if jd.get(k) is not None:
                if type(jd.get(k)) == list:                    
                    recarr[i][k] = (float(jd[k][0])+float(jd[k][1]))/2. # dangerous
                else:
                    recarr[i][k] = jd.get(k)

            

    return recarr


def readpickle(filename):
   '''
   rec = readpickle(filename)
   '''
   f = open(filename,  'r')
   rec = pickle.load(f)
   f.close()
   return rec

def readascii(filename='', lines=None, names='', comment='#',
              delimiter='', write_pickle=False, silent=False, verbose=False):
    '''
    read asciitable, return record array
    >> rec = readascii(filename='', names=['a', 'b'], comment='#',
                                  delimiter='', write_pickle=False, silent=False):                                  

    if header contains a line that equals the number of columns
    of the first row of data we use this to label the columns

    alternative input:
     lines = array of lines (eg, from f.open(filename).readlines())

    optinal keywords:
     names: names of cols (optional)
     write_pickple: name of the pickle, if True we write filename+.pickple
    '''
    
    if filename:
        try:         
            f = open(filename, 'r')
        except IOError:
            raise(IOError(filename))
    
    # read through header, find number of columns
    com_lines = []
    
    if filename:
        lines = f.readlines()

    if len(lines)==0:
        print('empty file')
        return None
    l=0
    line = lines[0]
    com_lines.append(line)

    
    while ((line[0] == comment) | (len(line) <= 1)) and (l<len(lines)):
        line = lines[l] #read untill we reach data
        com_lines.append(line)
        l+=1

    sline = line.strip(comment).split()
    if (delimiter !=''):
        sline = line.strip(comment).split(delimiter)
    ncols = len(sline)
    if filename:
        f.close() 
    
    if not silent: print('first line of data:', str(sline))
    if not silent: print('number of colums:', ncols)

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
                if not silent: print('using col names from header:', names)
                break
            
    if len(names) != ncols:        
        names= [] 
        for i in range(ncols):
            names.append('field-'+np.str(i))

    # Read file
    dd = np.recfromtxt(lines, delimiter=delimiter, names=names)

    # notes:
    # recform finds formats automatically, 
    # while genfromtxt and loadtxt do not

    # Make a pickle
    if write_pickle:
        if not(type(write_pickle) is str): fout = filename+'.pickle'
        else: fout = write_pickle
        if not silent: print('writing:', fout)
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
        print('give flinename= for output')
        return
    
    hdulist = pyfits.open(filename)
    tbdata = hdulist[ihdu].data
    if verbose: print('FITS dtype:', tbdata.dtype)
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

            if not(silent): print(coln)
            # Vizier output always has len(dt)==2
            new_dtype.append((coln, dt)) 


    newrec = np.empty(len(tbdata.field(0)), dtype=new_dtype)

    
    for i, coln in enumerate(newrec.dtype.names): 
        coldata = tbdata.field(tbdata.dtype.names[i]).copy()
        if verbose:
            print ('data:\t', coldata)
            print ('dtype:\t', newrec[coln].dtype)
            print ('name in rec:\t',coln)
        newrec[coln] = coldata

    return newrec


def readfits(filename='', ihdu=1):
    '''
    >> tbdata = readfits(filename='', ihdu=1)
    '''

    if not(filename) :
        print('give flinename= for output')
        return
    
    hdulist = pyfits.open(filename)
    tbdata = hdulist[ihdu].data
    cols = hdulist[ihdu].columns  

    hdulist.close()
    for colname in cols.names: print(colname)
       
    return tbdata

def rec2fits(rec=None, filename=''):
    '''
    write records array to fits table
    >> rec2fits(rec, filename)
    '''
    if not(filename) :
        print('give flinename= for output' )
        return

    # try to keep things working for older version of pyFITS <-- no longer possible with astropy.io.fits
    #if pyfits.__version__.split('.')[0]>=3:
    #    if pyfits.__version__.split('.')[1]>=3:
    #        tbhdu=pyfits.BinTableHDU.from_columns(rec) # the way to go
    #    else:
    #        tbhdu=pyfits.new_table(rec) # yields DeprecationWarning as of v3.3
    # if you have an very old Pyfits,
    # we have to things the messy way       
    # else:
    #     clist = []
    #     for dt in rec.dtype.descr:
    #         name = dt[0]
    #         form = dt[1]
    #         print name, form,
    #         fits_form = 'D'
    #         if len(form.split('i')) >1:
    #             fits_form = 'K'
    #         if len(form.split('S'))>1:
    #             fits_form='A'+form.split('S')[1]
    #         print fits_form
    #         cc = pyfits.Column(name=name, format=fits_form, array=np.array(rec[name]))
    #         clist.append(cc)

    #     cols = pyfits.ColDefs(clist) # make class that contains all colums
    #     tbhdu=pyfits.new_table(cols)
    

    tbhdu=pyfits.BinTableHDU.from_columns(rec) # the way to go
    
    tbhdu.writeto(filename, overwrite=True)
    print ('written:', filename)
    return

    
def writerec(rec=None, filename='', delimiter='\t', wtype='w'):
    '''
    write all collumns of a records array to a ascii file
    >> writerec(rec=rec, filename=filename, delimiter='\t', wtype='w')

    the column names are printed to the header
    (see rec2fits for writing to fits table)
    '''

    f = open(filename, wtype)
    print('writing', filename)
    
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

def writecols(cols=[],filename='',  names=[], delimiter='\t', wtype='w', format_str='{0:7.6e}'):
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
            print('column ',i, ' not of same length as first column. ', len(c), nrows)
            return

    f = open(filename, wtype)
    print('writing:\n', filename)

    # make the header
    if len(names):
        if len(names) != ncols:
            print('names not equal to number of collumns:', ncols, names)
            return
        f.writelines('# ')
        for n in names:
            f.writelines('{0:8}'.format(n)+delimiter)
        f.writelines('\n')
    
    
    # make the delimiter vector (last has to be empty)
    delimiter = np.repeat(delimiter, ncols)
    delimiter[-1] = ''

    #write body
    for i in range(nrows):
        for j, c in enumerate(cols[0:]):
            f.writelines((format_str.format(c[i])) + delimiter[j])
        f.writelines('\n')
    f.close()
    

def col_formatter(rcol, prec=7):
    '''
    helper function for rec2cds
    '''
    ndex = np.ceil(np.log10(max(abs((rcol)))))
    if np.isinf(ndex):
        ndex = 0 # all elements are zero 
    else: 
        ndex = int(ndex)
        
    sign_count = (min(rcol)<0)*1
    if np.issubdtype((rcol[0]), int) or type(rcol[0])==bool: 
        n_dig = max(1, ndex)+sign_count
        formatter = str(n_dig)
    else: 
        # for very large or small number use scintific notation 
        if (-ndex > (prec-3)) or (ndex>prec): 
            n_dig = prec+prec+sign_count+1
            formatter = '.'+str(int(prec))+'e'
        else: 
            n_dig = max(1,ndex)+sign_count+1+prec
            formatter = str(n_dig)+'.'+str(int(prec))+'f'

    ss='{0:'+formatter+'}' # more save: length of output string
    return formatter, len(ss.format(rcol[0]))


def rec2cds(rec, filename='', precision=None, unit=None, explanation=None):
    '''
    >> rec2cds(rec, filename='', precision=None, unit=None, explanation=None):

    convert records array to something like a CDS (fixed width) table

    optional input: 
     - precision: precision for floats (type: dictionary with column names or float)
     - unit, units of columns (type: dictionary)
     - explanation: describtion of columns, max 16 char  (type: dictionary)
    '''
    
    ncols = len(rec.dtype)
    formatter = np.repeat('__.__', ncols)
    n_dig = np.repeat(0, ncols)

    # same for all
    if np.isscalar(precision): 
        print('using '+str(precision), ' digits for all columns')
        for i, col in enumerate(rec.dtype.names):
            formatter[i], n_dig[i] =  col_formatter(rec[col], prec=precision)
            #print(col, '\t', n_dig[i], formatter[i])

    # dict for each col
    elif type(precision) == dict:
        for i, col in enumerate(rec.dtype.names):
            try: 
                prec = precision[col]
            except KeyError: 
                prec=7 
            formatter[i], n_dig[i] =  col_formatter(rec[col], prec=prec)
            #print(col, '\t', n_dig[i], formatter[i])

    # fixed 7 four digit precision if no input given
    else:
        print('using 7 digits for as default precision for floats')
        for i, col in enumerate(rec.dtype.names):
            formatter[i], n_dig[i] =  col_formatter(rec[col], prec=7)
            #print col, '\t', n_dig[i], formatter[i]

    # write, to terminal or file
    if filename: 
        print('writing to', filename)
        f = open(filename, 'w')

    # make the header:
    nbit = 1
    for i, col in enumerate(rec.dtype.names):

        # bits
        ss='{0:4.0f}-{1:4.0f}  '.format(nbit, nbit+n_dig[i])

        # units
        if unit: 
            try:
                ss+='{0:9} '.format(unit[col])
            except KerError:
                 ss+='{0:9} '.format('')
        # label
        ss+='{0:15} '.format(col)

        # explanation
        if explanation: 
            try:
                ss+='{0:35}'.format(explanation[col])
            except KerError:
                 ss+='{0:35}'.format('')
        
        print(ss)
        if filename:
            f.writelines(ss+'\n')
        
        nbit+=n_dig[i]+1
            
    # print body
    if filename: 
        f.writelines('\n')

    print('')
    
    for l, r in enumerate(rec):         
        line = ''
        for i, col in enumerate(rec.dtype.names):
            ss = '{0:'+formatter[i]+'} '
            line+=ss.format(r[col])
        if filename: 
            f.writelines(line+'\n')
            if np.mod(l, len(rec)/10)==0:
                print(line)
        else: 
            if np.mod(l, len(rec)/10)==0:
                print(line)

    if filename: 
        f.close()

    return 
