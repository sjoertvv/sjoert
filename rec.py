'''
functions for reccord arrays: slice, merge, append, initiate, etc. 
also some other tools for arrays/lists
2010  - started (Sjoert van Velzen)
'''
import numpy as np

def list2rec(dd, n=1, default=False):
    '''
    return empty records arrray from a list/tuple of pairs: (('col1', 0.),('col2', 'b') )
    >> rec = dict2rec(dict, n=10, default=False)

    n: number of rows
    set default=True to use value of dict as intial value (not empty)
    '''

    pre_dtype = []
    
    for el in dd:
        this_type  = np.array(el[1]).dtype
        print el[1], this_type
        pre_dtype+= [(el[0], this_type)]
    newrec = np.empty(n ,dtype=np.dtype(pre_dtype))

    if default:
        for el in dd:
            newrec[el[0]] = el[1]

    return newrec

def dict2rec(dd, n=1, default=False):
    '''
    return empty records arrray from dictionary keys
    >> rec = dict2rec(dict, n=10, default=False)

    n: number of rows
    set default=True to use value of dict as intial value (not empty)
    '''

    keys = dd.keys()
    pre_dtype = []
    for k in keys:
        this_type  = np.array(dd[k]).dtype
        print k, this_type
        pre_dtype+= [(k, this_type)]
    newrec = np.empty(n ,dtype=np.dtype(pre_dtype))

    if default:
        for k in keys:
            newrec[k] = dd[k]

    return newrec

def slice_cols(ori_rec, sub_cols):
    '''
    return new records array with subset cols from ori_rec
    >> new_rec = slice_cols(ori_rec, sub_cols)

    (there could be a more elegant way to do this,
    but this work fast enough for me)
    '''

    full_dtype = ori_rec.dtype.descr
    new_dtype = []

    for n in sub_cols:
        ii = np.where(np.asarray(ori_rec.dtype.names) == n)[0]
        if not(len(ii)):
            print 'error: sub_col not found in input rec:', n
            raise(IndexError)
        new_dtype.append(full_dtype[ii])

    # make new array and copy data from original
    new_rec = np.empty(len(ori_rec), dtype=new_dtype)
    for n in sub_cols:
        new_rec[n] = ori_rec[n]
    return new_rec


def merge_rec(d1,d2):
    '''
    append data from d1 to d2.
    >> new_rec = merge_rec(d1,d2)

    only collumns from d1 are used
    '''
    
    if len(d1.shape)==0:
        s1 = 1
    else:
        s1 = d1.shape[0] 
    if len(d2.shape)==0:
        s2 = 1
    else: 
        s2 = d2.shape[0]
    
    newrec = np.empty(s1+s2, dtype=d1.dtype)
    for key in d1.dtype.names:
        try:
            d2[key]
            newrec[key] = np.concatenate((d1[key], d2[key]))
        except ValueError:
            print 'waring '+ key+ ' not found in d2, putting this column to zero in output'
            newrec[key] =0.
    return newrec

def append_field(rec, name, arr, dtype=None):
    '''
    add data in to collumn with new name
    >> new_rec = append_field(rec, name, arr, dtype=None)

    only one field at the time
    arr can be scalar or array with len(rec)
    '''
    arr = np.asarray(arr)
    
    if dtype is None:
        dtype = arr.dtype
    newdtype = np.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = np.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec



def append_rows(rec, rows):
    '''
    append rows to recorc array
    >> new_rec = append_rows(rec, rows)

    rows should be an array of tuple(s)
    the order has to match the order of the orinigal rec
    '''
    if (type(rows[0]) != tuple):
        print 'Waring. rows are not an array of tuple(s).'
        
    d2 = np.array(rows, dtype=rec.dtype.descr)
    return merge_rec(rec, d2)

def match(a,b):
    '''
    return indices of matches in a and b
   >>  m_a, m_b = match(a,b):
 
    '''
    a = np.array(a)
    b = np.array(b)
    az = zip(range(len(a)),a)
    bz = zip(range(len(b)),b)

    ma = []
    mb = []

    for i, item in enumerate(a):
        inb = [xx[0] for xx in bz if xx[1] == item]
        #print inb
        for ib in inb:
            mb.append(ib)
            ma.append(i)
            
    return np.array(ma),np.array(mb)


