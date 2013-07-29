'''
functions for latex
make nice tables and convert floats to strings
'''

import numpy as np


def f2latex(flt, nd=1, debug=False):
    '''
    >> string = f2latex(flt, nd=1, debug=False)

    convert 1.67e23 to $1.7 times 10^{23}$
    nd is number of digits after decimal point
    '''

    def _doform(f, nd, debug):

        f = float(f)

        if f != 0:
            base = int(round(np.log10(np.abs(f))))
        else:
            base=0
            nd-=1
        
        num = f/10**base
        if debug: print 'num, base', num,base
        if (num < 1):
            base-=1
            num = f/10**base
            if debug: print 'new num base', num, base

#        if (abs(base)==1):
#            nd+=1
#            num = f
        form = '$%(1)0.'+str(nd)+'f'

        # eg 0.023 or 230
        if abs(base) > 1: form=form+'\\times 10^{'+str(int(base))+'}$'

            
        # eg 2.3
        if (base ==0):
            form = '$%(1)0.'+str(nd)+'f'+'$'
            num = f
        if (base ==1):
            form = '$%(1)0.'+str(nd-1)+'f'+'$'
            num = f
        # eg 0.23
        if base == -1:
            form = '$%(1)0.'+str(nd+1)+'f'+'$'
            num = f

        return form % {'1':num}

    try:
        out = []
        for f in flt:
           out.append(_doform(f,nd, debug))
    except TypeError:
        out = _doform(flt,nd, debug)
    return out
    



def rec2table(rec,  latex_names=None, dt_names=None, units=None, ndeci=1, filename=None):
    '''
    turn np.recarray into latex table
    >> rec2table(rec, units=None, filename=None)

    calls table function with recarray columns names
    optional input; 
    -- names: overwrites the recarray names
    -- units: second row of latex table
    -- filename: for writing 
    -- ndeci: number digits after decimal point (default=1)
    '''
    
    if not dt_names: 
        dt_names = rec.dtype.names

    if not latex_names: 
        latex_names = rec.dtype.names
    else: 
        latex_names = cols
    
    cols = []
    for n in dt_names:
        cols.append(rec[n])

    return table(cols, names=col_names, units=units, filename=filename, ndeci=ndeci)
    
def table(cols, names=None, units=None, ndeci=1, \
          filename=None, top=None, bottom=None):
    '''
    latex tabular from columns

    >> table([[42,2e8], ['spam', 'rabbit']], names=['numbers', 'names'])
     \hline 
     numbers & names \\
     \hline \hline
     $42$ & spam\\
     $2.0\times 10^{8}$ & rabbit\\

    input:
     - cols = list of columns 
    option input:
     - names = list of names (ie, first row of table)
     - units = list of units of columns (ie, second  row of table header)
     - ndeci=list for digits after decimal points 
     - filename='./table.tex'
     - top=str , top header  (eg, '\begin{table} \n ... ')
     - botom=str, bottom lines (eg, 'end{table}')
    '''
    if names:
        if len(cols) != len(names):
            print 'error, len(cols) != len(names):', len(cols), len(names)
            return

    if np.isscalar(ndeci): 
        ndeci = np.repeat(ndeci, len(cols))

    if len(ndeci) != len(cols):
        print 'error, len(ndeci) != len(cols)'
        return
            
    if filename:
        f = open(filename, 'w')

    #if given, add this first
    if top:
        f.writelines(top)

    # make the header
    if names:
        ss = '\hline \n'
        ss += names[0]    
        for n in names[1:]:
            ss += ' & '+n
        ss += ' \\\\ \n'

    if names and units:
        if (units[0]):
            ss += '('+units[0]+')'
        else:
            ss +=''
        for u in units[1:]:
            if u:
                ss += ' & '+ '('+u+')'
            else:
                ss +=' & '
        ss += ' \\\\ \n'

    if names:
        if filename:
            f.writelines(ss)
        else:
            print ss
        if filename:
            f.writelines('\\hline \\hline \n')
        else:
            print '\\hline \\hline'

    # make the body
    for i in range(len(cols[0])):

        cc = cols[0][i]
        if np.issubdtype(type(cc),str):
            ss = cc
        else:
           ss = f2latex(cc)
           
        for j in np.arange(1,len(cols)):
            cc = cols[j][i]
            if np.issubdtype(type(cc), str):
                ss +=' & '+ cc
            elif np.issubdtype(type(cc), int):
                 ss += ' & '+str(cc)            
            else:
                ss += ' & '+f2latex(cc, nd=ndeci[i])
                            
        ss += '\\\\'
        if filename:
            f.writelines(ss+'\n')
        else:
            print ss

    if bottom:
        f.writelines(bottom)
        
    return 
