'''
make html tables from arrays or strings and upload 
2011 - Sjoert van Velzen

example:

html_header('My page')
tab_header(['numbers', 'names'])
tab_row(['1', 'spam'])
tab_row(['2', 'rabbit'])
tab_end()
html_end()

'''

import numpy as np
import os

def html_header(title='', tracker=''):
    out= '<html>\n  <head> \n  <title>'+title+'</title>\n '+\
        tracker +'</head>\n  <body>\n'
    return out

def tab_header(sarr, width=[]):
    out = '  <table border=2 celpadding=5> \n'
    out += '<tr> \t'
    if (len(sarr) != len(width)):
        width = np.repeat('300', len(sarr))
    for i,s in enumerate(sarr):
        out+= '<th '+ 'width="'+width[i]+'"><u>'+s+'</u></th> \t'
    return out + '</tr> \n'
    
def tab_row(sarr):
    '''
    a row of html table
    >> tab_row(['spam', '12'])
    '''
    out = '   <tr> \t'
    for s in sarr:
        out+= '<td>'+s+'</td> \t'
    return out + '</tr> \n'

def tab_end():
    return ' </table>'
def html_end():
    return  '</body>\n</html>'

def upload(local_folder=None, remote_host='$HELIOS', remote_folder=None):
    '''
    >> upload(local_folder='./', remote_host='user@astro.edu, remote_folder=~/public_html')

    rsync local folder to server, some attempts are made to set proper permissions 
    '''
    
    if not(local_folder):
        print 'please give local folder to upload'
        return
    
    os.system('chmod 664 ' +local_folder+'*/*')
    os.system('chmod 664 ' +local_folder+'*.html')

    #  get the last bit of the local folder,
    #  this the base of the remote dir 
    if not(remote_folder):
        remf = local_folder.split('/')
        remote_folder = remf[len(remf)-2] +'/'        
    print 'copy from:', local_folder
    print 'to', remote_host+':'+remote_folder
    
    os.system('rsync -rv --progress '+local_folder+' '+ remote_host+':'+remote_folder)
    
    os.system('ssh ' + remote_host + ' "chmod  755 '+ \
          remote_folder+' "'  )
    os.system('ssh ' + remote_host + ' "chmod  644 '+ \
              remote_folder+'*.html "'  )
    os.system('ssh ' + remote_host + ' "chmod  644 '+ \
              remote_folder+'*.* "'  )
    os.system('ssh ' + remote_host + ' "chmod  644 '+ \
              remote_folder+'*/* "'  )
    
