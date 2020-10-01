

from math import radians

import numpy as np

catsHTM_path = '/Volumes/LaCie/catsHTM/'

def sdss_to_gaia(g,i, DR=2):
    gi = g-i

    if DR==1:
        return g -0.0912  -0.5310*gi -0.1042*gi**2   +0.0068*gi**3 # Jordi+10, DR1)
    else:
        return g -0.074189  -0.51409*gi -0.080607*gi**2   +0.0016001*gi**3 #(Evans+15 DR2)

def get_gaia(ra, dec, dist, catsHTM_path=catsHTM_path, verbose=False):
    '''
    >>> recarr, distarr = get_gaia(ra, dec, dist, catsHTM_path='/somewhere/')
    function to run catsHTM.cone_search and convert to np.array
     ra, dec in deg.
     dist in arcsec
    '''
    try:
        import catsHTM 
    except ImportError:
        print ('import catsHTM failed')
        print ('try: \n')
        print ('pip install catsHTM')
        
    # some files are corrupted, we have to catch the exception
    try:
        srcs, colnames, colunits = catsHTM.cone_search(
                                        'GAIADR2',
                                        np.radians(ra),
                                        np.radians(dec),
                                        dist,
                                        catalogs_dir=catsHTM_path)
        if len(srcs) == 0:
            if verbose:
                print('''no Gaia sources within {0:0.2f}"'''.format(dist))
            return None, None
    
    except OSError as ose:
        print("OSError from catsHTM %s"%str(ose))
        return None, None

    # convert output to numpy array
    gaia_arr = np.zeros(len(srcs[:,0]), dtype=[(k,'f8') for k in colnames])
    for i,k in enumerate(gaia_arr.dtype.names):
        gaia_arr[k] = srcs[:,i]

    gaia_match_dist = 3600*np.degrees(catsHTM.celestial.sphere_distance_fast(gaia_arr["RA"], gaia_arr["Dec"],
                                                                 np.radians(ra), np.radians(dec)))
    return gaia_arr, gaia_match_dist




