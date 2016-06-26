'''
functions from starutil_numpy (Dustin Lang),
plus stackoverflow
'''
import datetime
from datetime import datetime as dt
import time

import numpy as np

def jdtomjd(jd):
        return jd - 2400000.5
def mjdtojd(mjd):
        return mjd + 2400000.5
def mjdtodate(mjd):
        jd = mjdtojd(mjd)
        return jdtodate(jd)
def jdtodate(jd):
        unixtime = (jd - 2440587.5) * 86400. # in seconds
        return datetime.datetime.utcfromtimestamp(unixtime)
def mjdtoyear(mjd):
    if np.isscalar(mjd):
        return datetoyear(mjdtodate(mjd))
    else:
        out = np.zeros(len(mjd))
        for i in range(len(mjd)):
            out[i] = datetoyear(mjdtodate(mjd[i]))
        return out
def timedeltatodays(dt):
        return dt.days + (dt.seconds + dt.microseconds/1e6)/86400.
def datetomjd(d):
        d0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
        dt = d - d0
        # dt is a timedelta object.
        return timedeltatodays(dt)


def datetoyear(date):
    '''
    import is datetime object
    from http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    '''

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

