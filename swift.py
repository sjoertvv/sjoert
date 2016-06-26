import datetime, calendar
from time import datetoyear, datetomjd


wdict = dict([('V',5468.),('B',4392.),('U',3465.),('UVW1',2600.),('UVM2',2246.),('UVW2',1928.)])
met0 = datetime.datetime(2001,1,1, 0,0,0)
met0_unix = calendar.timegm(met0.utctimetuple())#time.mktime(met0.timetuple())

def swift2mjd(tl):
	'''
	time since mission to mjd
	'''
	dd = datetime.datetime.utcfromtimestamp(met0_unix+tl)
	return datetomjd(dd)
def swift2year(tl):
	'''
	time since mission to mjd
	'''
	dd = datetime.datetime.utcfromtimestamp(met0_unix+tl)
	return datetoyear(dd)



