import datetime, calendar
from simtime import datetoyear, datetomjd


wdict = dict([('V',5468.),('B',4392.),('U',3465.),('UVW1',2600.),('UVM2',2246.),('UVW2',1928.)])
met0 = datetime.datetime(2001,1,1, 0,0,0)
met0_unix = calendar.timegm(met0.utctimetuple())#time.mktime(met0.timetuple())

vega2ab_dict = dict([('V',-0.01),('B',-0.12),('U',1.02),('UVW1',1.53),('UVM2', 1.68),('UVW2',1.69)]) # based on ASASSN-15lh photometry


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



