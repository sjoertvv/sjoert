from matplotlib import pyplot as plt
import numpy as np
from math import pi
import sjoert
from scipy.special import gamma as Gamma
from numpy import log10, log, sqrt

testing = False # todo, moving testing of these functions to seperate script

# constants
sigma_T = 6.65e-25 # cm^2
c = 3e10 # cm/s
qe = 4.8e-10 # esu
me = 9.11e-28 # gram
k = 1.38e-16 # erg/K

# defauls
gamma_max, gamma_min = 1e4, 1.
N_arr = int(1e4)

def nu_s(B, gamma):
	return 3/2. * gamma**2 * qe * B / (2*pi*me*c) 
def gamma_e(B, nu):
	return sqrt(nu /(0.29*3/2. * qe * B / (2*pi*me*c) )) #0.29 to account for peak emitted power

def cooling(B, nu, Gamma_L=1):
	return 3/sigma_T * np.sqrt(2*np.pi*c*me*qe/(B**3*Gamma_L)) * nu**(-1/2.)

def gamma_e99(B, nu):
	return np.sqrt(nu/34e9/B)*100

def cooling_gamma(B, gamma):
	return 7.7e6 * (B)**-2 *(gamma/100.)**-1 # same as Ghisellini Eq 4.11

def cooling99(B, nu):
	gam = gamma_e99(B, nu)
	return 7.7e6 * (B)**-2 *(gam/100.)**-1 # same as Ghisellini Eq 4.11

def cooling_extra(B, nu):
	gam = gamma_e(B, nu)
	return 7.7e6 * (B)**-2 *(gam/100.)**-1 # same as Ghisellini Eq 4.11


if testing:
	print (cooling99(0.4, 16e9)/24/3600)
	print (cooling(0.4, 16e9)/24/3600)
	print (cooling_extra(0.4, 16e9)/24/3600)

	key = input()

z = 0.0205778 # 14li redshift
D = sjoert.stellar.lumdis(z, h=0.7)
def slab(B,nu, r, D=D, p=2.5, gamma_max=1e2, gamma_min=gamma_min, eps_e=1):
	
	#kap1 = 2.9e-18 * B**(4.25) * nu**(-3.25)
	#eps1 = 5.6e-20 * B**(3.75) * nu**(-0.75)

	kap1 = alpha_nu(B, nu, p, gamma_max, gamma_min, eps_e)
	eps1 = Ptot(B, nu, p, gamma_max, gamma_min, eps_e)

	return r**2 / (4*D**2) * eps/kap * ( 1-np.exp(-kap*r) )
	#return r**2 / (4*D**2) * eps1/kap1 * ( 1+np.exp(-kap1*r)*(r*kap1/6-1) )

# Larmor frequency Eq. 4.4.1
def nu_L(B):
	return qe*B/(2*pi*me*c) 

# R&L Eq. 6.36 
def eps(B, nu, p=2.5):	
	fe = 3**(p/2.) * (2.25/p**2.2 + 0.105)
	U_B = B**2/(8*pi) # magnetic field energy
	
	#return 3*sigma_T*c*K(B, p)*U_B / (16*pi**(1.5)*nu_L(B)) * (nu/nu_L(B))**(-(p-1)/2.)*fe  #Ghisellini course, Eq. 4.45, gives similar results

	omega = nu*(2*pi)

	return 2*pi*3**0.5 * qe**3 *K(B, p)/(2*pi * me * c**2 * (p+1)) * \
			(me * c * omega / (3 * qe * B))**(-(p-1)/2.) *\
			sina_avg(p) * Gamma(p/2+19/12.) * Gamma(p/4+1/12.)


def sina_avg(p):
	return pi**0.5 / 2. * Gamma((p+6)/4.) / Gamma((p+8)/4.)

def kappa(B, nu, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min):

	# Eq.  4.51 Ghisellini
	fk = 3**((p+1)/2.)*(1.8/p**0.7 + p**2/40.)
	return pi**0.5 * qe**2 * K(B, p, gamma_max, gamma_min) / (8*me*c) * nu_L(B)**((p+2)/2.) * nu**(-(p+4)/2.)* fk 

def electron_dist(p=2.5, gamma_max=gamma_max, gamma_min=gamma_min):
	return np.linspace(gamma_min, gamma_max)

# normalization of electorn distribution (in energy and gamma units)
def K(B, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min,eps_e=1.):
	return B**2/(8*pi)/eps_e * (-p+2)/(me*c**2) / (gamma_max**(-p+2)-gamma_min**(-p+2))

def Ke(B, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min, eps_e=1.):
	E_max, E_min = gamma_max*me*c**2, gamma_min*me*c**2
	return B**2/(8*pi)/eps_e * (-p+2) / (E_max**(-p+2)-E_min**(-p+2))

# electron energy density
def Ne(B, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min, eps_e=1.):
	E_max, E_min = gamma_max*me*c**2, gamma_min*me*c**2
	#return Ke(B, p, gamma_max, gamma_min, eps_e) * 1/(-p+1) * (E_max**(-p+1)-E_min**(-p+1))	
	return K(B, p, gamma_max, gamma_min, eps_e) * 1/(-p+1) * (gamma_max**(-p+1)-gamma_min**(-p+1))	

# Eq 6.31c
def Besselx(x, inf=1e9):
	chi = np.linspace(x, inf, 1e6)
	return x * np.trapz(scipy.special.kv(5/3., chi), chi)

# from Ghisellini
def Besselx_approx(x):
	return 4*pi/(3**0.5 * Gamma(1/3.)) *(x/2.)**(1/3.)*np.exp(-x)

# Eq. 6.18
def Pomega(B, x):
	return 3**0.5/(2*pi)*qe**3*B / (me*c**2) * Besselx_approx(x)

# Eq 6.21a
def Ptot(B, nu_in, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min, eps_e=1):
	#gamma_arr = np.logspace(log10(gamma_min), log10(gamma_max), N_arr)
	gamma_arr = np.linspace(gamma_min, gamma_max, N_arr)
	omega_arr = 3*gamma_arr**2*qe*B / (2*me*c) # Eq. 6.17c
	if np.isscalar(nu_in):
		nu_in = np.array([nu_in])
	
	norm = K(B, p, gamma_max, gamma_min, eps_e)* sina_avg(p)* 2*pi
	out = np.zeros(len(nu_in))
	for i, nu in enumerate(nu_in):
		x =  (nu*(2*pi)) / omega_arr 
		out[i] = np.trapz(Pomega(B, x)*gamma_arr**(-p), gamma_arr)
	return  norm*out


if testing:
	print ('R&L =1.99 :', Ptot(1, 1e9, p=1.99))
	print ('Giz p=1.99 :',eps(1, 1e9, p=1.99))
	print ('Heino p=2  :', 5.7e-19* 1/np.log(1e4)) # https://books.google.com/books?id=2RDkj8NipEkC&printsec=frontcover&source=gbs_ge_summary_r&cad=0#v=onepage&q=absorption&f=false
	print ('')
	print ('R&L =2.5    :', Ptot(1, 1e9, p=2.5))
	print ('Giz (p=2.5)  :', eps(1, 1e9, p=2.5))
	print ('Heino (p=2.5):',5.6e-20)

# R&L Eq. 6.52
def alpha_nu(B, nu_in, p=2.5, gamma_max=gamma_max, gamma_min=gamma_min, eps_e=1):	
	if np.isscalar(nu_in):
		nu_in = np.array([nu_in])

	out = np.zeros(len(nu_in))
	#gamma_arr = np.logspace(log10(gamma_min), log10(gamma_max), N_arr)
	gamma_arr = np.linspace(gamma_min, gamma_max, N_arr)
	E_arr = me*c**2 * gamma_arr
	omega_arr = 3*gamma_arr**2*qe*B / (2*me*c) # Eq. 6.17c

	norm = (p+2)*c**2 * sina_avg(p) * Ke(B, p, gamma_max, gamma_min, eps_e) * 2*pi
	for i, nu in enumerate(nu_in): 
		x =  (nu*(2*pi)) / omega_arr 
		out[i] = 1 / (8*pi*nu**2) *np.trapz(Pomega(B, x) *E_arr**(-p)/E_arr, E_arr) # without change of variables, gives same answer as next time 
			#Ke(B, p)*(me*c**2)**(-p) * np.trapz(Ptot(B, nu, p) *gamma_arr**(-p)/gamma_arr, gamma_arr) # after change of variables

	return norm*out


# # Ghisellini 4.448
# def kappa_alpha(B, nu_in, p=2.5):
# 	omega_arr = 3*gamma_arr**2*qe*B / (2*me*c) # Eq. 6.17c
# 	if np.isscalar(nu_in):
# 		nu_in = np.array([nu_in])
# 	out = np.zeros(len(nu_in))
# 	for i, nu in enumerate(nu_in):
		
# 		out[i] = 1/(8*pi*nu**2) Ke * np.trapz()
# 	return  out

# --extra--
# For Kino et al. https://arxiv.org/pdf/1403.0650.pdf
# this function doesn't have pitch angle averaging
def c1(p):
	return Gamma((3*p+2)/12.) * Gamma((3*p+22)/12.) # this is actually a very weak function of p
# Eq. 7
def alpha_nu_Kino(B, nu, p):
	return sqrt(3)*qe**3/(8*pi*me)*(3*qe/(2*pi*me**3*c**5)) * c1(p) \
			* B**(p+2) *Ke(p) * nu**(-(p+4)/2.)

if testing: 
	print ('')
	print ('Giz p=2.0  :',kappa(1, 1e9, 2.01))
	print ('R&L p=2.0  : ',alpha_nu(1, 1e9, 2.01))
	print ('Heino p=2  :', 4.5e-12 /log(gamma_max/gamma_min))


	nu = np.logspace(8.5, 11.5, 100)


	plt.clf()

	plt.plot(nu, slab(1.,nu, 1e15, p=2.01, gamma_max=1e3), '-', label='R&L')
	plt.plot(nu, slab(1.*2,nu, 1e15/2**2, p=2.01, gamma_max=1e3), '-', label='R&L')
	plt.plot(nu, slab(1.*4,nu, 1e15/4**2, p=2.01, gamma_max=1e3), '-', label='R&L')

	plt.legend(loc=0)
	plt.yscale('log')
	plt.xscale('log')

	plt.show()
	key = raw_input()



	plt.clf()
	plt.plot(nu, slab(1,nu, 1e15, p=2.01, gamma_max=1e2), '-', label='R&L, p=2.0, gamma_max=1e2')
	plt.plot(nu, slab(1,nu, 1e15, p=2.01, gamma_max=1e3), '-', label='R&L, p=2.0, gamma_max=1e3')
	plt.plot(nu, slab(1,nu, 1e15, p=2.01, gamma_max=1e4), '-', label='R&L, p=2.0, gamma_max=1e4')
	plt.plot(nu, slab(1,nu, 1e15, p=3.0, gamma_max=1e4), '-', label='R&L, p=3.0, gamma_max=1e4')

	plt.legend(loc=0)
	plt.yscale('log')
	plt.xscale('log')

	key = raw_input()



	plt.clf()
	plt.plot(nu, alpha_nu(1,nu, p=1.99), 'k-', label='R&L, p=2.0')
	plt.plot(nu, alpha_nu_Kino(1,nu, p=1.99),'k--', label='Kino, p=2.0')

	plt.plot(nu, alpha_nu(1,nu, p=3), 'r-', label='R&L, p=3')
	plt.plot(nu, alpha_nu_Kino(1,nu, p=3), 'r--', label='Kino, p=3')

	plt.legend(loc=0)
	plt.yscale('log')
	plt.xscale('log')

	key	=raw_input()

	#print (kappa(1, 1e9, p=2.5))
	#print (kappa(1, 1e9, p=3))

	p = 2.5
	# print (3**((p+1)/2.)*(1.8/p**0.7 + p**2/40.))
	# print (3**((p+1)/2.)*Gamma((3*p+22)/12.)*Gamma((3*3+2)/12.)*Gamma((p+6)/4.)/Gamma((p+8)/4.))

	# estimate B field from self-absorbed part and limit on size from EVN
	nu = 2e9
	Snu = 1e-3 # Jansky
	Inu = Snu*1e-23*4*sjoert.stellar.lumdis(z)**2/(1e17)**2 # estimate/lower limit
	Tb = Inu*c**2/(2*k*nu**2)
	# print ('Tb (/1e10)', Tb/1e10)
	# print ('B (G)', 1.4e12*nu*Tb**(-2))

