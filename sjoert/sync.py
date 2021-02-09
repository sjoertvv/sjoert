from matplotlib import pyplot as plt
import numpy as np
from math import pi
import sjoert
from scipy.special import gamma as Gamma
from numpy import log10, log, sqrt


'''
note, for all functions below eps_e= eps_e / eps_B

the formal/numerical equipartition is  ε_B /ε_e ≈ 6/11
'''


testing = False # todo, moving testing of these functions to seperate script

# constants
sigma_T = 6.65e-25 # cm^2
c = 3e10 # cm/s
qe = 4.8e-10 # esu
me = 9.11e-28 # gram
mp = 1.67e-24 # gram
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


if testing:
	print (cooling99(0.4, 16e9)/24/3600)
	print (cooling(0.4, 16e9)/24/3600)
	print (cooling_extra(0.4, 16e9)/24/3600)

	key = input()

def Beq(S_peak, D_L, nu_p, alpha=1, f=0.5,  who='Chevalier'):
	'''
	B_eq = Beq(S_peak, D_L, nu_p, z)

	Chevalier Eq. 14 (p=3), no z scaling
	'''
	return 0.58 * \
		alpha**(-4/19) * \
		(S_peak/1e3)**(-2/19.) * \
		(D_L/(1e6*sjoert.stellar.parsec))**(-4/19) * \
		(nu_p/5e9) * \
		(f/0.5)**(-4/19)


def Req(S_peak, D_L, nu_p, z, Gamma_bulk=1., 
		fA=1, fV=4/3, f=1, 
		p=3,
		epsilon_e=None,
		epsilon_B=None,
		alpha=6/11,
		verbose=False, who='Barniol'):
	'''
	R_eq = Req(S_peak, D_L, nu_p, z, fA=1, fV=4/3)

	Barniol-Duran+13 Eq. 21

	S_peak in mJy
	D_L = cm
	nu_p in GHz (observed)
	Gamma_bulk=1. is the bulk Lorentz factor
	
	alpha=6/11 	eps_B /eps_e is "true equipartition" (only used in Chevalier for now, 
				note the paper seems to have the inverse below Eq. 10, but that must be wrong)
	
	f=1  		filling factor of Chevalier, such that Volume=4/3 *f pi*R**3
	
	fV=4/3
	fA=1		the volume and area factors of Barniol-Duran (2013), fV=4/3 is the sperical Newtonian case

	epsilon_e=None 	the fraction of the proton energy that goes into electrons 
				 	this is used to estimate gamma_min, but only works when Gamma>~few
				 	gamma_min=2 is enforced by default
	epsilon_B=None  the fraction of total energy that goes into the magnetic field

	return radius in cm
	'''

	# Eq. 27
	if who=='Barniol':
		
		if epsilon_e and Gamma_bulk>1.1:
			xhi_e = (p-2)/(p-1) * epsilon_e * mp/me
			xhi_e = np.clip(xhi_e, 2/(Gamma_bulk-1))
		else:
			xhi_e = 2


		out= 1e17 * \
			(21.8 * (525)**(p-1))**(1/(13+2*p)) * \
			(xhi_e**((2-p)/(13+2*p))) * \
			(S_peak**((6+p)/(13+2*p))) * \
			((D_L/1e28)**(2*(p+6)/(13+2*p))) * \
			((nu_p/1e10)**(-1)) * \
			((1+z)**(-(19+3*p)/(13+2*p))) * \
			(fA**(-(5+p)/(13+2*p))) * \
			(fV**(-1/(13+2*p))) * \
			(Gamma_bulk**((p+8)/(13+2*p)))
			
		# increase by accounting for hot protons (section 4.2.2. of Barniol)
		if epsilon_e:
	
			out*= (1+1/epsilon_e)**(1/(13+2*p)) 

			if Gamma_bulk>1.1:
				out*=((Gamma_bulk-1)**((2-p)/(13+2*p)))
		
		# do the correction for out-of-equiparition systems
		if epsilon_B:
			eta = (epsilon_B/ (1-epsilon_B)) / (6/11) # note this this different from Barniol because this account for possibility of hot protons
		else:
			eta = 1

		if Gamma_bulk<1.1:
			out *= eta**(1/17)
		else:
			out *= eta**(1/12)

		# geomtrical correction for Newonian case
		if Gamma_bulk<1.1:
			out*=(4**(1/(13+2*p)))

		return out
		

	if who=='Alexander':
		return 4**(1/19) * \
		 3.2e15 * 	(S_peak)**(9/19) * \
					(D_L/1e26)**(18/19) * \
					(nu_p/1e10)**(-1) * \
					(1+z)**(-10/19) * \
					fA**(-8/19) *\
					fV**(-1/19)

	# Eq. 13 (p=3, spherical non-relativistic), no z depent from Alexander
	if who =='Chevalier':
		return 8.8e15 * \
			alpha**(-1/19) * \
			(S_peak/1e3)**(9/19) * \
			(D_L / (1e6*sjoert.stellar.parsec))**(18/19) * \
			(nu_p / 5e9)**-1 * \
			(f/0.5)**(-1/19) * \
			(1+z)**(-10/19) 

def Eeq(S_peak, D_L, nu_p, z, Gamma_bulk=1.,fA=1, fV=4/3, f=1, 
		p=3,
		epsilon_e=None,
		epsilon_B=None,
		alpha=6/11,
		verbose=False, who='Barniol'):
	'''
	E_eq = Req(S_peak, D_L, nu_p, z, fA=1, fV=4/3)

	Barniol-Duran+13 Eq. 25

	S_peak in mJy
	D_L = cm
	nu_p in GHz (observed)
	Gamma_bulk=1 is the bulk Lorentz factor
	
	f=1  filling factor of Chevalier, such that Volume=4/3 *f pi*R**3

	fV=4/3
	fA=1		the volume and area factors of Barniol-Duran (2013), fV=4/3 is the sperical Newtonian case

	alpha=6/11 	eps_B /eps_e is "true equipartition" (only used in Chevalier for now, 
				note the paper seems to have the inverse below Eq. 10, but that must be wrong)

	epsilon_e=None 	the fraction of the proton energy that goes into electrons 
				 	this can also be used to estimate gamma_min, but only works when Gamma>few
				 	gamma_min=2 is enforced by default
	epsilon_B=None  the fraction of total energy carried by the magnetic field (if None we assume system in equipartition)



	return energy in erg
	'''

	# Eq. 28
	if who=='Barniol':
		
		if epsilon_e and Gamma_bulk>1.1:
			xhi_e = (p-2)/(p-1) * epsilon_e * mp/me
			xhi_e = np.clip(xhi_e, 2/(Gamma_bulk-1))
		else:
			xhi_e = 2

		out= 1.3e48 * \
			((21.8)**(-2*(p+1)/(13+2*p))) * \
			(( 525**(p-1) * xhi_e**(2-p) )**(11/(13+2*p))) * \
			(S_peak**((14+3*p)/(13+2*p))) * \
			((D_L/1e28)**(2*(14+3*p)/(13+2*p))) * \
			((nu_p/1e10)**(-1)) * \
			((1+z)**(-(27+5*p)/(13+2*p))) * \
			(fA**(-3*(p+1)/(13+2*p))) * \
			(fV**(2*(p+1)/(13+2*p))) * \
			(Gamma_bulk**(-(5*p+16)/(13+2*p))) 

		# increase by accounting for hot protons
		if epsilon_e:
			out*= (1+1/epsilon_e)**(11/(13+2*p)) 
			if Gamma_bulk>1.1:
				out*=((Gamma_bulk-1)**(-11*(p-2)/(13+2*p)))

		# do the correction of radius for out-of-equiparition systems
		if epsilon_B:
			eta = (epsilon_B/ (1-epsilon_B)) / (6/11)
		else:
			eta = 1

		if Gamma_bulk<1.1:
			out *= (11/17)*eta**(-6/17) + (6/17) * eta**(11/17)
		else:
			out *= (11/17)*eta**(-5/12) + (6/17) * eta**(7/12)


		# geomtrical correction for Newonian case
		if Gamma_bulk<1.1:
			out*=(4**(11/(13+2*p)))


		return out

	# Eq. 25, give an absolute lower limit below energy from other electorn not accounted for
	if who=='Barniol_at_nu_p':
		out= 2.5e49 * (S_peak)**(20/17) * \
						(D_L/1e28)**(40/17) * \
						(nu_p/1e10)**(-1) * \
						(1+z)**(-37/17) * \
						Gamma**(-26/17) *\
						fA**(-9/17) *\
						fV**(6/17) 
		if Gamma_bulk<1.1:
			 out*=4**(11/17)
		return out

	if who=='Alexander':
		
		return 4**(11/19) * 1.9e46 * (S_peak)**(23/19) * \
							(D_L/1e26)**(46/19) * \
							(nu_p/1e10)**(-1) * \
							(1+z)**(-42/19) * \
							fA**(-12/19) *\
							fV**(8/19) 

	# adding z-dependence from Alexander
	if who=='Chevalier':
		#eta = 1/alpha/(6/11)
		return  (1+z)**(-42/19) * \
				4/3 * (f/0.5) * Req(S_peak, D_L, nu_p, z*0, f=f, alpha=alpha, who='Chevalier')**3 * \
				 Beq(S_peak, D_L, nu_p, f=f, alpha=alpha, who='Chevalier')**2 / (8*pi)


def slab(B,nu, r, D=3e27, p=2.5, gamma_max=1e3, gamma_min=1, eps_e=1):
	'''
	Snu = slab(B,nu, r, D=D, p=2.5, gamma_max=1e3, gamma_min=1, eps_e=1):

	the equation below works because Ptot is not per unit solid angle
	eps = Ptot/4pi
	and S_nu = Omega*I_nu = pi*(d/D)**2 I_nu 
	'''
	#kap1 = 2.9e-18 * B**(4.25) * nu**(-3.25)
	#eps1 = 5.6e-20 * B**(3.75) * nu**(-0.75)

	kap1 = alpha_nu(B, nu, p, gamma_max, gamma_min, eps_e)
	eps1 = Ptot(B, nu, p, gamma_max, gamma_min, eps_e)/(4*pi)

	d = r*2
	return np.pi * r**2 / (D**2) * eps1/kap1 * ( 1-np.exp(-kap1*d) )
	#return np.pi * r**2 / (D**2) * eps1/kap1 * ( 1+np.exp(-kap1*d)*(d*kap1/6-1) ) # from Flacke99

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
	'''
	radiated power per unit frequency per unit volume per unit time
	'''
	#gamma_arr = np.logspace(log10(gamma_min), log10(gamma_max), N_arr)
	gamma_arr = np.linspace(gamma_min, gamma_max, N_arr)
	omega_arr = 3*gamma_arr**2*qe*B / (2*me*c) # Eq. 6.17c
	if np.isscalar(nu_in):
		nu_in = np.array([nu_in])
	
	norm = K(B, p, gamma_max, gamma_min, eps_e)* sina_avg(p) * 2*pi
	out = np.zeros(len(nu_in))
	for i, nu in enumerate(nu_in):
		x =  (nu*(2*pi)) / omega_arr 
		out[i] = np.trapz(Pomega(B, x)*gamma_arr**(-p), gamma_arr)
	return  norm*out


if testing:
	print ('emission coeff or total power:')
	print ('R&L (interp) p=2 :', Ptot(1, 1e9, p=2.01, gamma_max=1e4))
	print ('R&L (gamma) p=2  :', eps(1, 1e9, p=2.01))
	print ('VLBI Book   p=2  :', 5.7e-19* 1/np.log(1e4)) # https://books.google.com/books?id=2RDkj8NipEkC&printsec=frontcover&source=gbs_ge_summary_r&cad=0#v=onepage&q=absorption&f=false, page 22
	print ('')
	print ('R&L (interp) p=2.5    :', Ptot(1, 1e9, p=2.5))
	print ('R&L (gamma) (p=2.5)  :', eps(1, 1e9, p=2.5))
	print ('Heino (p=2.5):',5.6e-20)
	print ('')
	print ('absorption:')
	print ('R&L (interp) p=2    :', alpha_nu(1, 1e9, p=2.01, gamma_max=1e4))
	print ('VLBI Book (p=2.0)	:',4.5e-12/np.log(1e4))


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
		out[i] = 1 / (8*pi*nu**2) *np.trapz(Pomega(B, x) *E_arr**(-p)/E_arr, E_arr) # without change of variables, gives same answer as next line 
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

	key = input()



	plt.clf()
	plt.plot(nu, alpha_nu(1,nu, p=1.99), 'k-', label='R&L, p=2.0')
	plt.plot(nu, alpha_nu_Kino(1,nu, p=1.99),'k--', label='Kino, p=2.0')

	plt.plot(nu, alpha_nu(1,nu, p=3), 'r-', label='R&L, p=3')
	plt.plot(nu, alpha_nu_Kino(1,nu, p=3), 'r--', label='Kino, p=3')

	plt.legend(loc=0)
	plt.yscale('log')
	plt.xscale('log')

	key	=input()

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

