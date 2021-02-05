import numpy as np
import pyccl as ccl
from .calculator_base import N5KCalculatorBase

from scipy.integrate import quad, nquad
from scipy.interpolate import interp1d

from pygsl.testing import sf

import time

import cppimport
funcs = cppimport.imp ("wrap_psico.wrap")

class N5KCalculatorPSICo(N5KCalculatorBase):
	name = 'PSICo'

	def BubbleSort_2vec (self, x, y):
		arr1inds = x.argsort()
		return x[arr1inds], y[arr1inds]  #ASCENDING ORDER

	def HH (self, z):
		return self.H0 * np.sqrt (self.Om * (1.+z)**3 + (1. - self.Om))

	def inv_HH (self, z):
		return  1 / np.sqrt (self.Om * (1.+z)**3 + (1. - self.Om))

	def Chiz (self, z):    
		f = lambda x : self.inv_HH (x)
		integ = quad (f, 0., z)
		return self.c * self.inv_H0 * integ[0]

	def hyp2f1 (self, a,b,c,x_hyp):
		prefac1 = sf.gamma (b-a) * sf.gamma (c) * pow (-x_hyp, -a) / sf.gamma (b) / sf.gamma (c-a);
		prefac2 = sf.gamma (a-b) * sf.gamma (c) * pow (-x_hyp, -b) / sf.gamma (a) / sf.gamma (c-b);
		inv_x_hyp = 1. / x_hyp
		hyp =  prefac1 * sf.hyperg_2F1 (a, a-c+1, a-b+1, inv_x_hyp) \
	    	+  prefac2 * sf.hyperg_2F1 (b, b-c+1, b-a+1, inv_x_hyp);
		return hyp

	def Dz (self, z):
    
		a_Dz = 1./3.
		b_Dz = 1.
		c_Dz = 11./6.
		x_Dz = (self.Om - 1.) / self.Om
    	
		if z == 0.:
			return 1.
        
		if np.abs (x_Dz) < 1.:
			den = sf.hyperg_2F1 (a_Dz, b_Dz, c_Dz, x_Dz);
		else: 
			den = self.hyp2f1 (a_Dz, b_Dz, c_Dz, x_Dz);
    	
		x1 = x_Dz / (1. + z)**3
		if np.abs (x1) < 1:
			num = sf.hyperg_2F1 (a_Dz, b_Dz, c_Dz, x1);
		else:
			num = self.hyp2f1 (a_Dz, b_Dz, c_Dz, x1);
	    	
		return num / den / (1. + z)
	
	def DChi (self, x):
		u = 2. * x * self.inv_Chimax - 1.
		dchi = 0.
		for i in range (self.LD):
			dchi = dchi + self.cD [i] * sf.legendre_Pl (i, u);
		return dchi
	
	def R (self, chi, l):
		arg = (l+0.5)/chi
		return self.Wgalaxy_noRSD (chi)**2 * self.Plin (arg) / chi**2
	
	def Cl_gal_limber (self, l):
		
		f = lambda x: self.R(x, l)
		clintg = quad (f, self.chimin_photoz, self.chimax_photoz)[0]
		
		return clintg / self.Dz (self.zav) / self.Dz (self.zav)
		
	def psicopy (self):
	
		#Parameters file
		self.c = 3e5
		h = 1.
		Omch2 = 0.03
		Ombh2 = 0.22
		ns = 0.95
		As = 2.21536e-9
	
		zCMB = 1100.
		zmax = 1500.
		self.bias = -1.98
	
		size_w = 1000
	
		self.zav = 0.5
		zm = 0.03
		sigmaz = 0.03
	
		lmin = 12
		lmax = 320
	
		#End of Parameters file
	
		self.H0 = 100. * h
		self.inv_H0 = 1. / self.H0
		Omc = Omch2 / (h**2)
		Omb = Ombh2 / (h**2)
		self.Om = Omc + Omb
	
		zbin_min = self.zav - 0.5 * zm * (1+self.zav)
		zbin_max = self.zav + 0.5 * zm * (1+self.zav)
	
		zphotoz_min = self.zav -  5. * sigmaz
		zphotoz_max = self.zav +  5. * sigmaz
	
		zarr = np.append (np.arange (0., 100., 0.1), np.arange (100., zmax+1, 1.))
	
		self.larr = self.get_ells ()
	
		Chiarr = np.array ([self.Chiz (z) for z in zarr]) 
		ChiCMB = self.Chiz (zCMB)
		Chimax = self.Chiz (zmax)
		self.inv_Chimax = 1. / Chimax
		zChi = interp1d (Chiarr, zarr)
	
		GLQ20 = np.loadtxt ("input_PSICo/GaussianQuadratureWeightsAndAbscissae_lmax20.dat")
		self.LD = 20
		uD = np.array  ([GLQ20[i][2] for i in range (self.LD)])
		wD = np.array  ([GLQ20[i][1] for i in range (self.LD)])
	
		uD, wD = self.BubbleSort_2vec (uD, wD)
	
		GLQ50 = np.loadtxt ("input_PSICo/GaussianQuadratureWeightsAndAbscissae_lmax50.dat")
	#Window Function GLQ
		GLChi = 50
		ChiChi = np.array  ([GLQ50[i][2] for i in range (GLChi)])
		wChi = np.array  ([GLQ50[i][1] for i in range (GLChi)])
	#ChiChi, wChi = self.BubbleSort_2vec (ChiChi, wChi)          #ATTENZIONE CHE, SE USATO NEL CODICE C++, 
	                                                        	#IL FATTO CHE NON VENGA ORDINATO POTREBBE ESSERE UN BUG
	    	
		self.cD = np.empty (self.LD+1)
		dataChi = np.empty (self.LD)
		for i in range (self.LD):
			dataChi [i] = self.Dz (zChi (Chimax * 0.5 * (uD [i] + 1.)))
		for l in range (self.LD+1):
			sumD = 0
			for i in range (self.LD):
				sumD = sumD + wD [i] * dataChi [i] * sf.legendre_Pl (l, uD [i])
			self.cD [l] = (2.*l+1.) * 0.5 * sumD
    	
		Wg1 = np.loadtxt ("input_PSICo/MICE_radial_selection.dat")
		self.Wg1_0 = np.array ([Wg1[i][0] for i in range (len(Wg1))])
		self.Wg1_1 = np.array ([Wg1[i][1] for i in range (len(Wg1))])
	
		self.Wg1_1 = self.HH (self.Wg1_0)/self.c * self.Wg1_1
		self.Wg1_0 = np.array ([self.Chiz (i) for i in self.Wg1_0])
		DChi_Wg1   = np.array ([self.DChi (i) for i in self.Wg1_0])
		#self.Wg1_1 = self.Wg1_1 * DChi_Wg1																#QUI
		self.Wg1_1 = self.Wg1_1
	
		self.Wgalaxy_noRSD = interp1d (self.Wg1_0, self.Wg1_1)
	
		Pk_lin_in = np.loadtxt ("input_PSICo/MICE_official_PS_cmb-conention_z0p5.dat")
		k = np.array ([Pk_lin_in [i][0] for i in range (len(Pk_lin_in))])
		Pk_lin = np.array ([Pk_lin_in [i][1] for i in range (len(Pk_lin_in))])
	
		klow = []
		Plow = []
	
	#print (k[0])
	#print (np.log (k[0])/np.log (10)-0.01) # == -3
	
		nklow = 0
		for kls in np.arange(-10., -3., 0.01):
			ten2k = np.power (10,kls)
			klow.append ([ten2k])
			Plow.append ([Pk_lin [0] * np.power (ten2k/k[0],ns)])
			nklow += 1
					     	
		self.karr = np.append (np.array(klow), k)
		self.Parr = np.append (np.array(Plow), Pk_lin)
	
		self.Plin = interp1d (self.karr, self.Parr)
	
		nk = len (self.karr)
	
		lims = np.where (self.karr**3 * self.Parr * np.power (self.karr,self.bias) > 100)
		limsx = lims[0][-1]
		limdx = lims[0][0]
		nksx = lims[0][-1]
		nkdx = lims[0][0]
	
		karrsx = []
		Parrsx = []
		indx=0
		for kk in self.karr [0:nkdx+1]:
			karrsx.append (kk)
			Parrsx.append (kk**3 * self.Parr [indx] * np.power (kk,self.bias))
			indx += 1
			    	
		karrdx = []
		Parrdx = []
		indx = 0
		for kk in self.karr [nksx:]:
			karrdx.append (kk)
			Parrdx.append (kk**3 * self.Parr [nksx + indx] * np.power (kk,self.bias))
			indx += 1
	
		karrsx = np.array (karrsx)
		Parrsx = np.array (Parrsx)
		karrdx = np.array (karrdx)
		Parrdx = np.array (Parrdx)
	
		Nint = 256
		N2 = int (Nint / 2)
		N = 256.
	    	
		kmin = 1e-5 #sarebbe min(k)?
		kmax = 1e3  #sarebbe max(k)?
	    	
		dkapa = np.log (kmax/kmin)
		dkapa_over_N = dkapa / N
	    	
		k0 = kmin
		chi0 = dkapa
	    	
		num = np.cdouble ()
	    	
		expdN = np.exp (dkapa_over_N)
		expdNm1 = expdN - 1.
	
		self.chimin_photoz = self.Chiz  (zphotoz_min)
		self.chimax_photoz = self.Chiz  (zphotoz_max)
	
		n = np.arange (Nint)
		chi = chi0 * np.exp (n * dkapa_over_N)
		k = k0 * np.exp (n * dkapa_over_N)
	
		kleft = 5e-5
		kright = 5e2
	
		in0 = k**3 * self.Plin (k) * np.power (k/k0, self.bias)
	
		left = (np.where (k < kleft))
		left_fact  = (k[left[0]] - kmin) / (kleft - kmin)
	
		right = (np.where (k > kright))
		right_fact =(kmax - k[right[0]]) / (kmax - kright) 
	
		Pin1 = np.empty(np.shape (in0))
	#in1 = np.empty(np.shape (in0), dtype='float64')
	
		Pin1 [left[0]] = in0 [left[0]] * (left_fact - 0.5 * np.sin (2. * np.pi * left_fact) / np.pi)
		Pin1 [left[0][-1]:right[0][0]] = in0 [left[0][-1]:right[0][0]] 
		Pin1 [right[0]] = in0 [right[0]] * (right_fact - 0.5 * np.sin (2. * np.pi * right_fact) / np.pi)
	
		over_em15 = np.where (Pin1>1e-15) #1e-15 è il limite dei double, oltre l'arrotondamento numerico ha effetto
		under_em15 = np.where (Pin1<1e-15) #1e-15 è il limite dei double, oltre l'arrotondamento numerico ha effetto
		Pin1 [under_em15] = 0.
	
		return 1.
	
##Limber with Python 100times slower than C++
#t0 =  time.time ()
#Clarr = np.empty (np.shape (self.larr))
#n = 0
#for l in self.larr:
#Clarr [n] = self.Cl_gal_limber (l)
#n += 1
#t1 =  time.time ()
#print ("Cl Limber ran in {} seconds".format (t1-t0))

	def setup(self):
		# Initialize cosmology
		check = self.psicopy ()
		if check != 1:
			print ("ERRORE CI FU!")


	def run (self):
		# Non-Limber
#print (self.larr)
		Clarr = funcs.Cl_NonLimber (self.larr.transpose (), 
																self.Wg1_0.transpose (), 
																self.Wg1_1.transpose (), 
																self.chimin_photoz, self.chimax_photoz, self.bias, 
																self.karr.transpose (),
																self.Parr.transpose ())
		np.savetxt ("test_n5k_ell.txt", self.larr)
		np.savetxt ("test_n5k_Cl.txt", Clarr.reshape (-1))

		# Limber gsl qags wrap
		ClarrLimb = funcs.Cl_Limber ( self.larr.transpose (), 
																	self.Wg1_0.transpose (), 
																	self.Wg1_1.transpose (), 
																	self.chimin_photoz, self.chimax_photoz, 
																	self.karr.transpose (),
																	self.Parr.transpose ())

		#ClarrLimb = ClarrLimb.reshape (-1) / self.Dz (self.zav) / self.Dz (self.zav) 					#QUI	
		ClarrLimb = ClarrLimb.reshape (-1)
		np.savetxt ("test_n5k_ell_2.txt", self.larr)
		np.savetxt ("test_n5k_ClLimb_gsl.txt", ClarrLimb)

#
#		# Limber scipy
#		ClarrLimb_scipy = np.empty (np.shape (self.larr))
#		n = 0
#		for l in self.larr:
#			ClarrLimb_scipy [n] = self.Cl_gal_limber (l)
#			n += 1
#		np.savetxt ("test_n5k_ell_3.txt", self.larr)
#		np.savetxt ("test_n5k_ClLimb_scipy.txt", ClarrLimb_scipy)
