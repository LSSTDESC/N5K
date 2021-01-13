import numpy as np
import pyccl as ccl
from .calculator_base import N5KCalculatorBase
import matterlib
from scipy.interpolate import CubicSpline as interp
import matplotlib.pyplot as plt

class N5KCalculatorMATTER(N5KCalculatorBase):
    name = 'matter'

    def setup(self):
        # Initialize cosmology
        print("Initializing lol")
        par = self.get_cosmological_parameters()
        self.cosmo = ccl.Cosmology(Omega_c=par['Omega_m']-par['Omega_b'],
                                   Omega_b=par['Omega_b'],
                                   h=par['h'], n_s=par['n_s'],
                                   A_s=par['A_s'], w0=par['w0'])
        dpk = self.get_pk()
        a = 1./(1+dpk['z'][::-1])
        self.cosmo._set_linear_power_from_arrays(a_array=a,
                                                 k_array=dpk['k'],
                                                 pk_array=dpk['pk_lin'][::-1][:])
        self.cosmo._set_nonlin_power_from_arrays(a_array=a,
                                                 k_array=dpk['k'],
                                                 pk_array=dpk['pk_lin'][::-1][:])#dpk['pk_nl'][::-1][:])
        # Initialize tracers
        if self.config.get('tracers_from_kernels', False):
            tpar = self.get_tracer_parameters()
            ker = self.get_tracer_kernels()
            a_g = 1./(1+ker['z_cl'][::-1])
            self.t_g = []
            for k, b in zip(ker['kernels_cl'], tpar['b_g']):
                t = ccl.Tracer()
                barr = np.full(len(a_g), b)
                t.add_tracer(self.cosmo,
                             (ker['chi_cl'], k),
                             transfer_a=(a_g, barr))
                self.t_g.append(t)
            self.t_s = []
            for k in ker['kernels_sh']:
                t = ccl.Tracer()
                t.add_tracer(self.cosmo,
                             kernel=(ker['chi_sh'], k),
                             der_bessel=-1, der_angles=2)
                self.t_s.append(t)
        else:
            nzs = self.get_tracer_dndzs()
            tpar = self.get_tracer_parameters()
            z_g = nzs['z_cl']
            z_s = nzs['z_sh']
            self.t_g = [ccl.NumberCountsTracer(self.cosmo, True, (z_g, n),
                                               bias=(z_g, np.full(len(z_g), b)))
                        for n, b in zip(nzs['dNdz_cl'], tpar['b_g'])]
            self.t_s = [ccl.WeakLensingTracer(self.cosmo, (z_s, n), True,
                                              (z_s, np.full(len(z_s), A_IA)))
                        for n, A_IA in zip(nzs['dNdz_sh'], tpar['A_IA'])]
        # The above should be (mostly) the same as in the calculator_ccl

        # 1) Find out what is the minimum and maximum relevant value of chi for each window
        Ntg = len(self.t_g)
        Nts = len(self.t_s)
        Nchi_test = 20000
        threshold = 1e-10
        age_test = ccl.comoving_radial_distance(self.cosmo,1e-4)*1.5 # Use the chi(z=10000.) as an upper bound of possible chi for the window functions. Factor 1.5 for safety
        chi_test = np.linspace(0,age_test,num=Nchi_test)
        self.chi_g_mins,self.chi_g_maxs = np.empty((2,Ntg),dtype="float64")
        self.chi_s_mins,self.chi_s_maxs = np.zeros((2,Nts,2),dtype="float64")
        for i,tg in enumerate(self.t_g):
           tg_test = tg.get_kernel(chi_test)
           maxtg = np.max(tg_test)
           mask = tg_test>threshold*maxtg
           self.chi_g_maxs[i] = chi_test[len(mask[0])-np.argmax(mask[0][::-1])-1]
           self.chi_g_mins[i] = chi_test[np.argmax(mask[0])]
        for i,ts in enumerate(self.t_s):
           ts_test = ts.get_kernel(chi_test)
           for itr in range(len(ts_test)):
             maxts = np.max(ts_test[itr])
             mask = ts_test[itr]>threshold*maxts
             self.chi_s_maxs[i][itr] = chi_test[len(mask)-np.argmax(mask[::-1])-1]
             imin = np.argmax(mask)
             self.chi_s_mins[i][itr] = (chi_test[imin] if imin>1 else 0.)
        # Now the minima and maxima in chi of the window functions have been found

        # 2) Define corresponding chi sampling
        Nchi_nonintegrated = 800
        Nchi_integrated = 1600
        self.chi_nonintegrated = [np.linspace(self.chi_g_mins[i],self.chi_g_maxs[i],num=Nchi_nonintegrated) for i in range(Ntg)]
        self.chi_integrated = [np.linspace(np.min(self.chi_s_mins[i]),np.max(self.chi_s_maxs[i]),num=Nchi_integrated) for i in range(Nts)]

        # 3) Get the Kernel and Transfer at this chi sampling
        self.kerfac_g = np.zeros((Ntg,Nchi_nonintegrated))
        for i in range(Ntg):
          kern = self.t_g[i].get_kernel(self.chi_nonintegrated[i])
          trans = self.t_g[i].get_transfer(0.,ccl.scale_factor_of_chi(self.cosmo,self.chi_nonintegrated[i]))
          growth_g = ccl.growth_rate(self.cosmo, ccl.scale_factor_of_chi(self.cosmo,self.chi_nonintegrated[i]))
          for itr in range(1):#range(kern.shape[0]):
            self.kerfac_g[i] += kern[itr]*trans[itr]
        self.kerfac_s = np.zeros((Nts,Nchi_integrated))
        for i in range(Nts):
          kern = self.t_s[i].get_kernel(self.chi_integrated[i])
          trans = self.t_s[i].get_transfer(0.,ccl.scale_factor_of_chi(self.cosmo,self.chi_integrated[i]))
          for itr in range(1):#range(kern.shape[0]):
            self.kerfac_s[i] += kern[itr]*trans[itr]
        
        # Shorthand notations:
        power = dpk['pk_lin'][::-1]
        Na_pk = len(power)

        # 4) Get the growth factor and pass it as well (sampled on same scale factor grid as the P(k))
        self.a_pk = a
        # self.growth = ccl.growth_rate(self.cosmo, self.a_pk) <-- This is not accurate enough. Instead, we replace it by pk_growth
        pk_growth = np.empty((Na_pk,))
        for i in range(Na_pk):
          pk_growth[i] = np.sqrt(np.mean(power[i]/power[-1]))
        self.growth = pk_growth
        growth_func = interp(self.a_pk,pk_growth)
        for i in range(Ntg):
          self.kerfac_g[i] *= growth_func(ccl.scale_factor_of_chi(self.cosmo,self.chi_nonintegrated[i]))
        for i in range(Nts):
          self.kerfac_s[i] *= growth_func(ccl.scale_factor_of_chi(self.cosmo,self.chi_integrated[i]))

        # 5) The kmin of the provided file is a bit too high. Here we extrapolate using k^(n_s) to reach lower k values (down to kmin = 1e-7/Mpc)
        kmin = 1e-7
        self.Nk_fft = 256
        Nk_small = int(np.log10(dpk['k'][0]/1e-7)/np.log10(dpk['k'][1]/dpk['k'][0])+1)
        assert(Nk_small > 10)
        ksmall = np.geomspace(1e-7,dpk['k'][0],endpoint=False,num=Nk_small)
        k_all = np.concatenate([ksmall,dpk['k']])
        self.tau_pk = ccl.comoving_radial_distance(self.cosmo,self.a_pk)
        self.k_pk = np.geomspace(kmin,dpk['k'][-1],num=self.Nk_fft)
        self.pk = np.empty((Na_pk,self.Nk_fft))
        self.deltaksq = np.empty((Na_pk,self.Nk_fft))
        for i in range(Na_pk):
          pk_all = np.concatenate([(ksmall/dpk['k'][0])**(par['n_s'])*power[i][0],power[i]])
          self.pk[i] = interp(k_all,pk_all)(self.k_pk)
          self.deltaksq[i] = interp(k_all,pk_all)(self.k_pk)*self.k_pk**3/(2.*np.pi**2)
        #Special note : we want Delta(k)=P(k)*k^3/(2*pi^2) instead of P(k) for the matterlib


        # 6) Pass everything to the matterlib
        self.ma = matterlib.Matter(ma_verbose=3)

        sfftcutoff = 100#50 (faster, but more inaccurate)
        stw = 50#25 (faster, but more inaccurate)
        sitw = 75
        self.ma.set((self.chi_nonintegrated,self.chi_integrated),(self.kerfac_g,self.kerfac_s),(self.a_pk,self.tau_pk,self.k_pk,self.deltaksq),self.growth,lmax=2000,size_fft_cutoff=sfftcutoff,tw_size=stw,integrated_tw_size=sitw,
        uses_separability=True)

    def run(self):
        # Compute power spectra
        self.ma.compute()
        
        # Get the Cl's for the given l array
        ls = self.get_ells()
        cls = self.ma.matter_cl(ls)

        # Store the results in the Cl's_xy array
        self.cls_gg = []
        self.cls_gs = []
        self.cls_ss = []
        for i1 in range(len(self.t_g)):
          for i2 in range(i1,len(self.t_g)):
            self.cls_gg.append(cls["dd"][(i1,i2)])
          for i2 in range(len(self.t_s)):
            self.cls_gs.append(cls["dl"][(i1,i2)])
        for i1 in range(len(self.t_s)):
          for i2 in range(i1,len(self.t_s)):
            self.cls_ss.append(cls["ll"][(i1,i2)])
        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)
        return
        
        
        

