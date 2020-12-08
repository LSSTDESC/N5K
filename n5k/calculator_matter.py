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
            print("FROM KERNELS")
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
            print("FROM TRACERS")
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
        # 1) Find out what is the minimum and maximum relevant value of chi for each window
        Ntg = len(self.t_g)
        Nts = len(self.t_s)
        Nchi_test = 20000
        threshold = 1e-10
        # growth rate = ccl.growth_rate(self.cosmo,a)
        age_test = ccl.comoving_radial_distance(self.cosmo,1e-4)*1.5 #Factor 1.5 for safety
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

        # 2) Define corresponding chi and scale factor sampling
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
            self.kerfac_g[i] += kern[itr]*trans[itr]*growth_g[itr]
        self.kerfac_s = np.zeros((Nts,Nchi_integrated))
        for i in range(Nts):
          kern = self.t_s[i].get_kernel(self.chi_integrated[i])
          trans = self.t_s[i].get_transfer(0.,ccl.scale_factor_of_chi(self.cosmo,self.chi_integrated[i]))
          for itr in range(1):#range(kern.shape[0]):
            self.kerfac_s[i] += kern[itr]*trans[itr]

        
        #import matplotlib.pyplot as plt
        #for i in range(Ntg):
        #  plt.plot(self.chi_nonintegrated[i],self.kerfac_g[i])
        ##maxidx = np.argmax(self.kerfac_g[0])
        ##maxval = self.kerfac_g[0][maxidx]
        ##overeidx = np.argmin(np.abs(maxval*np.exp(-1./2.)-self.kerfac_g[0]))
        ##chimax = self.chi_nonintegrated[0][maxidx]
        ##chiovere = self.chi_nonintegrated[0][overeidx]
        ##mu=chimax
        ##sig = np.abs(chimax-chiovere)
        ##plt.plot(self.chi_nonintegrated[0],self.kerfac_g[0])
        ##plt.plot(self.chi_nonintegrated[0],np.exp(-(self.chi_nonintegrated[0]-mu)**2/(2.*(sig)**2))*maxval)
        #plt.show()
        
        power = dpk['pk_lin'][::-1]
        Na_pk = len(power)

        # 4) Get the growth factor and pass it as well (sampled on same scale factor grid as the P(k))
        self.a_pk = a
        self.growth = ccl.growth_rate(self.cosmo, self.a_pk)
        pk_growth = np.empty((Na_pk,))
        for i in range(Na_pk):
          pk_growth[i] = np.sqrt(np.mean(power[i]/power[-1]))
        print(self.growth,pk_growth)
        self.growth = pk_growth
        growth_func = interp(self.a_pk,pk_growth)
        for i in range(Ntg):
          self.kerfac_g[i] *= growth_func(ccl.scale_factor_of_chi(self.cosmo,self.chi_nonintegrated[i]))
        for i in range(Nts):
          self.kerfac_s[i] *= growth_func(ccl.scale_factor_of_chi(self.cosmo,self.chi_integrated[i]))

        # 5) The kmin of the provided file is a bit too high. Here we extrapolate using k^(n_s) to reach lower k values
        kmin = 1e-7
        self.Nk_fft = 256
        Nk_small = int(np.log10(dpk['k'][0]/1e-7)/np.log10(dpk['k'][1]/dpk['k'][0])+1)
        assert(Nk_small > 10)
        ksmall = np.geomspace(1e-7,dpk['k'][0],endpoint=False,num=Nk_small)
        k_all = np.concatenate([ksmall,dpk['k']])
        self.tau_pk = ccl.comoving_radial_distance(self.cosmo,self.a_pk)
        self.k_pk = np.geomspace(kmin,dpk['k'][-1],num=self.Nk_fft)
        self.pk = np.empty((Na_pk,self.Nk_fft))
        for i in range(Na_pk):
          pk_all = np.concatenate([(ksmall/dpk['k'][0])**(par['n_s'])*power[i][0],power[i]])#np.concatenate([(ksmall/dpk['k'][0])**(par['n_s'])*dpk['pk_nl'][i][0],dpk['pk_nl'][i]])
          self.pk[i] = interp(k_all,pk_all)(self.k_pk)*self.k_pk**3/self.growth[i]**2
         
        #print(self.a_pk)
        ##import matplotlib.pyplot as plt
        ##for i in range(Na_pk):
        ##  plt.loglog(self.k_pk,self.pk[i]*self.k_pk**(-1.9)/self.growth[i]**2,color=plt.get_cmap("Greens")(float(i+1)/float(2*Na_pk)+0.5))
        ##plt.show()
        ##quit()
        #Special note : we want Delta(k)=P(k)*k^3 instead of P(k)
        

        # 6) Pass everything to the matterlib
        self.ma = matterlib.Matter(ma_verbose=3)

        self.ma.set((self.chi_nonintegrated,self.chi_integrated),(self.kerfac_g,self.kerfac_s),(self.a_pk,self.tau_pk,self.k_pk,self.pk),self.growth,lmax=2000,size_fft_cutoff=50,tw_size=25,integrated_tw_size=75)

        ##arr = np.loadtxt("window.dat").T
        ##import matplotlib.pyplot as plt
        ##plt.plot(14000.-arr[0],arr[1])
        ##plt.plot(self.chi_nonintegrated[0],self.kerfac_g[0])
        ##plt.show()
    def run(self):
        # Compute power spectra
        ls = self.get_ells()

        self.ma.compute()
        cls = self.ma.matter_cl(ls)
        
        fft_r,fft_i = self.ma.get_fft()
        fft = fft_r + 1j* fft_i
        N = self.Nk_fft
        tilt = 1.9
        kmin,kmax = 1e-7,1e2
        dlogk = np.log(kmax/kmin)
        Nfft = len(fft)
        nu_n_tilde_tilde = 2j*np.pi*np.arange(N)/dlogk * (N-1)/N 

        cn_tilde = np.array([np.sum(self.pk[-1]*self.k_pk**(-1.9)*np.exp(-2j*np.pi*np.arange(N)*n/N)) for n in range(N)])
        cn = cn_tilde/N*(kmin)**(-nu_n_tilde_tilde)

        #fnij = integrate()
        nzs = self.get_tracer_dndzs()
        tpar = self.get_tracer_parameters()
        z_g = nzs['z_cl']
        tr = ccl.NumberCountsTracer(self.cosmo, True, (z_g, nzs['dNdz_cl'][0, :]),
                                               bias=(z_g, np.full(len(z_g), tpar['b_g'][0])))
        
        comp_cl = ccl.angular_cl(self.cosmo, tr, tr, ls)
        print(comp_cl,cls['dd'][(0,0)])
        print(cls['dd'][(0,0)][-1]/comp_cl[-1])

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
        
        
        

