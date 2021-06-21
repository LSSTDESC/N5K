import numpy as np
import os.path
import pickle

from .calculator_base import N5KCalculatorBase
from .calculator_ccl import N5KCalculatorCCL
from scipy.special import spherical_jn
from .bessel_tools import bessel_zeros
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import simps

# maximum l at which we are going to use non-limber
lmax = 200

class N5KCalculatorDSBT(N5KCalculatorBase):
    name = 'DSBT'

    def setup(self):

        # We use the default CCL limber integrator
        self.limber = N5KCalculatorCCL("tests/config_ccl_limber.yml")
        self.limber.setup()

        # If transform matrix is already trained, we just load it
        if os.path.exists('cache_dsbt.npz'):
            cache = np.load('cache_dsbt.npz')
            self.T = cache['T']
            self.rln = cache['rln']
        else:
            # Retrieve relevant info
            pk = self.get_pk()
            ells = self.get_ells()

            # Here we are 
            ells = ells[ells<lmax]

            kernels = self.get_tracer_kernels()
            chi_kernel = kernels['chi_sh']

            k = pk['k']
            kmax = k[-1] # What's our maximum scale?

            # Ideally for the transform to work well this should be 1
            # but that's costly, so we can apply a factor that neglects small
            # scales... we won't get very accurate results at high ell
            res_factor = 50
            kmax = kmax / res_factor
            # Nmax at l=0, this is conservative
            nmax = int(chi_kernel[-1] * kmax / np.pi)

            # Retrieve once and for all the qln we need
            qln = bessel_zeros(nmax, ells)
            rln = qln / kmax

            # Precomputing transformation matrix
            pref_factor = np.sqrt(2*np.pi) /( kmax**3 * np.stack([spherical_jn(int(ells[i]+1), qln[i])
              for i in range(len(ells)) ])**2)
            # Compute Bessel term
            bessel = np.stack([spherical_jn(int(ells[i]), np.outer(rln[i], k))
                      for i in range(len(ells))])

            # And here we have our transform matrix
            T = np.expand_dims(pref_factor, -1) * bessel

            print("Precomputed transformation matrix of size", T.shape)

            self.T = T
            self.rln = rln
            # Export matrix to disk
            np.savez('cache_dsbt.npz', rln=rln, T=T)

    def run(self):
        pk = self.get_pk()
        kernels = self.get_tracer_kernels()
        background = self.get_background()
        ells = self.get_ells()

        k = pk['k']

        z2chi = interp1d(background['z'], background['chi'])
        pk_interp = interp2d(pk['k'], z2chi(pk['z']), pk['pk_nl'])
        kernels_interp_cl = [interp1d(kernels['chi_cl'], kernels['kernels_cl'][i]/kernels['chi_cl']**2,
                      fill_value=0., bounds_error=False ) for i in range(len(kernels['kernels_cl']))]
        kernels_interp_sh = [interp1d(kernels['chi_sh'], kernels['kernels_sh'][i]/kernels['chi_sh']**2,
                      fill_value=0., bounds_error=False ) for i in range(len(kernels['kernels_sh']))]

        # We first compute all the cls with non limber up to lmax
        ells = ells[ells<lmax]
        self.cls_gg = []
        self.cls_gs = []
        self.cls_ss = []
        for i1, t1 in enumerate(kernels_interp_cl):
            int1 = (self.T * np.stack([np.expand_dims(t1(self.rln[i]),-1) * np.sqrt(pk_interp(k, self.rln[i])) for i in
                                      range(len(ells))])).sum(axis=1)

            for t2 in kernels_interp_cl[i1:]:
                int2 = (self.T * np.stack([np.expand_dims(t2(self.rln[i]),-1) * np.sqrt(pk_interp(k, self.rln[i])) for i in
                                      range(len(ells))])).sum(axis=1)
                self.cls_gg.append(simps(int1*int2*k**2, k))

            for t2 in kernels_interp_sh:
                int2 = (self.T * np.stack([np.expand_dims(t2(self.rln[i]),-1) * np.sqrt(pk_interp(k, self.rln[i])) for i in
                                                      range(len(ells))])).sum(axis=1)
                self.cls_gs.append(simps(int1*int2*k**2, k))

        for i1, t1 in enumerate(kernels_interp_sh):
            int1 = (self.T * np.stack([np.expand_dims(t1(self.rln[i]),-1) * np.sqrt(pk_interp(k, self.rln[i])) for i in
                                      range(len(ells))])).sum(axis=1)
            for t2 in kernels_interp_sh[i1:]:
                int2 = (self.T * np.stack([np.expand_dims(t2(self.rln[i]),-1) * np.sqrt(pk_interp(k, self.rln[i])) for i in
                                      range(len(ells))])).sum(axis=1)
                self.cls_ss.append(simps(int1*int2*k**2, k))

        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)

        ells = self.get_ells()
        ells = ells[ells>=lmax]

        cls_gg = []
        cls_gs = []
        cls_ss = []
        for i1, t1 in enumerate(self.limber.t_g):
            for t2 in self.limber.t_g[i1:]:
                cls_gg.append(self.limber._get_cl(t1, t2, ells))
            for t2 in self.limber.t_s:
                cls_gs.append(self.limber._get_cl(t1, t2, ells))
        for i1, t1 in enumerate(self.limber.t_s):
            for t2 in self.limber.t_s[i1:]:
                cls_ss.append(self.limber._get_cl(t1, t2, ells))
        cls_gg = np.array(cls_gg)
        cls_gs = np.array(cls_gs)
        cls_ss = np.array(cls_ss)

        self.cls_gg = np.concatenate([self.cls_gg, cls_gg], axis=-1)
        self.cls_gs = np.concatenate([self.cls_gs, cls_gs], axis=-1)
        self.cls_ss = np.concatenate([self.cls_ss, cls_ss], axis=-1)