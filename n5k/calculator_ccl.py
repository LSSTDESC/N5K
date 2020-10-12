import numpy as np
import pyccl as ccl
from .calculator_base import N5KCalculatorBase


class N5KCalculatorCCL(N5KCalculatorBase):
    name = 'CCL'

    def setup(self):
        # Initialize cosmology
        self.cosmo = ccl.Cosmology(Omega_c=self.cosmo['OmegaM']-self.cosmo['OmegaB'],
                                   Omega_b=self.cosmo['OmegaB'],
                                   h=self.cosmo['h'], n_s=self.cosmo['n_s'],
                                   A_s=self.cosmo['A_s'], w0=self.cosmo['w0'])
        s8 = ccl.sigma8(self.cosmo)

        # Initialize tracers
        z_g, nz_g = self.get_nz_g()
        z_s, nz_s = self.get_nz_s()
        b_g = self.get_bias()
        A_IA = self.get_A_IA()
        self.t_g = [ccl.NumberCountsTracer(self.cosmo, True, (z_g, n),
                                           bias=(z_g, np.full(len(z_g), b)))
                    for n, b in zip(nz_g, b_g)]
        self.t_s = [ccl.WeakLensingTracer(self.cosmo, (z_g, n), True,
                                          (z_s, np.full(len(z_s), A_IA)))
                    for n in nz_s]

    def run(self):
        # Compute power spectra
        ls = self.get_ells()

        self.cls_gg = []
        self.cls_gs = []
        self.cls_ss = []
        for i1, t1 in enumerate(self.t_g):
            for t2 in self.t_g[i1:]:
                self.cls_gg.append(ccl.angular_cl(self.cosmo, t1, t2, ls))
            for t2 in self.t_s:
                self.cls_gs.append(ccl.angular_cl(self.cosmo, t1, t2, ls))
        for i1, t1 in enumerate(self.t_s):
            for t2 in self.t_s[i1:]:
                self.cls_ss.append(ccl.angular_cl(self.cosmo, t1, t2, ls))
        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)
