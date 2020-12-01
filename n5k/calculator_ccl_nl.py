import numpy as np
import pyccl as ccl
from .calculator_ccl import N5KCalculatorCCL


class N5KCalculatorCCLNonLimber(N5KCalculatorCCL):
    name = 'CCLNonLimber'

    def run(self):
        # Compute power spectra
        ls = self.get_ells()
        # Radial interval (in Mpc)
        dchi = self.config.get('d_chi', 5.)
        # Radial interval (in Mpc)
        l_nonlimber = self.config.get('l_nonlimber', 100)

        self.cls_gg = []
        self.cls_gs = []
        self.cls_ss = []
        for i1, t1 in enumerate(self.t_g):
            for t2 in self.t_g[i1:]:
                self.cls_gg.append(ccl.angular_cl(self.cosmo, t1, t2, ls,
                                                  l_limber=l_nonlimber,
                                                  limber_integration_method='spline',
                                                  dchi_nonlimber=dchi))
            for t2 in self.t_s:
                self.cls_gs.append(ccl.angular_cl(self.cosmo, t1, t2, ls,
                                                  l_limber=l_nonlimber,
                                                  limber_integration_method='spline',
                                                  dchi_nonlimber=dchi))
        for i1, t1 in enumerate(self.t_s):
            for t2 in self.t_s[i1:]:
                self.cls_ss.append(ccl.angular_cl(self.cosmo, t1, t2, ls,
                                                  l_limber=l_nonlimber,
                                                  limber_integration_method='spline',
                                                  dchi_nonlimber=dchi))
        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)
