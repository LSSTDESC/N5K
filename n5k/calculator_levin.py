import numpy as np
from .calculator_base import N5KCalculatorBase

import levinpower

class N5KCalculatorLevin(N5KCalculatorBase):
    name = 'Levin'

    def setup(self):
        # Initialize cosmology
        pk = self.get_pk()
        kernels = self.get_tracer_kernels()
        background = self.get_background()

        number_count = kernels["kernels_cl"].shape[0]

        precompute_splines = self.config.get('precompute_splines', False)
        ell = self.get_ells().astype(int)
        self.levin_calculator = levinpower.LevinPower(
                          precompute_splines,
                          ell,
                          number_count, 
                          background["z"], background["chi"],
                          kernels["chi_cl"],
                          np.concatenate((kernels["kernels_cl"].T,
                                          kernels["kernels_sh"].T), axis=1),
                          pk["k"], pk["z"],
                          pk["pk_lin"].flatten(),
                          pk["pk_nl"].flatten())

    def run(self):
        # Compute power spectra
        parallelize_ell = self.config.get('parallelize_ell', False)
        self.cls_gg, self.cls_gs, self.cls_ss = \
            self.levin_calculator.compute_C_ells(parallelize_ell)

        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)
