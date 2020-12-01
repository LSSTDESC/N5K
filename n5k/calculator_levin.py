import numpy as np
import pyccl as ccl
from .calculator_base import N5KCalculatorBase

import levinpower

class N5KCalculatorLevin(N5KCalculatorBase):
    name = 'Levin'

    def setup(self):
        # Initialize cosmology
        par = self.get_cosmological_parameters()

        pk = self.get_pk()
        kernels = self.get_tracer_kernels()
        background = self.get_background()

        number_count = kernels["kernels_cl"].shape[0]
        self.levin_calculator = levinpower.LevinPower(number_count, 
                          background["z"], background["chi"], 
                          kernels["chi_cl"], np.concatenate((kernels["kernels_cl"].T, kernels["kernels_sh"].T), axis=1), 
                          pk["k"], pk["z"], pk["pk_lin"].flatten(), pk["pk_nl"].flatten())


    def run(self):
        # Compute power spectra
        ell = self.get_ells().astype(int)

        self.cls_gg, self.cls_gs, self.cls_ss = self.levin_calculator.compute_C_ells(ell)

        self.cls_gg = np.array(self.cls_gg)
        self.cls_gs = np.array(self.cls_gs)
        self.cls_ss = np.array(self.cls_ss)
