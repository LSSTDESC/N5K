import numpy as np
import pyccl as ccl
from .calculator_base import N5KCalculatorBase
from scipy.interpolate import interp1d
from fftlog_fang import fftlog_fang
import time

class N5KCalculatorFKEM(N5KCalculatorBase):
		name = 'FKEM' # Fang, Krause, Eifler, MacCrann

		def setup(self):
				# Initialize cosmology
				par = self.get_cosmological_parameters()
				dpk = self.get_pk()
				a = 1./(1+dpk['z'][::-1])
				self.cosmo = ccl.CosmologyCalculator(Omega_c=par['Omega_m']-par['Omega_b'],
																						 Omega_b=par['Omega_b'],
																						 h=par['h'], n_s=par['n_s'],
																						 A_s=par['A_s'], w0=par['w0'],
																						 pk_linear={'a': a,
																												'k': dpk['k'],
																												'delta_matter:delta_matter': dpk['pk_lin'][::-1][:]},
																						 pk_nonlin={'a': a,
																												'k': dpk['k'],
																												'delta_matter:delta_matter': dpk['pk_nl'][::-1][:]})

				self.cosmo_nonlin_part = ccl.CosmologyCalculator(Omega_c=par['Omega_m']-par['Omega_b'],
																						 Omega_b=par['Omega_b'],
																						 h=par['h'], n_s=par['n_s'],
																						 A_s=par['A_s'], w0=par['w0'],
																						 pk_linear={'a': a,
																												'k': dpk['k'],
																												'delta_matter:delta_matter': dpk['pk_lin'][::-1][:]},
																						 pk_nonlin={'a': a,
																												'k': dpk['k'],
																												'delta_matter:delta_matter': (dpk['pk_lin'][::-1][:])})

																						 # pk_nonlin={'a': a,
																							# 					'k': dpk['k'],
																							# 					'delta_matter:delta_matter': (dpk['pk_nl'][::-1][:]-dpk['pk_lin'][::-1][:])})


				chi_min = ccl.comoving_radial_distance(cosmo=self.cosmo, a=1./(1.+0.002))
				chi_max = ccl.comoving_radial_distance(cosmo=self.cosmo, a=1./(1.+4.))
				Nchi = 800
				self.chi_logspace_arr = np.logspace(np.log10(chi_min),np.log10(chi_max),num=Nchi, endpoint=True)
				self.dlnr = np.log(chi_max/chi_min)/(Nchi-1.)
				a_arr = ccl.scale_factor_of_chi(self.cosmo, self.chi_logspace_arr)
				growfac_arr = ccl.growth_factor(self.cosmo, a_arr)

				self.plin0 = dpk['pk_lin'][0]
				self.k_in = dpk['k']
				self.pklin_interp = interp1d(np.log(self.k_in), np.log(self.plin0), fill_value='extrapolate')

				# Initialize tracers
				if self.config.get('tracers_from_kernels', False):
						tpar = self.get_tracer_parameters()
						ker = self.get_tracer_kernels()
						a_g = 1./(1+ker['z_cl'][::-1])
						self.t_g = []
						self.t_g_logspace = []


						for k in ker['kernels_cl']:
								t = ccl.Tracer()
								barr = np.ones_like(a_g)

								t.add_tracer(self.cosmo,
														 (ker['chi_cl'], k),
														 transfer_a=(a_g, barr))
								self.t_g.append(t)

								ker_at_chi = interp1d(np.log(ker['chi_cl']), k, bounds_error=False, fill_value=0.)
								ker_logspace = ker_at_chi(np.log(self.chi_logspace_arr)) * self.chi_logspace_arr * growfac_arr
								self.t_g_logspace.append(ker_logspace)

						self.t_s = []
						self.t_s_logspace = []
						for k in ker['kernels_sh']:
								t = ccl.Tracer()
								t.add_tracer(self.cosmo,
														 kernel=(ker['chi_sh'], k),
														 der_bessel=-1, der_angles=2)
								self.t_s.append(t)

								ker_at_chi = interp1d(np.log(ker['chi_sh']), k,bounds_error=False, fill_value=0.)
								ker_logspace = ker_at_chi(np.log(self.chi_logspace_arr)) * self.chi_logspace_arr * growfac_arr
								self.t_s_logspace.append(ker_logspace)

				else:
						print("calculator_fang.py: Nonlimber not supported!")
						sys.exit()
						nzs = self.get_tracer_dndzs()
						tpar = self.get_tracer_parameters()
						z_g = nzs['z_cl']
						z_s = nzs['z_sh']
						self.t_g = [ccl.NumberCountsTracer(self.cosmo, True,
																							 (z_g, nzs['dNdz_cl'][:, ni]),
																							 bias=(z_g,
																										 np.full(len(z_g), b)))
												for ni, b in zip(range(0, 10),
																				 tpar['b_g'])]
						self.t_s = [ccl.WeakLensingTracer(self.cosmo,
																							(z_s, nzs['dNdz_sh'][:, ni]),
																							True)
												for ni in range(0, 5)]


		def _get_cl_limber(self, t1, t2, ls):
				return ccl.angular_cl(self.cosmo, t1, t2, ls,
															limber_integration_method='spline')
		def _get_cl_limber_nonlin_part(self, t1, t2, ls):
				return ccl.angular_cl(self.cosmo, t1, t2, ls, limber_integration_method='spline') - ccl.angular_cl(self.cosmo_nonlin_part, t1, t2, ls,limber_integration_method='spline')

		def _get_cls_nonlimber_lin(self, fchi1, fchi2, ls):
				nu = 1.01
				myfftlog1 = fftlog_fang(self.chi_logspace_arr, fchi1, nu=nu, N_extrap_low=0, N_extrap_high=0, c_window_width=0.25, N_pad=112)
				myfftlog2 = fftlog_fang(self.chi_logspace_arr, fchi2, nu=nu, N_extrap_low=0, N_extrap_high=0, c_window_width=0.25, N_pad=112)
				Nell = ls.size
				cls_nonlimber_lin = np.zeros(Nell)
				# t0 = time.time()
				for i in np.arange(Nell):
						ell = ls[i]
						k, fk1 = myfftlog1.fftlog(ell)
						k, fk2 = myfftlog2.fftlog(ell)
						cls_nonlimber_lin[i] = np.sum(fk1 * fk2 * k**3 * np.exp(self.pklin_interp(np.log(k))) ) * self.dlnr * 2./np.pi
				# t1 = time.time()
				# print(Nell, t1-t0)
				return cls_nonlimber_lin


		def _get_cl_nonlimber_fang(self, t1, t2, ls, i1, i2, cl_type):
				ls_nonlim = ls[ls<200]
				cls_limber_nonlin_part = self._get_cl_limber_nonlin_part(t1, t2, ls_nonlim)
				if(cl_type=='gg'):
					cls_nonlimber_lin = self._get_cls_nonlimber_lin(self.t_g_logspace[i1], self.t_g_logspace[i1+i2], ls_nonlim)
					cls_nonlimber = cls_limber_nonlin_part + cls_nonlimber_lin
					cls_limber = self._get_cl_limber(t1, t2, ls[ls>=200])
					return np.hstack((cls_nonlimber, cls_limber))
				# elif(cl_type=='gs'):
				# 	cls_nonlimber_lin =  self._get_cls_nonlimber_lin(self.t_g_logspace[i1], self.t_s_logspace[i2], ls_nonlim)
				# else:
				# 	cls_nonlimber_lin =  self._get_cls_nonlimber_lin(self.t_s_logspace[i1], self.t_s_logspace[i1+i2], ls_nonlim)
				else:
					return self._get_cl_limber(t1, t2, ls)
				# cls_nonlimber = cls_limber_nonlin_part + cls_nonlimber_lin
				# print(cls_nonlimber.size)
				# print(ls_nonlim.size)
				# exit()
				# cls_limber = self._get_cl_limber(t1, t2, ls[ls>=200])
				# return np.hstack((cls_nonlimber, cls_limber))

		def run(self):
				# Compute power spectra
				ls = self.get_ells()

				self.cls_gg = []
				self.cls_gs = []
				self.cls_ss = []
				for i1, t1 in enumerate(self.t_g):
						for i2, t2 in enumerate(self.t_g[i1:]):
								self.cls_gg.append(self._get_cl_nonlimber_fang(t1, t2, ls, i1, i2, 'gg'))
						for i2, t2 in enumerate(self.t_s):
								self.cls_gs.append(self._get_cl_nonlimber_fang(t1, t2, ls, i1, i2, 'gs'))
				for i1, t1 in enumerate(self.t_s):
						for i2, t2 in enumerate(self.t_s[i1:]):
								self.cls_ss.append(self._get_cl_nonlimber_fang(t1, t2, ls, i1, i2, 'ss'))
				self.cls_gg = np.array(self.cls_gg)
				self.cls_gs = np.array(self.cls_gs)
				self.cls_ss = np.array(self.cls_ss)
