"""
.. module:: matterlib
    :synopsis: Python wrapper around Matter module of CLASS
.. moduleauthor:: Nils Schöneberg <schoeneberg@physik.rwth-aachen.de>
"""
from math import exp,log
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i

# Import the .pxd containing definitions
from mattery cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass


cdef class Matter:
    """
    from matter import Matter
    """
    cdef matters ma
    cpdef int computable

    # Called at the end of a run, to free memory
    def clean(self):
        free(self.ma.sampled_sources)
        free(self.ma.k_sampling)
        free(self.ma.logk_sampling)
        free(self.ma.tw_sampling)
        free(self.ma.integrated_tw_sampling)
        free(self.ma.exp_integrated_tw_sampling)
        free(self.ma.tw_weights)
        free(self.ma.integrated_tw_weights)
        free(self.ma.sampled_sources)
        #matter_free(&self.ma)

    def __init__(self,ma_verbose=0):
        """
        """
        self.computable = False
        self.ma.uses_intxi_logarithmic = True
        self.ma.matter_verbose = ma_verbose

    def set(self,chi,kfac,pk,growth,lmax=50):
        """
        """
        chi_ni,chi_i = chi
        ntr_ni = len(chi_ni)
        ntr_i = len(chi_i)
        nchi_ni = len(chi_ni[0])
        nchi_i = len(chi_i[0])
        kfac_g,kfac_s = kfac
        a_pk,tau_pk, k_pk, pk = pk

        self.ma.size_fft_input = len(k_pk)
        self.ma.size_fft_cutoff = 100
        assert(self.ma.size_fft_cutoff < self.ma.size_fft_input)
        self.ma.tw_size = nchi_ni
        self.ma.integrated_tw_size = nchi_i
        self.ma.tau0 = 14000.

        print("Assinging tw sampling")

        self.ma.tw_size = 25
        self.ma.integrated_tw_size = 75
        ntw_ni = self.ma.tw_size
        ntw_i = self.ma.integrated_tw_size
        self.ma.tw_min = <double*>malloc(sizeof(double)*ntr_ni)
        self.ma.tw_max = <double*>malloc(sizeof(double)*ntr_ni)
        self.ma.num_windows = max(ntr_ni,ntr_i)
        self.ma.tw_sampling = <double*>malloc(sizeof(double)*self.ma.num_windows*ntw_ni)
        self.ma.integrated_tw_sampling = <double*>malloc(sizeof(double)*self.ma.num_windows*ntw_i)
        self.ma.exp_integrated_tw_sampling = <double*>malloc(sizeof(double)*self.ma.num_windows*ntw_i)
        self.ma.tw_weights = <double*>malloc(sizeof(double)*self.ma.num_windows*ntw_ni)
        self.ma.integrated_tw_weights = <double*>malloc(sizeof(double)*self.ma.num_windows*ntw_i)
        self.ma.ptw_sampling = <double*>malloc(sizeof(double)*self.ma.num_windows*nchi_ni)
        self.ma.ptw_integrated_sampling = <double*>malloc(sizeof(double)*self.ma.num_windows*nchi_i)
        tw_ni = [self.ma.tau0-np.linspace(chi_ni[nwd][0],chi_ni[nwd][-1],num=ntw_ni)[::-1] for nwd in range(ntr_ni)]
        tw_i = [self.ma.tau0-np.linspace(chi_i[nwd][0],chi_i[nwd][-1],num=ntw_i)[::-1] for nwd in range(ntr_i)]
        tw_ni_weights = [self.weights(tw_ni[i]) for i in range(ntr_ni)]
        tw_i_weights = [self.weights(tw_i[i]) for i in range(ntr_i)]
        # Can cut off the tw_i already earlier due to tilt tw^(1-nu) factor

        for nwd in range(ntr_ni):
          for itw in range(ntw_ni):
            self.ma.tw_sampling[nwd*ntw_ni+itw] = tw_ni[nwd][itw]
            self.ma.tw_weights[nwd*ntw_ni+itw] = tw_ni_weights[nwd][itw]
          self.ma.tw_max[nwd]=tw_ni[nwd][-1]
          self.ma.tw_min[nwd]=tw_ni[nwd][0]
        for nwd in range(ntr_i):
          for itw in range(ntw_i):
            if self.ma.uses_intxi_logarithmic:
                self.ma.exp_integrated_tw_sampling[nwd*ntw_i+itw] = self.ma.tau0-tw_i[nwd][itw]
                self.ma.integrated_tw_sampling[nwd*ntw_i+itw] = np.log(self.ma.tau0-tw_i[nwd][itw]+1.0e-5) ## TODO :: remove this factor
                self.ma.integrated_tw_weights[nwd*ntw_i+itw] = tw_i_weights[nwd][itw]
            else:
                self.ma.integrated_tw_sampling[nwd*ntw_i+itw] = tw_i[nwd][itw]
                self.ma.integrated_tw_weights[nwd*ntw_i+itw] = tw_i_weights[nwd][itw]
        for nwd in range(ntr_i,ntr_ni): # TODO :: remove
          for itw in range(ntw_i):
            if self.ma.uses_intxi_logarithmic:
                self.ma.exp_integrated_tw_sampling[nwd*ntw_i+itw] = 0.
            self.ma.integrated_tw_sampling[nwd*ntw_i+itw] = 0.
            self.ma.integrated_tw_weights[nwd*ntw_i+itw] = 0.
        for nwd in range(ntr_ni):
          for ichi in range(nchi_ni):
            self.ma.ptw_sampling[nwd*nchi_ni+ichi] = self.ma.tau0-chi_ni[nwd][::-1][ichi]
        for nwd in range(ntr_i):
          for ichi in range(nchi_i):
            self.ma.ptw_integrated_sampling[nwd*nchi_i+ichi] = self.ma.tau0-chi_i[nwd][::-1][ichi]
        for nwd in range(ntr_i,ntr_ni):
          for ichi in range(nchi_i):
            self.ma.ptw_integrated_sampling[nwd*nchi_i+ichi] = 0.

        self.ma.tau_size = len(a_pk)
        self.ma.tau_sampling = <double*>malloc(sizeof(double)*len(a_pk))
        self.ma.sampled_sources = <double*>malloc(sizeof(double)*len(a_pk)*len(k_pk))
        for ia in range(len(a_pk)):
          self.ma.tau_sampling[ia] = self.ma.tau0-tau_pk[ia]
          for ik in range(len(k_pk)):
            self.ma.sampled_sources[ik*len(a_pk)+ia] = pk[ia,ik]
        self.ma.k_sampling = <double*>malloc(sizeof(double)*len(k_pk))
        self.ma.logk_sampling = <double*>malloc(sizeof(double)*len(k_pk))
        for ik in range(len(k_pk)):
          self.ma.k_sampling[ik] = k_pk[ik]
          self.ma.logk_sampling[ik] = np.log(k_pk[ik])
        self.ma.deltalogk = self.ma.logk_sampling[len(k_pk)-1]-self.ma.logk_sampling[0]
        self.ma.uses_separability = 1
        self.ma.non_diag = self.ma.num_windows -1
        self.ma.has_cls = True
        self.ma.l_lss_max = lmax
        self.ma.bias = 1.9
        self.ma.l_logstep = 1.12
        self.ma.l_linstep = 40
        self.ma.has_unintegrated_windows = 1
        self.ma.has_integrated_windows = 1
        self.ma.uses_limber_approximation = 0
        self.ma.t_size = 250
        self.ma.radtp_size_total = 2
        self.ma.ptw_window = <double**>malloc(sizeof(double*)*self.ma.radtp_size_total)
        for i in range(self.ma.radtp_size_total):
          size = (nchi_i if self.matter_is_integrated(i) else nchi_ni)
          self.ma.ptw_window[i] = <double*>malloc(sizeof(double)*size*self.ma.num_windows)
          if self.matter_is_integrated(i):
            for nwd in range(ntr_i):
              for ichi in range(nchi_i):
                self.ma.ptw_window[i][nwd*nchi_i+ichi] = kfac_s[nwd][::-1][ichi]
            for nwd in range(ntr_i,ntr_ni):
              for ichi in range(nchi_i):
                self.ma.ptw_window[i][nwd*nchi_i+ichi] = 0.
          else:
            for nwd in range(ntr_ni):
              for ichi in range(nchi_ni):
                self.ma.ptw_window[i][nwd*nchi_ni+ichi] = kfac_g[nwd][::-1][ichi]
        self.ma.growth_factor_tau = <double*>malloc(sizeof(double)*self.ma.tau_size)
        for ia in range(len(a_pk)):
          self.ma.growth_factor_tau[ia] = growth[ia]
        
        self.computable = True

    def weights(self,x):

        warray = np.full((len(x),),(x[-1]-x[0])/(len(x)-1))
        warray[0]*=0.5
        warray[len(x)-1]*=0.5

        return warray
        
    def matter_is_integrated(self,radtp):
      if radtp == 0:
        return True
      else:
        return False

    def compute(self):
        """
        """
        if not self.computable:
          return False
        cdef ErrorMsg errmsg

        if matter_init(&(self.ma)) == _FAILURE_:
            raise CosmoComputationError(self.ma.error_message)

        return

    def matter_cl(self, ells, nofail=False):
        def index_symmetric_matrix(a,b,N):
            if(a <= b):
                return b+N*a-(a*(a+1))//2
            else:
                return a+N*b-(b*(b+1))//2
        """
        matter_density_cl(ells, nofail=False)

        Return a dictionary of the primary number count/shear C_l for the matter structure

        Parameters
        ----------
        ells : array
            Define the l array for which the C_l will be returned
        nofail: bool, optional
            Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : numpy array of numpy.ndarrays
            Array that contains the list (in this order) of self correlation of
            1st bin, then successive correlations (set by non_diagonal) to the
            following bins, then self correlation of 2nd bin, etc. The array
            starts at index_ct_dd.
        """
        cdef int nl = len(ells)
        cdef double **dcl = <double**> calloc(self.ma.cltp_grid_size,sizeof(double*))
        for index_cltp_grid in range(self.ma.cltp_grid_size):
            dcl[index_cltp_grid] = <double*> calloc(self.ma.window_size[index_cltp_grid]*nl, sizeof(double))

        lmaxR = self.ma.l_lss_max

        if (not self.ma.has_cltp_nc) and (not self.ma.has_cltp_sh):
            raise CosmoSevereError("No density Cl computed with matters struct")
        if ells[-1] > lmaxR:
            if nofail:
                self._pars_check("l_max_lss",ells[-1])
                self._pars_check("output",'nCl')
                self.compute()
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d, but l requested was %d"%(lmaxR,ells[-1]))

        cl = {}
        spectra = []
        if self.ma.has_cltp_nc:
            spectra.append('dd')
            if self.ma.has_cltp_sh:
                spectra.append('dl')
        if self.ma.has_cltp_sh:
            spectra.append('ll')

        indices_cltp_grid = {}
        if self.ma.has_cltp_nc:
            index_cltp_grid = index_symmetric_matrix(self.ma.cltp_index_nc,self.ma.cltp_index_nc,self.ma.cltp_size)
            indices_cltp_grid['dd']=index_cltp_grid
            if self.ma.has_cltp_sh:
                index_cltp_grid = index_symmetric_matrix(self.ma.cltp_index_nc,self.ma.cltp_index_sh,self.ma.cltp_size)
                indices_cltp_grid['dl']=index_cltp_grid
        if self.ma.has_cltp_sh:
            index_cltp_grid = index_symmetric_matrix(self.ma.cltp_index_sh,self.ma.cltp_index_sh,self.ma.cltp_size)
            indices_cltp_grid['ll']=index_cltp_grid


        ell_view = <double*>malloc(nl*sizeof(double))
        for il in range(nl):
          ell_view[il] = ells[il]
        if matter_cl_at_l(&self.ma, ell_view, nl, dcl) == _FAILURE_:
            raise CosmoSevereError(self.ma.error_message)
        for elem in spectra:
            cl[elem] = {}
            wd_size_elem = self.ma.window_size[indices_cltp_grid[elem]]
            index = 0
            for index_wd1 in range(self.ma.num_windows):
              for index_wd2 in range(self.ma.window_index_start[indices_cltp_grid[elem]][index_wd1],self.ma.window_index_end[indices_cltp_grid[elem]][index_wd1]+1):
                cl[elem][(index_wd1,index_wd2)] = np.zeros(nl,dtype="float64")
                for il,ellval in enumerate(ells):
                  cl[elem][(index_wd1,index_wd2)][il] = dcl[indices_cltp_grid[elem]][index+il*wd_size_elem]
                index+=1

        cl['ell'] = ells

        for index_cltp_grid in range(self.ma.cltp_grid_size):
            free(dcl[index_cltp_grid])
        free(dcl)

        return cl
