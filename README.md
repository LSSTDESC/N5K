# DRAFT

# N5K: Non-local No-Nonsense Non-Limber Numerical Knockout

In this challenge, you are asked to provide a code to compute a set of power spectra for a 3x2pt (cosmic shear, galaxy-galaxy lensing, and galaxy clustering) analysis without the use of the Limber approximation.

The challenge entries will be evaluated on the basis of accuracy and efficiency. The code which can accomplish this task fastest and within the accuracy requirements of an LSST Y10 cosmological analysis will win the challenge. The winning code will be incorporated for use as the non-Limber integration tool for the [Core Cosmology Library](https://github.com/LSSTDESC/CCL/).

The accuracy and speed of the entry codes will be determined on the basis of a 3x2pt analysis setup following that used for the LSST Y10 scenario in the [LSST DESC Science Requirements Document v1](https://arxiv.org/pdf/1809.01669.pdf). We include in this repo several of the required data inputs from this calculation, derived from those included with the Science Requirements Document [Data Products Release](https://zenodo.org/record/2662127#.X2NtDobTWEA) (LICENSE!).

There are 10 source redshift bins and 10 lens redshift bins. Their redshift distributions are found in the files X.

The $\ell$ values at which you should compute the spectra are given in the file 'input/ell-values'

You should output the data vector elements in the same order as that used for the 3x2pt elements in the Science Requirements Document Data Product Release, described here in section C of the file 'input/README.md' (sourced from the SRD Data Products Release) and accompanying 'input/gglensing_zbin_y10' (IS THIS THE UP TO DATE VERSION?). The resulting data vector will have 1000 elements.

Your code should use the Core Cosmology Library (CCL - link) to compute all other required inputs. You may either call CCL directly through its python interface or write its outputs to file to read in to your code. We provide a jupyter notebook which replicates the CCL calls used to obtain the quantities used as input in obtaining the benchmarks against which we compare. You may use these to save the outputs or to duplicate these commands. (Also provide these saved to file? Should we just provide a single fully calculated set of kernels (for each combination of dNdz's) to be integrated??)

## How to enter

Make a pull request to this repository which includes (WHAT). 

All entrants will have the opportunity to be an author on the resulting paper.

## Deadline

The challenge will close on December 15, 2020 (CONFIRM).

## FAQ

** Can I participate if I am not a member of DESC? **

Yes, you can, and you can be an author on the paper (we have a special exemption from the pub policy for this challenge).

