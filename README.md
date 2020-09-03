# DRAFT

# N5K: Non-local No-Nonsense Non-Limber Numerical Knockout

In this challenge, you are asked to compute a set of power spectra for a 3x2pt (cosmic shear, galaxy-galaxy lensing, and galaxy clustering) analysis without the use of the Limber approximation.

The input quantities follow those used in the Y10 LSST DESC Science Requirements Document v1 (link). There are 10 source redshift bins and 10 lens redshift bins. You should compute the data vector elements in the order given by the SRD (add further information), resulting in a data vector of 1200 (?) elements.

We provide dNdz's (from the SRD - link or include directly?) for each redshift bin of lens and sources as well as the cosmological parameters to be used.

Your code should use the Core Cosmology Library (CCL - link) to compute all other required inputs. You may either call CCL directly through its python interface or write its outputs to file to read in to your code. We provide a jupyter notebook which replicates the CCL calls used to obtain the quantities used as input in obtaining the benchmarks against which we compare. You may use these to save the outputs or to duplicate these commands. (Also provide these saved to file? Should we just provide a single fully calculated set of kernels (for each combination of dNdz's) to be integrated??)

## How to enter

Make a pull request to this repository which includes a file containing the Non-Limber-integrated data vector and a record of the timing for how long the task took (how to do this??). 

All entrants will have the opportunity to be an author on the resulting paper.

## Deadline

The challenge will close on December 15, 2020 (??)

## FAQ

** Can I participate if I am not a member of DESC? **

Yes, you can, and you can be an author on the paper (we have a special exemption from the pub policy for this challenge).




