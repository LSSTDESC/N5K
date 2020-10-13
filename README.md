# DRAFT

# N5K: Non-local No-Nonsense Non-Limber Numerical Knockout

In this challenge, you are asked to compute a set of power spectra for a 3x2pt (cosmic shear, galaxy-galaxy lensing, and galaxy clustering) analysis without the use of the Limber approximation.

The challenge entries will be evaluated on the basis of accuracy and speed. The code which can accomplish this task fastest and within the accuracy requirements of an LSST Y10 cosmological analysis will win the challenge. The winning code will be incorporated for use as the non-Limber integration tool for the [Core Cosmology Library](https://github.com/LSSTDESC/CCL/).

All entrants will have the opportunity to be an author on the resulting paper.

## How to enter

The challenge asks you to compute the $C-\ell$s required for a 3x2pt analysis setup similar to the LSST Y10 scenario in the [LSST DESC Science Requirements Document v1](https://arxiv.org/pdf/1809.01669.pdf). The 'input' folder of this repo contains some required inputs for this calculation, derived in some cases from those included with the Science Requirements Document [Data Products Release](https://zenodo.org/record/2662127#.X2NtDobTWEA):
- kernels for the 10 number counts tracers and 5 weak lensing tracers as a function of z  (N5K/input/kernels.npz)
- linear and nonlinear version of the matter power spectrum as a function of k and z (N5K/input/Pk.npz)

N5K/n5k/calculator\_base.py contains a base class N5KCalculatorBase. Write a subclass of N5KCalculatorBase which contains methods setup() (to set up your nonlimber calculation) and run() (to run it). N5K/n5k/calculator\_ccl.py contains an example of what this would look like doing the calculation using CCL's current (Limber) $C_\ell$ calculation tools.

Make a pull request to this repository which includes your new subclass. 

## Deadline

The challenge will close on December 15, 2020 (CONFIRM).

## FAQ

** Can I participate if I am not a member of DESC? **

Yes, you can, and you can be an author on the paper (we have a special exemption from the pub policy for this challenge).

## License information

We make use of both modified and original data products from the LSST Dark Energy Science Collaboration (DESC) Science Requirements Document v1 Released Data Products (citation below), which are subject to the [Creative Commons Attribution Share Alike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/legalcode) license.

(License for this repo?)

The LSST Dark Energy Science Collaboration, Mandelbaum, Rachel, Eifler, Tim, Hlozek, Renee, Collett, Thomas, Gawiser, Eric, … Troxel, M. A. (2018). The LSST Dark Energy Science Collaboration (DESC) Science Requirements Document v1 Released Data Products (Version 1.0.1) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.2662127

