# N5K: Non-local No-Nonsense Non-Limber Numerical Knockout

In this challenge, you are asked to compute a set of power spectra for a 3x2pt (cosmic shear, galaxy-galaxy lensing, and galaxy clustering) analysis without the use of the Limber approximation.

The challenge entries will be evaluated on the basis of accuracy, speed, and integrability with the [Core Cosmology Library](https://github.com/LSSTDESC/CCL/) (CCL). CCL is python (with some heavy-lifting done in C under the hood); a code which is in python or has a python wrapper will satisfy the integrability criteria. Given this, the code which can accomplish the challenge task fastest and within the accuracy requirements of an LSST Y10 cosmological analysis will win the challenge. The winning code will be incorporated for use as the non-Limber integration tool for CCL.

All entrants will have the opportunity to be an author on the resulting paper, which will describe the challenge, compare the different methods, detail the results, and discuss lessons learned.

## How to enter

The challenge asks you to compute the angular spectra required for a 3x2pt analysis setup similar to the LSST Y10 scenario in the [LSST DESC Science Requirements Document v1](https://arxiv.org/pdf/1809.01669.pdf). The 'input' folder of this repo contains some required inputs for this calculation:
- kernels for the 10 number counts tracers and 5 weak lensing tracers as a function of comoving radial distance chi  [N5K/input/kernels.npz](input/kernels.npz).
- Linear and non-linear matter power spectrum as a function of k and z [N5K/input/Pk.npz](input/Pk.npz).
- Background radial comoving distance and normalized expansion rate (H(z)/H(0)), in case you need them: [N5K/input/background.npz](input/background.npz).
- dN/dz's for the 10 number counts tracers and 5 weak lensing tracers as a function of z [N5K/input/dNdzs.npz](input/dNdzs.npz).

For the purposes of the challenge, we ignore intrinsic alignments, redshift-space distortions, and magnification. The kernels are thus number counts and cosmic shear only.

[calculator_base.py](n5k/calculator_base.py) contains a base class `N5KCalculatorBase`. Write a subclass of `N5KCalculatorBase` which contains methods `setup()` (to set up your nonlimber calculation), run `run()` (to run it) and `teardown()` (to clean up any dirt you've created). Only the `run()` method will be timed (but don't be cheeky!). [calculator_ccl.py](n5k/calculator_ccl.py) contains an example of what this would look like doing the calculation using CCL's current (Limber) calculation tools. If you need to modify the challenge machinery (e.g. the base class itself or other provided files), you must make a separate pull request to do this (i.e. don't do it in your challenge entry pull request).

Specifically, the non-Limber integral to be computed for each element of the angular power spectrum is:

<img src="https://render.githubusercontent.com/render/math?math=C_\ell = \frac{2}{\pi} \int_0^\infty d\chi_1 K(\chi_1) \int_0^\infty d\chi_2 K(\chi_2) \int_0^\infty dk \, k^2 P_\delta(k,z_1,z_2)j_\ell(k \chi_1)j_\ell(k \chi_2)">

where K are the kernels and <img src="https://render.githubusercontent.com/render/math?math=P_\delta"> is the non-linear matter power spectrum. You should assume that <img src="https://render.githubusercontent.com/render/math?math=P_\delta(k,z_1,z_2) = \sqrt{P_\delta(k,z_1)P_\delta(k,z_2)}">.

Make a pull request to this repository which includes your new subclass as well as a script which creates the conda or virtualenv environment in which your entry can run. Remember, if you need to modify the provided common code, like the base class, make a separate PR!

If you choose to use given dN/dz's instead of the precomputed full kernels, it is your responsibility to ensure other required cosmological factors are correctly computed using the parameters defined in the base class.


## Deadline

The challenge will close on January 15, 2021.

## FAQ

**Can I participate if I am not a member of DESC?**

Yes, you can, and you can be an author on the paper (we have a special exemption from the Publication Policy for this challenge). 

We request that all non-DESC entrants agree to abide by the [DESC Code of Conduct](https://lsstdesc.org/assets/pdf/policies/LSST_DESC_Professional_Conduct.pdf) throughout their participation in the challenge (as DESC entrants are also expected to do).

**What is the accuracy level I should be aiming for?**

The required accuracy level will be calculated as that which introduces a spurious <img src="https://render.githubusercontent.com/render/math?math=\chi^2"> which is less than 1 for an LSST Y10 3x2pt analysis as defined in the SRD, computed using a Gaussian covariance matrix. The precise tolerable error for a representative set of <img src="https://render.githubusercontent.com/render/math?math=\ell"> will be posted before the challenge deadline for comparison but is not yet available.
