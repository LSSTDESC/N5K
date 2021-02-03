# N5K: Non-local No-Nonsense Non-Limber Numerical Knockout

**Deadline extension:** To allow time to properly address some recently flagged issues with the accuracy benchmarks, we are extending the final challenge deadline to **February 5th, 2021**.

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

Specifically, the non-Limber integral to be computed for each element of the angular power spectrum is,for clustering:
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)j_\ell(k&space;\chi_1)j_\ell(k&space;\chi_2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)j_\ell(k&space;\chi_1)j_\ell(k&space;\chi_2)" title="" /></a>

for cosmic shear:
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\frac{(\ell&space;&plus;2)!}{(\ell-2)!}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)\frac{j_\ell(k&space;\chi_1)}{(k\chi_1)^2}\frac{j_\ell(k&space;\chi_2)}{(k&space;\chi_2)^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\frac{(\ell&space;&plus;2)!}{(\ell-2)!}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)\frac{j_\ell(k&space;\chi_1)}{(k\chi_1)^2}\frac{j_\ell(k&space;\chi_2)}{(k&space;\chi_2)^2}" title="{}" /></a>

and for galaxy-galaxy lensing (assuming quantities labeled with '1' are associated with shear and '2' with number counts):
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\sqrt{\frac{(\ell&space;&plus;2)!}{(\ell-2)!}}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)\frac{j_\ell(k&space;\chi_1)}{(k\chi_1)^2}j_\ell(k&space;\chi_2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{200}&space;C_\ell&space;=&space;\frac{2}{\pi}&space;\sqrt{\frac{(\ell&space;&plus;2)!}{(\ell-2)!}}&space;\int_0^\infty&space;d\chi_1&space;K(\chi_1)&space;\int_0^\infty&space;d\chi_2&space;K(\chi_2)&space;\int_0^\infty&space;dk&space;\,&space;k^2&space;P_\delta(k,z_1,z_2)\frac{j_\ell(k&space;\chi_1)}{(k\chi_1)^2}j_\ell(k&space;\chi_2)" title="" /></a>

where K is in each case the kernel for the appropriate tracer (as provided in `input/kernels.npz`) and <img src="https://render.githubusercontent.com/render/math?math=P_\delta"> is the non-linear matter power spectrum. You should assume that <img src="https://render.githubusercontent.com/render/math?math=P_\delta(k,z_1,z_2) = \sqrt{P_\delta(k,z_1)P_\delta(k,z_2)}">.

Make a pull request to this repository which includes your new subclass as well as a script which creates the conda or virtualenv environment in which your entry can run. Remember, if you need to modify the provided common code, like the base class, make a separate PR!

If you choose to use given dN/dz's instead of the precomputed full kernels, it is your responsibility to ensure other required cosmological factors are correctly computed using the parameters defined in the base class.


## Deadline

To enter, you must have made a draft entry pull request on this repository by **January 15, 2021** (17:00 PST). 

You will then have time to finalise your entry. Final updates to entry pull requests must be made by **February 5th, 2021** (17:00 PST) (extended from January 22nd, 2021).


## FAQ

**Can I participate if I am not a member of DESC?**

Yes, you can, and you can be an author on the paper (we have a special exemption from the Publication Policy for this challenge). 

We request that all non-DESC entrants agree to abide by the [DESC Code of Conduct](https://lsstdesc.org/assets/pdf/policies/LSST_DESC_Professional_Conduct.pdf) throughout their participation in the challenge (as DESC entrants are also expected to do).

**What is the accuracy level I should be aiming for?**

The required accuracy level will be calculated as that which introduces a spurious <img src="https://render.githubusercontent.com/render/math?math=\chi^2"> which is less than 1 for an LSST Y10 3x2pt analysis as defined in the SRD, computed using a Gaussian covariance matrix. The precise tolerable error for a representative set of <img src="https://render.githubusercontent.com/render/math?math=\ell"> will be posted before the challenge deadline for comparison but is not yet available.

**Can I precompute some stuff and put that calculation in my `setup()` method?**

The thing we care about is how long it takes at each step in a sampling process. So you should put in your `run()` method anything that would need to be rerun at each step in e.g. an MCMC chain. Precomputations which are independent of the parameters that would get varied in that process can go in `setup()`.

**How many cores will the calculation be run on for timing?**

The thing we care about is being able to do this calculation at each step in a parameter inference sampling process. This usually means a single core for a single theory calculation, of which the Non-Limber integration forms one step. So assume a single core.

**Can I assume the availability of GPUs?**

If you can convince us that your GPU method can be easily integrated with the rest of the DESC parameter inference pipeline and can be run on GPU resources that DESC will definitely have access to, go for it. Otherwise, no.

