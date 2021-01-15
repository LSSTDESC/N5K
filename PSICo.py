import numpy as np

import WindowFunction
import Survey
import PowerSpectrum
import PreliminaryFunctions
import Cl


wf = WindowFunction.build (par_fname)
survey = Survey.build (par_fname)
lmps = PowerSpectrum.LinearMatterPowerSpectrum.build (par_fname)

cl = np.array ([])
treshold = PreliminaryFunctions.parse ("treshold")
for l in range (lmin, lmax, 1):
	limber = Cl.Limber.build (l, ws, survey, lmps)
	if (need_nonLimber == 1):
		nonLimber = Cl.NonLimber.build (l, ws, survey, lmps)
		cl = np.append (cl, [nonLimber])
	if nonlimber/limber - 1. <= treshold:
		need_nonlimber = 0
		cl = np.append (cl, [Limber])

outfname = PreliminaryFunctions.parse ("outfname")
Cl.save_Cl_to_text_file (outfname, cl)
