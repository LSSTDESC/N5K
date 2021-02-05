import numpy as np
import cppimport
funcs = cppimport.imp ("wrap")

l = np.arange (10)
print (l)
print (funcs.Cl_auto_gal (l))
