import numpy as np
import os
import matplotlib.pyplot as plt 

cltype = 'gg'
for index in range(55):
# index = 10
	bench = np.load('tests/benchmarks_nl_cl%s.npz'%(cltype))

	fang = np.load('tests/ccl_nonlim_fang_cl%s.npz'%(cltype))

	limber = np.load('tests/ccl_limber_cl%s.npz'%(cltype))

	plt.plot(bench['ls'], bench['cls'][index], label='bench')
	plt.plot(fang['ls'], fang['cls'][index], label='fang')
	plt.plot(limber['ls'], limber['cls'][index], label='limber')
	plt.xscale('log')
	# plt.yscale('log')
	plt.legend()
	plt.show()