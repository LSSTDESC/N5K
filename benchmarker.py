import sys
sys.path.append("fftlogx/")
import numpy as np
import time
import n5k

conf = sys.argv[1]
name = sys.argv[2]
lmax = int(sys.argv[3])

conf_test = {'output_prefix': 'outputs/tester',
             'benchmark_prefix': 'tests/benchmarks_nl'}
ct = n5k.N5KCalculatorTester(conf_test)
ct.setup()
sn = ct.compare(name, conf)
fn_sn = conf_test['output_prefix']+'_comp_'+name+'.npz'
f_sn = np.load(fn_sn)
ids = f_sn['ls'] <= lmax
sn_lmax = np.sqrt(np.sum(f_sn['sn_per_l'][ids]**2))
print(f"Total S/N = {sn}")
print(f"S/N to ell={lmax} = {sn_lmax}")
