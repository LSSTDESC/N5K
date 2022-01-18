import sys
sys.path.append("fftlogx/")
import numpy as np
import time
import n5k


def time_run(cls, config):
    c = cls(config)
    c.setup()
    niter = 5
    ts = np.zeros(niter)
    for i in range(niter):
        t0 = time.time()
        c.run()
        tf = time.time()
        ts[i] = tf-t0
    tmean = np.mean(ts)
    terr = np.std(ts)/np.sqrt(niter)
    c.write_output()
    c.teardown()
    print('%s: t=(%f+-%f)s'%(cls.name,tmean,terr))
    print(ts)


conf = sys.argv[1]
calc = n5k.n5k_calculator_from_name(sys.argv[2])
time_run(calc, conf)
