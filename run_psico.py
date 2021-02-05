import numpy as np
import time
import n5k

def time_run(cls, config):
	c = cls(config)
	c.setup()
	t0 = time.time()
	c.run()
	tf = time.time()
#	c.write_output()
	c.teardown()
	print(f'{cls.name}: t={tf-t0}s')

for cls, config in zip([n5k.N5KCalculatorPSICo],
                       ['tests/config.yml']):
	time_run(cls, config)
