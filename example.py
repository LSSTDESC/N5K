import numpy as np
import time
import n5k

c = n5k.N5KCalculatorCCLNonLimber('tests/config_ccl_nonlimber.yml')
c.setup()
c.run()
c.write_output()
c.teardown()
exit(1)

print([c.__name__ for c in n5k.N5KCalculatorBase.__subclasses__()])
print(n5k.N5KCalculatorCCL.__subclasses__())
# The evaluation of challenge entries will also include accuracy and
# integratbility - this script is just to show an example of timing
# an entry.
exit(1)
def time_run(cls, config):
    c = cls(config)
    c.setup()
    t0 = time.time()
    c.run()
    tf = time.time()
    c.write_output()
    c.teardown()
    print(f'{cls.name}: t={tf-t0}s')

for cls, config in zip([n5k.N5KCalculatorBase,
                        n5k.N5KCalculatorCCL,
                        n5k.N5KCalculatorCCL,
                        n5k.N5KCalculatorCCLNonLimber],
                       ['tests/config.yml', 'tests/config.yml',
                        'tests/config_ccl_kernels.yml',
                        'tests/config_ccl_nonlimber.yml']):
    time_run(cls, config)
