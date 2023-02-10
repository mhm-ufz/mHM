import matplotlib.pyplot as plt

import mhm

# download test domain 1
mhm.download_test(path="example_domain")
# run the downloaded example
mhm.model.init(cwd="example_domain")
mhm.model.run()
runoff = mhm.get_runoff()

# finalize will deallocate all variables.
mhm.model.finalize()

plt.plot(runoff[:, 0])
plt.show()
