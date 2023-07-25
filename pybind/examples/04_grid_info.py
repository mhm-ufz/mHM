from pathlib import Path

import mhm

here = Path(__file__).parent
mhm.model.set_verbosity(level=1)
# use single test domain 1
mhm.model.init(cwd=here / ".." / ".." / "test_domain")
mhm.run.prepare()

print("horizons :", mhm.get.number_of_horizons())
print("ncols, nrows, ncells, xll, yll, cell_size, no_data")
print("L0 :", mhm.get.l0_domain_info(domain=1))
print("L1 :", mhm.get.l1_domain_info(domain=1))
print("L11:", mhm.get.l11_domain_info(domain=1))
print("L2 :", mhm.get.l2_domain_info(domain=1))

print()
shp = mhm.get.l2_domain_shape(domain=1)
msk = mhm.get.l2_domain_mask(*shp, domain=1)
print("L2 :", mhm.get.l2_domain_size(domain=1))
print("L2 :", mhm.get.l2_domain_shape(domain=1))
print("L2 :", mhm.get.l2_domain_mask(*shp, domain=1))
print(msk.shape)

mhm.model.finalize()
