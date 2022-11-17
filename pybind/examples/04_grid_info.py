import matplotlib.pyplot as plt

import mhm

# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd="../mhm",
)
mhm.run.prepare()  # global_parameters(:, 3) by default
# mhm.run.prepare_domain()
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
