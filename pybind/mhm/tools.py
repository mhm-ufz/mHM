import numpy as np

from .wrapper import get


def get_runoff():
    shp = get.runoff_shape()
    return get.runoff(*shp)


def get_variable(name, index=1, indexing="xy"):
    name = name.upper()  # convention
    grid = name.split("_")[0].lower()
    if grid not in ["l0", "l1", "l11", "l2"]:
        raise ValueError(f"Unknown variable: {name}")
    n = getattr(get, grid + "_domain_size")()
    shp = getattr(get, grid + "_domain_shape")()
    # mask in mHM is the opposite in numpy
    mask = np.array(getattr(get, grid + "_domain_mask")(m=shp[0], n=shp[1]), dtype=bool)
    # ncols, nrows, ncells, xll, yll, cell_size, no_data
    grid_info = getattr(get, grid + "_domain_info")()
    var = getattr(get, grid + "_variable")(name=name, n=n, idx=index)
    # reshaping
    mask = mask.ravel(order="F")
    output = np.ma.empty_like(mask, dtype=float)
    output.fill_value = grid_info[-1]
    output.mask = ~mask
    output[mask] = var
    output = output.reshape((grid_info[1], grid_info[0]), order="F")
    return output.T if indexing == "ij" else output
