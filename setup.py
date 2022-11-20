"""Python bindings for mHM."""
import os

from skbuild import setup

# maybe overwrite the default version
mhm_build_type = os.getenv("MHM_BUILD_TYPE", "Release")
forces_path = os.getenv("MHM_BUILD_FORCES_PATH", "")
# init cmake args
cmake_args = [
    f"-DCMAKE_BUILD_TYPE={mhm_build_type}",
    "-DBUILD_MHM_PYBIND=ON",
]

print(f"## mHM Python setup: build-type '{mhm_build_type}'")

# you can set MHM_BUILD_PARALLEL=0 or MHM_BUILD_PARALLEL=1
if int(os.getenv("MHM_BUILD_PARALLEL", "0")):
    cmake_args += ["-DCMAKE_WITH_OpenMP=ON"]
    print("## mHM Python setup: OpenMP used by env-var.")
if forces_path:
    cmake_args += [f"-DCPM_forces_SOURCE={forces_path}"]
    print(f"## mHM Python setup: using forces path '{forces_path}'")

setup(
    packages=["mhm"],
    package_dir={"": "pybind"},
    cmake_install_dir="pybind/mhm",
    cmake_args=cmake_args,
    zip_safe=False,
)
