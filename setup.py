"""Python bindings for mHM."""
import os
import shutil
import tarfile
import tempfile
from pathlib import Path
from urllib.request import urlretrieve

import setuptools
import skbuild

FORCES_URL = "https://git.ufz.de/chs/forces/-/archive/{branch}/forces-{branch}.tar.gz"


def _download_tar(url, path):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tar_file = Path(tmp_dir) / "temp.tar.gz"
        tar_dir = Path(tmp_dir) / "temp"
        urlretrieve(url, tar_file)
        with tarfile.open(tar_file, "r:gz") as tar:
            tar.extractall(path=tar_dir)
        # move sub-folder content to desired path
        first_sub = tar_dir / os.listdir(tar_dir)[0]
        shutil.copytree(first_sub, path, ignore_dangling_symlinks=True)


class sdist(setuptools.command.sdist.sdist):
    """Custom sdist command to include FORCES in source distribution."""

    def run(self):
        print("## mHM Python setup: adding FORCES to sdist")
        here = Path(__file__).parent
        ver_file = here / "version_forces.txt"
        forces_dir = here / "forces"
        forces_ver = ver_file.read_text().strip()
        forces_url = FORCES_URL.format(branch=forces_ver)
        # remove potentially existing dir
        shutil.rmtree(forces_dir, ignore_errors=True)
        # download forces
        _download_tar(url=forces_url, path=forces_dir)
        # run sdist
        super().run()
        # remove forces dir again
        shutil.rmtree(forces_dir)


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

entry_points = {"console_scripts": ["mhm-download = mhm.download:cli"]}
# env var to control the installation of a console script for mHM
if int(os.getenv("MHM_BUILD_PY_SCRIPT", "1")):
    entry_points["console_scripts"].append("mhm = mhm.cli:mhm")
    cmake_args += ["-DBUILD_MHM_DRIVER=ON"]
    print("## mHM Python setup: creating console script for mHM driver")
else:
    cmake_args += ["-DBUILD_MHM_DRIVER=OFF"]
    print("## mHM Python setup: no console script for mHM driver")

skbuild.setup(
    packages=["mhm"],
    package_dir={"": "pybind"},
    cmake_install_dir="pybind/mhm",
    cmake_args=cmake_args,
    zip_safe=False,
    entry_points=entry_points,
    cmdclass={"sdist": sdist},
)
