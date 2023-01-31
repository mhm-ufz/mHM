"""Python bindings for mHM."""
import os
import shutil
import subprocess
import tempfile

import setuptools
import skbuild


class sdist(setuptools.command.sdist.sdist):
    """Custom sdist command to include FORCES in source distribution."""

    def run(self):
        print("## mHM Python setup: adding FORCES to sdist")
        here = os.getcwd()
        forces_dir = os.path.join(here, "forces")
        # remove potentially existing dir
        shutil.rmtree(forces_dir, ignore_errors=True)
        # download forces by calling cmake (CPM will download correct version)
        cmake = skbuild.constants.CMAKE_DEFAULT_EXECUTABLE
        with tempfile.TemporaryDirectory(dir=here) as tmp_dir:
            subprocess.run(
                (cmake, f"-B{tmp_dir}", f"-DCPM_SOURCE_CACHE={here}"),
                capture_output=True,
                cwd=here,
            )
        # CPM stores forces in folder with git-hash as name
        if not os.path.exists(forces_dir) or len(os.listdir(forces_dir)) != 1:
            shutil.rmtree(forces_dir, ignore_errors=True)
            raise RuntimeError("mHM Python sdist: could not download FORCES.")
        hash_dir = os.path.join(forces_dir, os.listdir(forces_dir)[0])
        # copy contents to 'forces' top-folder
        shutil.copytree(
            hash_dir, forces_dir, copy_function=shutil.move, dirs_exist_ok=True
        )
        shutil.rmtree(hash_dir)
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

# env var to control the installation of a console script for mHM
if int(os.getenv("MHM_BUILD_PY_SCRIPT", "1")):
    entry_points = {"console_scripts": ["mhm = mhm.cli:mhm"]}
    cmake_args += ["-DBUILD_MHM_DRIVER=ON"]
    print("## mHM Python setup: creating console script for mHM driver")
else:
    entry_points = {}
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
