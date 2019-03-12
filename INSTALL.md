Dependencies:
=============
For Windows, some Linux distributions and soon also MacOS
specific installation instructions can be found below

- a fortran compiler
- make (a tool to compile a program)
- cmake (version >= 3.5) (a tool to create a system dependent makefile)
- fitting netcdf-fortran libraries (libraries for the usage of the data format netcdf on which mhm depends)
- (optional, but makes things much easier) git

git is a version-control system. If you want to contribute to a project, it is highly recommended to
use git. You can use git to download (clone) the project to your local pc and have a look at the history or
synchronize it without copying the whole repository again. But you can also download the project folder without
git on your local pc.

System dependend installation instructions:
===========================================
### Windows:
cygwin is an environment with a terminal that allows to compile and
run programs of Unix-like systems

https://cygwin.com/

you can find further instructions to install cygwin on the webpage, as well as
instructions on how to install further dependencies after the installation

After the installation of cygwin and its dependencies mHM will be installed
using cygwin. All commands and the execution of mHM only run in that environment.

install cygwin and during the installation choose the dependencies

gcc-fortran (the fortran compiler)

make

cmake (version >= 3.5)

libnetcdf-fortran-devel

git (optional) (you can also install git without cygwin, and it is not required to install git at all, see above.)

While installing cygwin you will have to choose a mirror. A mirror is a server
on the internet where the files for the installation come from. Choose something
located near to you and if you have no idea, choose the first. In the next step you
can find the dependencies in a window appearing during the installation, after
changing the category to full. Here you
can enter the search term (e.g. make) in the search field. You can choose a
version. The newest one, if not marked as experimental, is usually a good choice.

Note for people at the ufz: Install cygwin locally, do not choose a location on the
network for the installation

cygwin sometimes creates a new home directory for you. You may find it in

`C:\cygwin64\home\$username`

Depending on the setup it may also be your home directory.

### Ubuntu, Mint and other apt-get based systems with matching repositories:
`sudo apt-get install git # (optional)`

`sudo apt-get install gfortran netcdf-bin libnetcdf-dev libnetcdff-dev cmake`

### Archlinux:
`sudo pacman -S git # (optional)`

`sudo pacman -S gcc-libs netcdf-fortran cmake`

### Module systems:
load modules gcc, netcdf-fortran, cmake

or load modules intel, netcdf-fortran, cmake

they will have system specific names, environments, etc. Use module spider to find the right packages and the
right dependencies, potentially use corresponding wiki pages

#### on eve (the cluster at the ufz):

from the source directory use a script provided in moduleLoadScripts, for example for the gnu 7.3 compiler:

`source moduleLoadScripts/gfortran73`

### MacOS:
to be added

Specific setups:
================
The hints here are for the replacement of the step "cmake .." in the installation instruction.

You can skip this part and continue with "Installation", if you do not have a module system
setup (like on clusters) or if you have not installed all packages with a package manager
or cygwin

### Module systems:
The executable can be build in a way that it runs independend of loaded modules in the end. The
module system though adds paths in the backround the user should not care about too much, so
the setup is a workaround. (This would be the case with any other building tool aswell.)
It should be stable, anyway. If you want to have a module independend build in the end, instead of
just executing

`cmake ..`

either run

`cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND:STRING=ON ..`

or

`cmake -C ../CMakeCacheFiles/eve ..`

or change the variable `CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND` with ccmake to ON after running cmake ..

### None standard locations for the netcdf-library (e.g. standard setup Macs in CHS):
find, where the nf-config file is located

for example:

`find / -iname "*nf-config*" 2>/dev/null`

(searches the root directory / for a file with a name containing the string "nf-config", not
taking into account upper and lower case, but writes the error messages like "permission denied" into
the void)

then, instead of running `cmake ..` if not using the standard compiler, set the fortran compiler variable to the wished compiler, e.g.

`export FC=gfortran`

then either run

`cmake -DCMAKE_NETCDF_DIR:STRING=/path/to/nf-config/of/used/compiler`

or copy ../CMakeCacheFiles/specificSetup to somewhere

`cp ../CMakeCacheFiles/specificSetup .`

(though, if you want to keep it, rather choose a place outside the build repository), edit it: add
there the path to your nf-config file, and after editing, run

`cmake -C specificSetup ..`

or change the variable `CMAKE_NETCDF_DIR` to the path to the nf-config file with ccmake after running cmake ..

Installation:
=============
change to a directory where you want to store the source code

clone the corresponding mHM repository into a folder, either using git (if installed)

`git clone -b nag_compilation https://git.ufz.de/mhm/mhm.git mhm-nag_compilation/`

for cloning it into a folder `mhm-nag_compilation`,

or download and unpack it using the download link on <https://git.ufz.de/mhm/mhm/tree/nag_compilation>
you can find it, clicking onto the cloud symbol with the arrow in it.

create a build directory where you want to store the build, e.g. inside the git source directory

`cd mhm-nag_compilation`

`mkdir build`

change into the build directory

`cd build`

execute cmake with the path to the git source directory

Note: for specific setups like on
module systems or when the netcdf libraries are not located where the package manager would do so
and they are not saved in environment variables (i.e. classical MacOS setups in CHS) have
a look to "Specific setups", above

`cmake ..`

if everything worked well a Makefile was created with the corresponding paths. Execute make

`make`

if this also worked fine an executable was created, which has to be moved or copied to the git source directory

`cd ..`

`cp build/mhm .`

on Windows it is called mhm.exe. Run instead

`cp build/mhm.exe .`

now you might execute mHM.

`./mhm`

note for the development of the cmake setup: one could automatically link the executable with the cmake code inside the git source directory which is not done for two reasons:
1. it should never be commited, nothing should be build inside the source directory which we did not do by hand
2. the directory where mHM is executed usually is not the source directory but the directory where you want to run
   your tests. In case of the test setup it is the same, usually it is not

Trouble shooting:
=================
On brew/homebrew setup MacOS systems there is no working nf-config by now. Execute

`nf-config --all`

and if it says something like "is not implemented yet" the issue is not solved yet but on my tracklist

in any other case feel free to write an email to <maren.kaluza@ufz.de>

cmake is far from being my main task, so it will probably take a while until I can track a problem. I would
be happy having bug reports, anyhow.
