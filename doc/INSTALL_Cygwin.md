# Cygwin details

[TOC]

[Cygwin](https://cygwin.com/) is an environment with a terminal that allows to compile and
run programs of Unix-like systems. You can find further instructions to install cygwin on the webpage, as well as
instructions on how to install further dependencies after the installation.

After the installation of cygwin and its dependencies mHM will be installed
using cygwin. All commands and the execution of mHM only run in that environment.

Install cygwin by executing the cygwin setup and choose the following dependencies:

- gcc-fortran (the fortran compiler)
- make
- cmake (version >= 3.18)
- libnetcdf-fortran-devel
- libhdf5-devel
- libgfortran
- gfortran

While installing cygwin you will have to choose a mirror. A mirror is a server
on the internet where the files for the installation come from. Choose any server
located near your city and when in doubt, choose the first one in the list.
In the next step you can find all available packages provided by cygwin, set
the view to "full". In the search panel you can filter the packages
by the dependencies listed above (e.g. make). When you choose a
version, the newest one is usually a good choice if not marked as experimental.

*Note for UFZ members:* Install cygwin locally (on your own computer), do not choose a location on the
network for the installation.

Some cygwin versions create a new home directory for you. You may check e.g. here:

    C:\cygwin64\home\$username

As from December 2019, step-by-step guidelines, how to install all cygwin libraries can be viewed in [this youtube video](https://youtu.be/FGJOcYEzbP4)
created by Mehmet Cüneyd Demirel (Istanbul Technical University).

Step-by-step mHM compilation in CYGWIN platform

1) Change directory to mHM folder. Use single quote e.g. 'D:/mhm/' if there is space in the path.
        cd mhm

2) Make a sub-directory inside mHM folder e.g. build
        mkdir build

3) Change directory to build subfolder
        cd build

4) Execute `cmake` with the path to the Git source directory as parameter. (If you followed the instructions above, the path is `..` )
        cmake ..

Other `cmake` options:

To avoid memory issues, allocate stack memory during cmake

    cmake -DCMAKE_Fortran_FLAGS="-Wl,--stack,12485760" ..

Then use cmake option below. Note that memory dump is a common issue for cygwin users when compiling with OpenMP. For memory allocation please also use this line below.

    cmake -DCMAKE_Fortran_FLAGS="${CMAKE_Fortran_FLAGS} -Wl,--stack,12485760" -DMHM_WITH_OpenMP=ON -DCMAKE_BUILD_TYPE=Release ..

4) Execute `make`
        make

5) If all went well then mhm.exe must be created inside the build/app folder. If you are in build folder then copy mhm.exe to upper folder.
        cp ./app/mhm.exe ..

6) Change directory to upper level and then call mhm
        cd ..
        ./mhm examples/domain_01


## Troubleshooting

If libraries are not found, the problem can be:
- you accidentally tried to use commands within the cmd shell and not within the cygwin shell
- your cygwin setup might be broken. Try uninstalling, following the instructions <https://cygwin.com/faq/faq.html#faq.setup.uninstall-service>, and reinstalling again. During
reinstallation only install the required dependencies.
