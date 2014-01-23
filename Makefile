# -*- Makefile -*-
#
# PURPOSE
#     CHS Makefile for Fortran, C and mixed projects
#
# CALLING SEQUENCE
#     make [options] [VARIABLE=VARIABLE ...] [targets]
#
#     Variables can be set on the command line [VAR=VAR] or in the section SWITCHES below.
#
#     If $(PROGNAME) is given, an executable will be compiled.
#     If $(LIBNAME)  is given, a library will be created.
#
#     Sources are in $(SRCPATH), which can be several directories separated by whitespace.
#
#     Fortran 90 file endings: .f90, .F90, .f95, .F95, .f03, .F03, .f08, .F08
#     Fortran 77 file endings: .f,   .F,   .for, .FOR, .f77, .F77, .ftn, .FTN
#     C file endings:          .c,   .C
#
# TARGETS
#     all (default), check (=test), clean, cleanclean, cleancheck (=cleantest=checkclean=testclean),
#     dependencies (=depend), html, pdf, latex, doxygen, info
#
# OPTIONS
#     All make options such as -f makefile. See 'man make'.
#
# VARIABLES
#     All variables defined in this makefile.
#     This makefile has lots of conditional statements depending on variables.
#     If the variable works as a switch then the condition checks for variable = true,
#     i.e. ifeq ($(variable),true)
#     otherwise the variable can have any other value.
#     See individual variables in section SWITCHES below or type 'make info'.
#
#     Variables can be empty for disabling a certain behaviour,
#     e.g. if you do not want to use IMSL, set:  imsl=no  or  imsl=
#
#     For main variables see 'make info'.
#
# DEPENDENCIES
#    This makefile uses the following files:
#        $(MAKEDPATH)/make.d.sh, $(CONFIGPATH)/$(system).$(compiler), $(CONFIGPATH)/$(system).alias
#    The default $(MAKEDPATH) and $(CONFIGPATH) is make.config
#    The makefile can use doxygen for html and pdf automatic documentation. It is then using:
#        $(DOXPATH)/doxygen.config
#    If this is not available, it uses the perl script f2html for html documentation:
#        $(CONFIGPATH)/f2html, $(CONFIGPATH)/f2html.fgenrc
#
# RESTRICTIONS
#    Not all packages work with or are compiled for all compilers.
#    The static switch is maintained like a red-headed stepchild. Libraries might be not ordered correctly
#    if static linking and --begin/end-group is not supported.
#
# EXAMPLE
#    make release=debug compiler=intel11 imsl=vendor mkl=mkl95 PROGNAME=prog
#
# NOTES
#    Further information is given in the README, for example
#    on the repository of the makefile, further reading, how to add a new compiler on a system, or
#    how to add a new system.
#
# LICENSE
#    This file is part of the UFZ makefile project.
#
#    The UFZ makefile project is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    The UFZ makefile project is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with the UFZ makefile project. If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2011-2013 Matthias Cuntz, Juliane Mai, Stephan Thober
#
# Written Matthias Cuntz & Juliane Mai, UFZ Leipzig, Germany, Nov. 2011 - mc (at) macu.de

SHELL = /bin/bash

#
# --- SWITCHES -------------------------------------------------------
#

# . is current directory, .. is parent directory
SRCPATH    := ./src ./lib      # where are the source files; use test_??? to run a test directory
PROGPATH   := .       # where shall be the executable
CONFIGPATH := make.config # where are the $(system).$(compiler) files
MAKEDPATH  := make.config # where is the make.d.sh script
DOXPATH    := .           # where is doxygen.config
CHECKPATH  := .           # path for $(CHECKPATH)/test* and $(CHECKPATH)/check* directories if target is check
#
PROGNAME := mhm # Name of executable
LIBNAME  := #mhm.a #libminpack.a # Name of library
#
# Options
# Systems: eve and personal computers such as mcimac for Matthias Cuntz' iMac; look in $(MAKEDPATH) or type 'make info'
system   := eve
# Compiler: intelX, gnuX, nagX, sunX, where X stands for version number, e.g. intel13;
#   look at $(MAKEDPATH)/$(system).alias for shortcuts or type 'make info'
compiler := intel
# Releases: debug, release
release  := release
# Netcdf versions (Network Common Data Form): netcdf3, netcdf4, [anything else]
netcdf   := netcdf4
# LAPACK (Linear Algebra Pack): true, [anything else]
lapack   :=
# MKL (Intel's Math Kernel Library): mkl, mkl95, [anything else]
mkl      :=
# Proj4 (Cartographic Projections Library): true, [anything else]
proj     := 
# IMSL (IMSL Numerical Libraries): vendor, imsl, [anything else]
imsl     :=
# OpenMP parallelization: true, [anything else]
openmp   :=
# Linking: static, shared, dynamic (last two are equal)
static   := shared

# The Makefile sets the following variables depending on the above options:
# FC, FCFLAGS, F90FLAGS, DEFINES, INCLUDES, LD, LDFLAGS, LIBS
# flags, defines, etc. will be set incremental. They will be initialised with
# the following EXTRA_* variables. This allows for example to set an extra compiler
# option or define a preprocessor variable such as: EXTRA_DEFINES := -DNOGUI -DDPREC
#
#
# Specific notes
# If you encounter during linking error messages such as
#     ... relocation truncated to fit: R_X86_64_PC32 against ...
# then you ran out of memory address space, i.e. some hard-coded numbers in the code got too big.
# Check that you have set the 64-bit addressing model in the F90FLAGS and LDFAGS: -m64
# On *nix systems, you can set the addressing model with -mcmodel=medium (F90FLAGS and LDFAGS) for gfortran and intel.
# Intel might also need -shared-intel at the LDFLAGS, i.e.
#     EXTRA_F90FLAGS := -mcmodel=medium
#     EXTRA_LDFLAGS  := -mcmodel=medium -shared-intel
#
#
# Specific notes on optimisation and debugging
# INTEL optimisation: -fast (=-ipo -O3 -static)
#     -fast             Multifile interprocedure optimization
# INTEL debug: -fpe=0 -fpe-all=0 -no-ftz -ftrapuv
#     -fpe=0 -fpe-all=0  errors on all floating point exceptions except underflow.
#     -no-ftz            catches then also all underflows.
#     -ftrapuv           sets undefined numbers to arbitrary values so that floating point exceptions kick in.
# SUN optimisation: -xipo=2
#     -xipo=n 0 disables interprocedural analysis, 1 enables inlining across source files,
#             2 adds whole-program detection and analysis.
# SUN debug: -ftrap=%all, %none, common, [no%]invalid, [no%]overflow, [no%]underflow, [no%]division, [no%]inexact.
#     -ftrap=%n  Set floating-point trapping mode.
# NAG debug: -C=undefined -C=intovf
#     -C=undefined is also checking 0-strings. Function nonull in UFZ mo_string_utils will stop with error.
#     -C=undefined must be used on all routines, i.e. also on netcdf for example.
#                  This means that all tests do not work which use netcdf and/or lapack.
#     -C=intovf    check integer overflow, which is intentional in UFZ mo_xor4096.
EXTRA_FCFLAGS  :=
EXTRA_F90FLAGS :=
EXTRA_DEFINES  := 
EXTRA_INCLUDES :=
EXTRA_LDFLAGS  :=
EXTRA_LIBS     :=
EXTRA_CFLAGS   :=

#
# --- CHECK 0 ---------------------------------------------------
#

# Check available switches
ifeq (,$(findstring $(release),debug release))
    $(error Error: release '$(release)' not found: must be in 'debug release')
endif

ifneq ($(netcdf),)
    ifeq (,$(findstring $(netcdf),netcdf3 netcdf4))
        $(error Error: netcdf '$(netcdf)' not found: must be in 'netcdf3 netcdf4')
    endif
endif

ifeq (,$(findstring $(static),static shared dynamic))
    $(error Error: static '$(static)' not found: must be in 'static shared dynamic')
endif

#
# --- PATHES ------------------------------------------------
#

# Make absolute pathes from relative pathes - there should be no space nor comment at the end of the next lines
SRCPATH    := $(abspath $(SRCPATH:~%=${HOME}%))
PROGPATH   := $(abspath $(PROGPATH:~%=${HOME}%))
CONFIGPATH := $(abspath $(CONFIGPATH:~%=${HOME}%))
MAKEDPATH  := $(abspath $(MAKEDPATH:~%=${HOME}%))
DOXPATH    := $(abspath $(DOXPATH:~%=${HOME}%))
CHECKPATH  := $(abspath $(CHECKPATH:~%=${HOME}%))

# Program names
# Only Prog or Lib
ifeq (,$(strip $(PROGNAME)))
    ifeq (,$(strip $(LIBNAME)))
        $(error Error: PROGNAME or LIBNAME must be given.)
    else
        islib   := True
        LIBNAME := $(PROGPATH)/$(strip $(LIBNAME))
    endif
else
    ifeq (,$(strip $(LIBNAME)))
        islib    := False
        PROGNAME := $(PROGPATH)/$(strip $(PROGNAME))
    else
        $(error Error: only one of PROGNAME or LIBNAME can be given.)
    endif
endif

MAKEDSCRIPT  := make.d.sh
MAKEDEPSPROG := $(MAKEDPATH)/$(MAKEDSCRIPT)

#
# --- CHECK 1 ---------------------------------------------------
#

systems := $(shell ls -1 $(CONFIGPATH) | sed -e "/$(MAKEDSCRIPT)/d" -e '/f2html/d' | cut -d '.' -f 1 | sort | uniq)
ifeq (,$(findstring $(system),$(systems)))
    $(error Error: system '$(system)' not found: known systems are $(systems))
endif

#
# --- ALIASES ---------------------------------------------------
#

# Include compiler alias on specific systems, e.g. nag for nag53
icompiler := $(compiler)
ALIASINC  := $(CONFIGPATH)/$(system).alias
ifeq (exists, $(shell if [ -f $(ALIASINC) ] ; then echo 'exists' ; fi))
    include $(ALIASINC)
endif

#
# --- CHECK 2 ---------------------------------------------------
#
compilers := $(shell ls -1 $(CONFIGPATH) | sed -e "/$(MAKEDSCRIPT)/d" -e '/f2html/d' -e '/alias/d' | grep $(system) | cut -d '.' -f 2 | sort | uniq)
gnucompilers := $(filter gnu%, $(compilers))
nagcompilers := $(filter nag%, $(compilers))
ifeq (,$(findstring $(icompiler),$(compilers)))
    $(error Error: compiler '$(icompiler)' not found: configured compilers for system $(system) are $(compilers))
endif

#
# --- DEFAULTS ---------------------------------------------------
#

# These variables will be used to compile
FC       :=
FCFLAGS  := $(EXTRA_FCFLAGS)
F90      :=
F90FLAGS := $(EXTRA_F90FLAGS)
CC       :=
CFLAGS   := $(EXTRA_CFLAGS)
DEFINES  := $(EXTRA_DEFINES)
INCLUDES := $(EXTRA_INCLUDES)
# and link, and therefore set below
LD       :=
LDFLAGS  := $(EXTRA_LDFLAGS)
LIBS     := $(EXTRA_LIBS) $(addprefix -L,$(SRCPATH))
AR       := ar
ARFLAGS  := -ru
RANLIB   := ranlib

#
# --- COMPILER / MACHINE SPECIFIC --------------------------------
#

# Set path where all the .mod, .o, etc. files will be written, set before include $(MAKEINC)
OBJPATH := $(addsuffix /.$(strip $(icompiler)).$(strip $(release)), $(SRCPATH))

# Include the individual configuration files
MAKEINC := $(addsuffix /$(system).$(icompiler), $(abspath $(CONFIGPATH:~%=${HOME}%)))
#$(info "MAKEINC: "$(MAKEINC))
ifneq (exists, $(shell if [ -f $(MAKEINC) ] ; then echo 'exists' ; fi))
    $(error Error: '$(MAKEINC)' not found.)
endif
include $(MAKEINC)

# Always use -DCFORTRAN for mixed C and Fortran compilations
DEFINES  += -DCFORTRAN

# Mac OS X is special, there is (almost) no static linking.
# MAC OS X does not work with -rpath. Set DYLD_LIBRARY_PATH if needed.
iOS := $(shell uname -s)
istatic := $(static)
ifneq (,$(findstring $(iOS),Darwin))
    istatic := dynamic
endif

# Start group for cyclic search in static linking
iLIBS :=
ifeq ($(istatic),static)
    iLIBS += -Bstatic -Wl,--start-group
else
    ifneq (,$(findstring $(iOS),Darwin))
        iLIBS += -Wl,-dynamic
    else
        iLIBS += -Bdynamic
    endif
endif

# --- COMPILER ---------------------------------------------------
ifneq (,$(findstring $(icompiler),$(gnucompilers)))
    ifneq (exists, $(shell if [ -d "$(GFORTRANDIR)" ] ; then echo 'exists' ; fi))
        $(error Error: GFORTRAN path '$(GFORTRANDIR)' not found.)
    endif
    GFORTRANLIB ?= $(GFORTRANDIR)/lib
    iLIBS       += -L$(GFORTRANLIB) -lgfortran
    RPATH       += -Wl,-rpath,$(GFORTRANLIB)
endif

# --- IMSL ---------------------------------------------------
ifneq (,$(findstring $(imsl),vendor imsl))
    ifneq (exists, $(shell if [ -d "$(IMSLDIR)" ] ; then echo 'exists' ; fi))
        $(error Error: IMSL path '$(IMSLDIR)' not found.)
    endif
    IMSLINC ?= $(IMSLDIR)/include
    IMSLLIB ?= $(IMSLDIR)/lib

    INCLUDES += -I$(IMSLINC)
    ifneq ($(ABSOFT),)
        INCLUDES += -p $(IMSLINC)
    endif
    DEFINES  += -DIMSL

    ifeq (,$(findstring $(iOS),Darwin))
        iLIBS     += -z muldefs
        ifneq ($(istatic),static)
            iLIBS += -i_dynamic
        endif
    endif

    ifneq (,$(findstring $(iOS),Darwin))
        iLIBS += -L$(IMSLLIB) -limsl -limslscalar -limsllapack -limslblas -limsls_err -limslmpistub -limslsuperlu
    else
        ifeq ($(imsl),imsl)
            iLIBS += -L$(IMSLLIB) -limsl -limslscalar -limsllapack_imsl -limslblas_imsl -limsls_err -limslmpistub -limslsuperlu
        else
            iLIBS += -L$(IMSLLIB) -limsl -limslscalar -limsllapack_vendor -limslblas_vendor -limsls_err -limslmpistub -limslsuperlu -limslhpc_l
        endif
    endif
    RPATH += -Wl,-rpath,$(IMSLLIB)
endif

# --- OPENMP ---------------------------------------------------
iopenmp=
ifeq ($(openmp),true)
    ifneq (,$(findstring $(icompiler),$(gnucompilers)))
        iopenmp := -fopenmp
    else
        iopenmp := -openmp
    endif
    DEFINES += -DOPENMP
endif
F90FLAGS += $(iopenmp)
FCFLAGS  += $(iopenmp)
CFLAGS   += $(iopenmp)
LDFLAGS  += $(iopenmp)
# IMSL needs openmp during linking in any case
ifneq ($(openmp),true)
    ifneq (,$(findstring $(imsl),vendor imsl))
        ifneq (,$(findstring $(icompiler),$(gnucompilers)))
            LDFLAGS += -fopenmp
        else
            LDFLAGS += -openmp
        endif
    endif
endif

# --- MKL ---------------------------------------------------
ifneq (,$(findstring $(mkl),mkl mkl95))
    ifeq ($(mkl),mkl95) # First mkl95 then mkl for .mod files other then intel
        ifneq (exists, $(shell if [ -d "$(MKL95DIR)" ] ; then echo 'exists' ; fi))
            $(error Error: MKL95 path '$(MKL95DIR)' not found.)
        endif
        MKL95INC ?= $(MKL95DIR)/include
        MKL95LIB ?= $(MKL95DIR)/lib

        INCLUDES += -I$(MKL95INC)
        ifneq ($(ABSOFT),)
            INCLUDES += -p $(MKL95INC)
        endif
        DEFINES  += -DMKL95

        iLIBS += -L$(MKL95LIB) -lmkl_blas95_lp64 -lmkl_lapack95_lp64
        RPATH += -Wl,-rpath,$(MKL95LIB)
        ifneq ($(ABSOFT),)
            F90FLAGS += -p $(MKL95INC)
        endif
    endif

    ifneq (exists, $(shell if [ -d "$(MKLDIR)" ] ; then echo 'exists' ; fi))
        $(error Error: MKL path '$(MKLDIR)' not found.)
    endif
    MKLINC ?= $(MKLDIR)/include
    MKLLIB ?= $(MKLDIR)/lib

    INCLUDES += -I$(MKLINC)
    ifneq ($(ABSOFT),)
        INCLUDES += -p $(MKLINC)
    endif
    DEFINES  += -DMKL

    iLIBS += -L$(MKLLIB) -lmkl_intel_lp64 -lmkl_core #-lpthread
    ifneq (,$(findstring $(imsl),vendor imsl))
       iLIBS += -lmkl_intel_thread #-lpthread
    else
        ifeq ($(openmp),true)
            iLIBS += -lmkl_intel_thread #-lpthread
        else
            iLIBS += -lmkl_sequential #-lpthread
        endif
    endif
    RPATH += -Wl,-rpath,$(MKLLIB)
endif

# --- NETCDF ---------------------------------------------------
ifneq (,$(findstring $(netcdf),netcdf3 netcdf4))
    ifneq (exists, $(shell if [ -d "$(NCDIR)" ] ; then echo 'exists' ; fi))
        $(error Error: NETCDF path '$(NCDIR)' not found.)
    endif
    NCINC ?= $(strip $(NCDIR))/include
    NCLIB ?= $(strip $(NCDIR))/lib

    INCLUDES += -I$(NCINC)
    ifneq ($(ABSOFT),)
        INCLUDES += -p $(NCINC)
    endif
    DEFINES += -DNETCDF

    iLIBS += -L$(NCLIB)
    RPATH += -Wl,-rpath,$(NCLIB)
    ifeq (libnetcdff, $(shell ls $(NCLIB)/libnetcdff.* 2> /dev/null | sed -n '1p' | sed -e 's/.*\(libnetcdff\)/\1/' -e 's/\(libnetcdff\).*/\1/'))
        iLIBS += -lnetcdff
    endif
    iLIBS += -lnetcdf

    ifeq (exists, $(shell if [ -d "$(NCDIR2)" ] ; then echo 'exists' ; fi))
        NCINC2 ?= $(strip $(NCDIR2))/include
        NCLIB2 ?= $(strip $(NCDIR2))/lib

        INCLUDES += -I$(NCINC2)
        ifneq ($(ABSOFT),)
            INCLUDES += -p $(NCINC2)
        endif

        iLIBS += -L$(NCLIB2)
        RPATH += -Wl,-rpath,$(NCLIB2)
        ifeq (libnetcdff, $(shell ls $(NCLIB2)/libnetcdff.* 2> /dev/null | sed -n '1p' | sed -e 's/.*\(libnetcdff\)/\1/' -e 's/\(libnetcdff\).*/\1/'))
            iLIBS += -lnetcdff
        endif
    endif

    # other libraries for netcdf4, ignored for netcdf3
    ifeq ($(netcdf),netcdf4)
        iLIBS += -L$(HDF5LIB) -lhdf5_hl -lhdf5 -L$(SZLIB) -lsz -lz
        RPATH += -Wl,-rpath,$(SZLIB) -Wl,-rpath,$(HDF5LIB)
        ifneq ($(CURLLIB),)
            iLIBS += -L$(CURLLIB) -lcurl
            RPATH += -Wl,-rpath,$(CURLLIB)
        endif
   endif
endif

# --- PROJ --------------------------------------------------
ifeq ($(proj),true)
    ifneq (exists, $(shell if [ -d "$(PROJ4DIR)" ] ; then echo 'exists' ; fi))
        $(error Error: PROJ4 path '$(PROJ4DIR)' not found.)
    endif
    PROJ4LIB ?= $(PROJ4DIR)/lib
    iLIBS    += -L$(PROJ4LIB) -lproj
    RPATH    += -Wl,-rpath=$(PROJ4LIB)

    ifneq (exists, $(shell if [ -d "$(FPROJDIR)" ] ; then echo 'exists' ; fi))
        $(error Error: FPROJ path '$(FPROJDIR)' not found.)
    endif
    FPROJINC ?= $(FPROJDIR)/include
    FPROJLIB ?= $(FPROJDIR)/lib

    INCLUDES += -I$(FPROJINC)
    ifneq ($(ABSOFT),)
        INCLUDES += -p $(FPROJINC)
    endif
    DEFINES  += -DFPROJ
    iLIBS    += -L$(FPROJLIB) -lfproj4 $(FPROJLIB)/proj4.o
    RPATH    += -Wl,-rpath,$(FPROJLIB)
endif

# --- LAPACK ---------------------------------------------------
ifeq ($(lapack),true)
    # Mac OS X uses frameworks
    ifneq (,$(findstring $(iOS),Darwin))
        iLIBS += -framework veclib
    else
        ifneq (exists, $(shell if [ -d "$(LAPACKDIR)" ] ; then echo 'exists' ; fi))
            $(error Error: LAPACK path '$(LAPACKDIR)' not found.)
        endif
        LAPACKLIB ?= $(LAPACKDIR)/lib
        iLIBS     += -L$(LAPACKLIB) -lblas -llapack
        RPATH     += -Wl,-rpath,$(LAPACKLIB)
    endif
    DEFINES += -DLAPACK
endif

# --- DOXYGEN ---------------------------------------------------
ISDOX := True
ifneq (,$(filter doxygen html latex pdf, $(MAKECMDGOALS)))
    ifeq (exists, $(shell if [ -f $(DOXPATH)/"doxygen.config" ] ; then echo 'exists' ; fi))
        ifneq ($(DOXYGENDIR),)
            ifneq (exists, $(shell if [ -f $(strip $(DOXYGENDIR))/"doxygen" ] ; then echo 'exists' ; fi))
                $(error Error: doxygen not found in $(strip $(DOXYGENDIR)).)
            else
                DOXYGEN := $(strip $(DOXYGENDIR))/"doxygen"
            endif
        else
            ifneq (, $(shell which doxygen))
                DOXYGEN := doxygen
            else
                $(error Error: doxygen not found in $PATH.)
            endif
        endif
        ifneq ($(DOTDIR),)
            ifneq (exists, $(shell if [ -f $(strip $(DOTDIR))/"dot" ] ; then echo 'exists' ; fi))
                $(error Error: dot not found in $(strip $(DOTDIR)).)
            else
                DOTPATH := $(strip $(DOTDIR))
            endif
        else
            ifneq (, $(shell which dot))
                DOTPATH := $(dir $(shell which dot))
            else
                DOTPATH :=
            endif
        endif
        ifneq ($(TEXDIR),)
            ifneq (exists, $(shell if [ -f $(strip $(TEXDIR))/"latex" ] ; then echo 'exists' ; fi))
                $(error Error: latex not found in $(strip $(TEXDIR)).)
            else
                TEXPATH := $(strip $(TEXDIR))
            endif
        else
            ifneq (, $(shell which latex))
                TEXPATH := $(dir $(shell which latex))
            else
                $(error Error: latex not found in $PATH.)
            endif
        endif
        ifneq ($(PERLDIR),)
            ifneq (exists, $(shell if [ -f $(strip $(PERLDIR))/"perl" ] ; then echo 'exists' ; fi))
                $(error Error: perl not found in $(strip $(PERLDIR)).)
            else
                PERLPATH := $(strip $(PERLDIR))
            endif
        else
            ifneq (, $(shell which perl))
                PERLPATH := $(dir $(shell which perl))
            else
                $(error Error: perl not found in $PATH.)
            endif
        endif
    else
        ISDOX += False
        ifneq (,$(filter doxygen latex pdf, $(MAKECMDGOALS)))
            $(error Error: no doxygen.config found in $(DOXPATH).)
        endif
    endif
endif

#
# --- FINISH SETUP ---------------------------------------------------
#

ifeq ($(release),debug)
    DEFINES += -DDEBUG
endif

# Mac OS X is special, there is (almost) no static linking; otherwise close static group
ifeq ($(istatic),static)
    iLIBS += -Wl,--end-group
endif

# The NAG compiler links via gcc so that one has to give -Wl twice and double commas for the option
# i.e. instead of  -Wl,rpath,/path   you need   -Wl,-Wl,,rpath,,/path
ifneq (,$(findstring $(icompiler),$(nagcompilers)))
    comma  := ,
    iiLIBS := $(subst -Wl,-Wl$(comma)-Wl,$(subst $(comma),$(comma)$(comma),$(iLIBS)))
    iRPATH := $(subst -Wl,-Wl$(comma)-Wl,$(subst $(comma),$(comma)$(comma),$(RPATH)))
else
    iiLIBS := $(iLIBS)
    iRPATH := $(RPATH)
endif
LIBS += $(iiLIBS)
# Only Linux and Solaris can use -rpath in executable
ifeq (,$(findstring $(iOS),Darwin))
    LIBS += $(iRPATH)
endif

# some targets should not compiler the code first, e.g. producing documentation
# but some targets should not recompile but be aware of the source files, e.g. clean
iphony    := False
iphonyall := False
ifneq (,$(strip $(MAKECMDGOALS)))
    ifneq (,$(findstring /$(strip $(MAKECMDGOALS))/,/check/ /test/ /html/ /latex/ /pdf/ /doxygen/))
        iphony := True
    endif
    ifneq (,$(findstring $(strip $(MAKECMDGOALS))/,/check/ /test/ /html/ /latex/ /pdf/ /doxygen/ /cleancheck/ /cleantest/ /checkclean/ /testclean/ /info/ /clean/ /cleanclean/))
        iphonyall := True
    endif
endif


# ASRCS contain Fortran 90 source dir informations
ifeq (False,$(iphony))
    SRCS := $(wildcard $(addsuffix /*.f90, $(SRCPATH))) $(wildcard $(addsuffix /*.F90, $(SRCPATH)))  $(wildcard $(addsuffix /*.f95, $(SRCPATH))) $(wildcard $(addsuffix /*.F95, $(SRCPATH))) $(wildcard $(addsuffix /*.f03, $(SRCPATH))) $(wildcard $(addsuffix /*.F03, $(SRCPATH))) $(wildcard $(addsuffix /*.f08, $(SRCPATH))) $(wildcard $(addsuffix /*.F08, $(SRCPATH)))
endif
# source files but all with .o
OSRCS := $(patsubst %.f90,%.o,$(patsubst %.F90,%.o,$(patsubst %.f95,%.o,$(patsubst %.F95,%.o,$(patsubst %.f03,%.o,$(patsubst %.F03,%.o,$(patsubst %.f08,%.o,$(patsubst %.F08,%.o,$(SRCS)))))))))
# object files
OBJS  := $(join $(dir $(OSRCS)), $(addprefix .$(strip $(icompiler)).$(strip $(release))/,$(notdir $(OSRCS))))
# dependency files
DOBJS := $(OBJS:.o=.d)
# g90 debug files of NAG compiler are in current directory or in source directory
GOBJS := $(addprefix $(CURDIR)/,$(patsubst %.o,%.g90,$(notdir $(OBJS)))) $(patsubst %.o,%.g90,$(OSRCS))


# Same for Fortran77 files
ifeq (False,$(iphony))
    FSRCS := $(wildcard $(addsuffix /*.f, $(SRCPATH))) $(wildcard $(addsuffix /*.F, $(SRCPATH)))  $(wildcard $(addsuffix /*.for, $(SRCPATH))) $(wildcard $(addsuffix /*.FOR, $(SRCPATH))) $(wildcard $(addsuffix /*.f77, $(SRCPATH))) $(wildcard $(addsuffix /*.F77, $(SRCPATH))) $(wildcard $(addsuffix /*.ftn, $(SRCPATH))) $(wildcard $(addsuffix /*.FTN, $(SRCPATH)))
endif
# source files but all with .o
FOSRCS := $(patsubst %.f,%.o,$(patsubst %.F,%.o,$(patsubst %.for,%.o,$(patsubst %.FOR,%.o,$(patsubst %.f77,%.o,$(patsubst %.F77,%.o,$(patsubst %.ftn,%.o,$(patsubst %.FTN,%.o,$(FSRCS)))))))))
# object files
FOBJS  := $(join $(dir $(FOSRCS)), $(addprefix .$(strip $(icompiler)).$(strip $(release))/,$(notdir $(FOSRCS))))
# dependency files
FDOBJS := $(FOBJS:.o=.d)
# g90 debug files of NAG compiler are in current directory or in source directory
FGOBJS := $(addprefix $(CURDIR)/,$(patsubst %.o,%.g90,$(notdir $(FOBJS)))) $(patsubst %.o,%.g90,$(FOSRCS))


# Same for C files with ending .c
ifeq (False,$(iphony))
    CSRCS := $(wildcard $(addsuffix /*.c, $(SRCPATH))) $(wildcard $(addsuffix /*.C, $(SRCPATH)))
endif
COSRCS := $(patsubst %.c,%.o,$(patsubst %.C,%.o,$(CSRCS)))
# object files
COBJS  := $(join $(dir $(COSRCS)), $(addprefix .$(strip $(icompiler)).$(strip $(release))/,$(notdir $(COSRCS))))
# dependency files
CDOBJS := $(COBJS:.o=.d)

# Libraries in source path
ifeq (False,$(iphony))
    LSRCS := $(wildcard $(addsuffix /*.a, $(SRCPATH))) $(wildcard $(addsuffix /*.so, $(SRCPATH)))  $(wildcard $(addsuffix /*.dylib, $(SRCPATH)))
endif
LOSRCS := $(patsubst %.a,,$(patsubst %.so,,$(patsubst %.dylib,,$(LSRCS))))
LOBJS  := $(addprefix -L,$(dir $(SRCPATH))) $(addprefix -l, $(patsubst lib%, %, $(notdir $(LOSRCS))))

# The Absoft compiler needs that ABSOFT is set to the Absoft base path
ifneq ($(ABSOFT),)
    export ABSOFT
endif
ifneq ($(LDPATH),)
    empty:=
    space:= $(empty) $(empty)
    export LD_LIBRARY_PATH=$(subst $(space),$(empty),$(LDPATH))
endif

# Link with the fortran compiler if fortran code
ifneq ($(SRCS)$(FSRCS),)
    LD := $(F90)
else
    LD := $(CC)
endif

INCLUDES += $(addprefix -I,$(OBJPATH))

#
# --- TARGETS ---------------------------------------------------
#

#.SUFFIXES: .f90 .F90 .f95 .F95 .f03 .F03 .f08 .F08 .f .F .for .FOR .ftn .FTN .c .C .d .o .a .so .dylib
.SUFFIXES:

.PHONY: clean cleanclean cleantest checkclean testclean cleancheck html latex pdf doxygen check test info

all: $(PROGNAME) $(LIBNAME)

# Link program
$(PROGNAME): $(OBJS) $(FOBJS) $(COBJS)
	@echo "Linking program"
	$(LD) $(LDFLAGS) -o $(PROGNAME) $(OBJS) $(FOBJS) $(COBJS) $(LIBS) $(LOBJS)

# Link library
$(LIBNAME): $(DOBJS) $(FDOBJS) $(CDOBJS) $(OBJS) $(FOBJS) $(COBJS)
	@echo "Linking library"
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJS) $(FOBJS) $(COBJS)
	$(RANLIB) $(LIBNAME)

# Get dependencies
$(DOBJS):
	@dirname $@ | xargs mkdir -p 2>/dev/null
	@nobj=$$(echo $(DOBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(SRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	$(MAKEDEPSPROG) $$src .$(strip $(icompiler)).$(strip $(release)) $(SRCS) $(FSRCS)

$(FDOBJS):
	@dirname $@ | xargs mkdir -p 2>/dev/null
	@nobj=$$(echo $(FDOBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(FSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	obj=$$(echo $(FOBJS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	dobj=$$(echo $(FOBJS) | tr ' ' '\n' | sed -n $${nobj}p | sed 's|\.o[[:blank:]]*$$|.d|') ; \
	echo "$$obj $$dobj : $$src" > $$dobj

$(CDOBJS):
	@dirname $@ | xargs mkdir -p 2>/dev/null
	@nobj=$$(echo $(CDOBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(CSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	pobj=$$(dirname $@) ; psrc=$$(dirname $$src) ; \
	gcc $(DEFINES) -MM $$src | sed "s|.*:|$(patsubst %.d,%.o,$@) $@ :|" > $@

# Compile
$(OBJS):
ifneq (,$(findstring $(icompiler),gnu41 gnu42))
	@nobj=$$(echo $(OBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(SRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	echo $(F90) -E -x c $(DEFINES) $(INCLUDES) $(F90FLAGS) $${src} > .tmp.gf3 ; \
	$(F90) -E -x c $(DEFINES) $(INCLUDES) $(F90FLAGS) $${src} > .tmp.gf3
	$(F90) $(DEFINES) $(INCLUDES) $(F90FLAGS) $(MODFLAG)$(dir $@) -c .tmp.gf3 -o $@
	rm .tmp.gf3
else
	@nobj=$$(echo $(OBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(SRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	echo $(F90) $(DEFINES) $(INCLUDES) $(F90FLAGS) $(MODFLAG)$(dir $@) -c $${src} -o $@ ; \
	$(F90) $(DEFINES) $(INCLUDES) $(F90FLAGS) $(MODFLAG)$(dir $@) -c $${src} -o $@
endif

$(FOBJS):
ifneq (,$(findstring $(icompiler),gnu41 gnu42))
	@nobj=$$(echo $(FOBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(FSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	echo $(FC) -E -x c $(DEFINES) $(INCLUDES) $(FCFLAGS) $$src > .tmp.gf3 ; \
	$(FC) -E -x c $(DEFINES) $(INCLUDES) $(FCFLAGS) $$src > .tmp.gf3
	$(FC) $(DEFINES) $(INCLUDES) $(FCFLAGS) -c .tmp.gf3 -o $@
	rm .tmp.gf3
else
	@nobj=$$(echo $(FOBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(FSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	echo $(FC) $(DEFINES) $(INCLUDES) $(FCFLAGS) -c $$src -o $@ ; \
	$(FC) $(DEFINES) $(INCLUDES) $(FCFLAGS) -c $$src -o $@
endif

$(COBJS):
	@nobj=$$(echo $(COBJS) | tr ' ' '\n' | grep -n -w -F $@ | sed 's/:.*//') ; \
	src=$$(echo $(CSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	echo $(CC) $(DEFINES) $(INCLUDES) $(CFLAGS) -c $${src} -o $@ ; \
	$(CC) $(DEFINES) $(INCLUDES) $(CFLAGS) -c $${src} -o $@

# Helper Targets
clean:
	rm -f $(DOBJS) $(FDOBJS) $(CDOBJS) $(OBJS) $(FOBJS) $(COBJS) $(addsuffix /*.mod, $(OBJPATH))
ifeq (False,$(islib))
	rm -f "$(PROGNAME)"
endif
	rm -f $(GOBJS) $(FGOBJS)
	rm -f *make_check_test_file

cleanclean: clean
	rm -rf $(addsuffix /.*.r*, $(SRCPATH)) $(addsuffix /.*.d*, $(SRCPATH))
	rm -rf "$(PROGNAME)".dSYM $(addsuffix /html, $(SRCPATH))
	@if [ -f $(DOXPATH)/"doxygen.config" ] ; then rm -rf $(PROGPATH)/latex ; fi
	@if [ -f $(DOXPATH)/"doxygen.config" ] ; then rm -rf $(PROGPATH)/html ; fi
ifeq (True,$(islib))
	rm -f "$(LIBNAME)"
endif

cleancheck:
	for i in $(shell ls -d $(CHECKPATH)/test* $(CHECKPATH)/check*) ; do \
	    $(MAKE) SRCPATH=$$i cleanclean ; \
	done

cleantest: cleancheck

checkclean: cleancheck

testclean: cleancheck

check:
ifeq (True,$(islib))
	$(error Error: check and test must be done with PROGNAME not LIBNAME.)
endif
	for i in $(shell ls -d $(CHECKPATH)/test* $(CHECKPATH)/check*) ; do \
	    rm -f "$(PROGNAME)" ; \
	    j=$${i/minpack/maxpack} ; \
	    ldextra= ; \
	    if [ $$i != $$j ] ; then \
	    	 $(MAKE) -s MAKEDPATH=$(MAKEDPATH) SRCPATH="$$i"/../../minpack PROGPATH=$(PROGPATH) \
	    	      CONFIGPATH=$(CONFIGPATH) PROGNAME= LIBNAME=libminpack.a system=$(system) \
	    	      release=$(release) netcdf=$(netcdf) static=$(static) proj=$(proj) \
	    	      imsl=$(imsl) mkl=$(mkl) lapack=$(lapack) compiler=$(compiler) \
	    	      openmp=$(openmp) > /dev/null ; \
                 ldextra="-L. -lminpack" ; \
	    fi ; \
	    $(MAKE) -s MAKEDPATH=$(MAKEDPATH) SRCPATH=$$i PROGPATH=$(PROGPATH) \
	         CONFIGPATH=$(CONFIGPATH) PROGNAME=$(PROGNAME) system=$(system) \
	         release=$(release) netcdf=$(netcdf) static=$(static) proj=$(proj) \
	         imsl=$(imsl) mkl=$(mkl) lapack=$(lapack) compiler=$(compiler) \
	         openmp=$(openmp) EXTRA_LDFLAGS="$$ldextra" > /dev/null \
	    && { $(PROGNAME) 2>&1 | grep -E '(o\.k\.|failed)' ;} ; status=$$? ; \
	    if [ $$status != 0 ] ; then echo "$$i failed!" ; fi ; \
	    $(MAKE) -s SRCPATH=$$i cleanclean ; \
	    if [ $$i != $$j ] ; then \
	    	 $(MAKE) -s SRCPATH="$$i"/../../minpack PROGNAME= LIBNAME=libminpack.a cleanclean ; \
	    fi ; \
	done

test: check

depend: dependencies

dependencies:
	@for i in $(DOBJS) ; do \
	    nobj=$$(echo $(DOBJS) | tr ' ' '\n' | grep -n -w -F $${i} | sed 's/:.*//') ; \
	    src=$$(echo $(SRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    obj=$$(echo $(OBJS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    if [ $${src} -nt $${obj} ] ; then rm $${i} ; fi ; \
	done
	@for i in $(FDOBJS) ; do \
	    nobj=$$(echo $(FDOBJS) | tr ' ' '\n' | grep -n -w -F $${i} | sed 's/:.*//') ; \
	    src=$$(echo $(FSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    obj=$$(echo $(FOBJS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    if [ $${src} -nt $${obj} ] ; then rm $${i} ; fi ; \
	done
	@for i in $(CDOBJS) ; do \
	    nobj=$$(echo $(CDOBJS) | tr ' ' '\n' | grep -n -w -F $${i} | sed 's/:.*//') ; \
	    src=$$(echo $(CSRCS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    obj=$$(echo $(COBJS) | tr ' ' '\n' | sed -n $${nobj}p) ; \
	    if [ $${src} -nt $${obj} ] ; then rm $${i} ; fi ; \
	done
	@rm -f $(addsuffix /$(MAKEDSCRIPT).dict, $(OBJPATH))

doxygen:
	@cat $(DOXPATH)/"doxygen.config" | \
	     sed -e "/^PERL_PATH/s|=.*|=$(PERLPATH)|" | \
	     sed -e "/^DOT_PATH/s|=.*|=$(DOTPATH)|" | env PATH=${PATH}:$(TEXPATH) $(DOXYGEN) -

html:
	@if [ $(ISDOX) == True ] ; then \
	    cat $(DOXPATH)/"doxygen.config" | \
	        sed -e "/^PERL_PATH/s|=.*|=$(PERLPATH)|" | \
	        sed -e "/^DOT_PATH/s|=.*|=$(DOTPATH)|" | env PATH=${PATH}:$(TEXPATH) $(DOXYGEN) - ; \
	else \
	    for i in $(SRCPATH) ; do \
	        $(CONFIGPATH)/f2html -f $(CONFIGPATH)/f2html.fgenrc -d $$i/html $$i ; \
	    done ; \
	fi

latex: pdf

pdf: doxygen
	@cd latex ; env PATH=${PATH}:$(TEXPATH) make pdf

info:
	@echo "CHS Makefile"
	@echo ""
	@echo "Config"
	@echo "system   = $(system)"
	@echo "compiler = $(compiler)"
	@echo "release  = $(release)"
	@echo "netcdf   = $(netcdf)"
	@echo "lapack   = $(lapack)"
	@echo "mkl      = $(mkl)"
	@echo "proj     = $(proj)"
	@echo "imsl     = $(imsl)"
	@echo "openmp   = $(openmp)"
	@echo "static   = $(static)"
	@echo ""
	@echo "Files/Pathes"
	@echo "SRCPATH    = $(SRCPATH)"
	@echo "PROGPATH   = $(PROGPATH)"
	@echo "CONFIGPATH = $(CONFIGPATH)"
	@echo "MAKEDPATH  = $(MAKEDPATH)"
	@echo "CHECKPATH  = $(CHECKPATH)"
	@echo "PROGNAME   = $(basename $(PROGNAME))"
	@echo "LIBNAME    = $(basename $(LIBNAME))"
	@echo "FILES      = $(SRCS) $(FORSRCS) $(CSRCS) $(LASRCS)"
	@echo ""
	@echo "Programs/Flags"
	@echo "FC        = $(FC)"
	@echo "FCFLAGS   = $(FCFLAGS)"
	@echo "F90       = $(F90)"
	@echo "F90FLAGS  = $(F90FLAGS)"
	@echo "CC        = $(CC)"
	@echo "CFLAGS    = $(CFLAGS)"
	@echo "DEFINES   = $(DEFINES)"
	@echo "INCLUDES  = $(INCLUDES)"
	@echo "LD        = $(LD)"
	@echo "LDFLAGS   = $(LDFLAGS)"
	@echo "LIBS      = $(LIBS)"
	@echo "AR        = $(AR)"
	@echo "ARFLAGS   = $(ARFLAGS)"
	@echo "RANLIB    = $(RANLIB)"
	@echo ""
	@echo "Configured compilers on $(system): $(compilers)"
ifeq (exists, $(shell if [ -f $(ALIASINC) ] ; then echo 'exists' ; fi))
	@echo ""
	@echo "Compiler aliases for $(system)"
	@sed -n '/ifneq (,$$(findstring $$(compiler)/,/endif/p' $(ALIASINC) | \
	 sed -e '/endif/d' -e 's/icompiler ://' | \
	 sed -e 's/ifneq (,$$(findstring $$(compiler),//' -e 's/))//' | \
	 paste - - | tr -d '\t' | tr -s ' '
endif
	@echo ""
	@echo "Targets"
	@echo "all (default)  Compile program or library"
	@echo "check          Run all checks in $(CHECKPATH)/test* and $(CHECKPATH)/check* directories"
	@echo "checkclean     alias for cleancheck"
	@echo "clean          Clean compilation of current compiler and release"
	@echo "cleancheck     Cleanclean all test directories $(CHECKPATH)/test* and $(CHECKPATH)/check*"
	@echo "cleanclean     Clean compilations of all compilers and releases"
	@echo "cleantest      alias for cleancheck"
	@echo "depend         alias for dependencies"
	@echo "dependencies   Redo dependencies"
	@echo "doxygen        Run doxygen html with $(DOXPATH)/doxygen.config"
	@echo "html           Run either doxygen html with $(DOXPATH)/doxygen.config or f2html of Jean-Marc Beroud"
	@echo "info           Prints info about current settings and possibilities"
	@echo "latex          alias for pdf"
	@echo "pdf            Run doxygen PDF with $(DOXPATH)/doxygen.config"
	@echo "test           alias for check"
	@echo "testclean      alias for cleancheck"
	@echo ""
	@echo "All possibilities"
	@echo "system      $(shell ls -1 $(CONFIGPATH) | sed -e '/"$(MAKEDSCRIPT)"/d' -e '/f2html/d' | cut -d '.' -f 1 | sort | uniq)"
	@echo "compiler    $(shell ls -1 $(CONFIGPATH) | sed -e '/"$(MAKEDSCRIPT)"/d' -e '/f2html/d' -e '/alias/d' | cut -d '.' -f 2 | sort | uniq)"
	@echo "release     debug release"
	@echo "netcdf      netcdf3 netcdf4 [anything else]"
	@echo "lapack      true [anything else]"
	@echo "mkl         mkl mkl95 [anything else]"
	@echo "proj        true [anything else]"
	@echo "imsl        vendor imsl [anything else]"
	@echo "openmp      true [anything else]"
	@echo "static      static shared (=dynamic)"

# All dependencies created by perl script make.d.sh
ifeq (False,$(iphonyall))
    $(info Checking dependencies ...)
    -include $(DOBJS) $(FDOBJS) $(CDOBJS)
endif
