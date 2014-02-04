#!/bin/bash
#
# -------------------------------------------------------------------------------------------
# BRIEF
# -------------------------------------------------------------------------------------------
#
# Runs mHM with different settings specified in namelists.
#
#     $ ./check.sh   
#
# This script has to be run successfully without failure statements before
# a new feature is merged back to the trunk of mHM.  
#
# -------------------------------------------------------------------------------------------
# DETAILS
# -------------------------------------------------------------------------------------------
#
# To run all the test cases type 

#      $ cd check/
#      $ ./check.sh

# on your command line. If test case was running properly, i.e. is 
# consistent with the pre-defined output, command line output will be 

#      $ case #[id] o.k. 

# otherwise, if the test scenario failed either during runtime or 
# because of inconsistent output 

#      $ case #[id] failed

# will be written on standard output. Additionally, you will find a file
# mhm.output in each specific test case folder case_[id]/ containing the 
# mHM print-outs. Inspecting this file might help you to find the reason 
# of failure. Only if all test cases pass the make test, a new feature
# of mHM is allowed to be merged back to the trunk which consists of two
# steps:

# (1) Update your branch to current status of the trunk, i.e.
#     include features already apparent in the trunk but not your branch

#      $ cd my-branch/
#      $ svn merge https://svn.ufz.de/svn/mhm/trunk
#      $ svn commit -m "mhm: branches: my-branch: I have the current trunk \
#                       in my branch and solved all conflicts"

# (2) Merge your branch to the trunk

#      $ cd trunk/
#      $ svn merge --reintegrate https://svn.ufz.de/svn/mhm/branches/my-branch
#      $ svn commit -m "mhm: trunk: successfully integrated feature \
#                       <this-is-the-feature>"

# Please be aware that the test cases only include switches specified 
# within a namelist. This means that you should run "./check.sh" with 
# the Makefile setting for

#     (1a) compiler = Gnu
#     (1b) compiler = Nag
#     (1c) compiler = Intel
#     (2a) openmp   = true
#     (2b) openmp   = false

# The test cases tested so far are:

#    case 0 :: DEFAULT

#    case 1 :: RESTART READING AND WRITING
#       restart_flag_states_read   = .TRUE.      (default: .FALSE.)
#       restart_flag_states_write  = .TRUE.      (default: .TRUE. )
#       restart_flag_config_read   = .TRUE.      (default: .FALSE.)
#       restart_flag_config_write  = .TRUE.      (default: .TRUE. )

#    case 2 :: ROUTING DEACTIVATED
#       processCase(8) = 0                       (default: 1)

#    case 3 :: OPTIMIZATION ACTIVATED
#       optimize       = .TRUE.                  (default: .FALSE.)
#       opti_method    = 1                       (default: 3      )
#       nIterations    = 6                       (default: 400    )

# -------------------------------------------------------------------------------------------
#
#
# Copyright 2014 Juliane Mai
#
# License
# This file is part of the UFZ CHS mesoscale hydrologic model mHM.
#
# The UFZ CHS mesoscale hydrologic model mHM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The UFZ CHS mesoscale hydrologic model mHM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# -------------------------------------------------------------------------------------------
#
set -e
#
# Perform a cleanup if script is interupted
trap cleanup 1 2 3 6
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$
#
# ---------------------------------------------------------------------------------------------------------------------
# ./run_mucm.sh -p -c 1.24562587653 iter_1/
#
function usage () {
    printf "${pprog} [directory]\n"
    printf "Checks mHM in for different scenarios specified in namelists.\n"
    printf "\n"
    printf "Input\n"
    printf "    None. \n"
    printf "\n"
    printf "Options\n"
    printf "    None.\n"
    printf "\n"
    printf "Example\n"
    printf "    ${pprog} \n"
}
#
# cleanup at end wnd when interupted
# function cleanup ()
# {
#   # things to remove or clean
# }

mhm_exe='mhm_for_check'

# Compile mHM
cd ..
make PROGNAME=${mhm_exe}

cd check
case_dirs=$(\ls -d case_*)
#case_dirs=$(\ls -d case_0)

for case_dir in ${case_dirs} ; do
    echo ''
    echo '//////////////////////////////////////////////////////////////////////////////////////////////////////////////////'
    echo '                                           '${case_dir}' will be checked ...'
    echo '//////////////////////////////////////////////////////////////////////////////////////////////////////////////////'
    echo ''
    cd ${case_dir}

    # remove old exe and link recently compiled
    if [ -f ${mhm_exe} ] ; then
	\rm ${mhm_exe}
    fi
    ln -s ../../${mhm_exe}

    # run mhm
    echo '   mhm is running ...'
    ./${mhm_exe} > 'mhm_terminal.out'

    # mhm was running till the end?
    found=$(grep 'mHM: Finished!' mhm_terminal.out)
    if [[ $found == 'mHM: Finished!' ]] ; then 
	echo '   mhm finished'
	
        # compare outputs
        # TODO    

        # cleanup 
        # \rm mhm_terminal.out

	echo ''
        echo '   '${case_dir}' o.k.'
        cd ..
    else
	echo '   mhm was aborted'
	echo ''
        echo '   '${case_dir}' failed'
    fi

done

echo ''

# Clean up mhm directory
#make clean
#\rm mhm_for_check

exit 0
