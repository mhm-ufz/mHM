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
# To run all the test cases on (s)ystem "eve" with (c)ompiler "nag53" in (r)elease mode "debug"
# type 

#      $ cd check/
#      $ ./check.sh -s eve -c nag53 -r debug

# on your command line. If test case was running properly, i.e. is 
# consistent with the pre-defined output, command line output will be 

#      $     ###########################
#      $     #     case_[id]  o.k.     #
#      $     ###########################
#      $     using
#      $        System:   eve
#      $        Compiler: nag53
#      $        Release:  debug

# otherwise, if the test scenario failed either during runtime or 
# because of inconsistent output 

#      $     ###########################
#      $     #     case_[id]  failed   #
#      $     ###########################
#      $     using
#      $        System:   eve
#      $        Compiler: nag53
#      $        Release:  debug

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
# different compilers and in release/debug. These settings can be made 
# by using the optional arguments of check.sh, i.e. "-c", "-r", "-s".

# The test cases specified so far are:

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
# eve 
#     :: gnu45   :: release  --> ok
#     :: gnu46   :: release  --> ok
#     :: intel11 :: release  --> ok
#     :: intel12 :: release  --> ok
#     :: nag53   :: release  --> 

#     :: gnu45   :: debug    --> ok
#     :: gnu46   :: debug    --> ok
#     :: intel11 :: debug    --> ok
#     :: intel12 :: debug    --> ok
#     :: nag53   :: debug    --> ok

set -e
set -o pipefail
#
# Perform a cleanup if script is interupted
# trap cleanup 1 2 3 6  # output and error messages will be lost
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$
#
# ---------------------------------------------------------------------------------------------------------------------
#
function usage () {
    printf "${pprog} [directory]\n"
    printf "Checks mHM for different scenarios specified in namelists (case_* folders). \n"
    printf "The (s)ystem, (c)ompiler, and (r)elease settings used by the Makefile can be specified optionally.\n"
    printf "\n"
    printf "Input\n"
    printf "    None.\n"
    printf "\n"
    printf "Options\n"
    printf "    -h     Prints this help screen.\n"
    printf "    -c     Compilers used to generate a mHM executable, e.g. 'gnu' \n"
    printf "           (DEFAULT: all found on specified system)\n"
    printf "    -r     Release or debug mode, e.g. 'release'\n"
    printf "           (DEFAULT: debug)\n"
    printf "    -s     System on which mHM will be compiled, e.g. 'eve' "\n
    printf "           (DEFAULT: eve). \n"
    printf "\n"
    printf "Example\n"
    printf "    ${pprog} -c gnu45 -s eve -r debug\n"
}

# cleanup at end wnd when interupted
function cleanup ()
{
    cd ${isdir}
    if [ -f 'mhm_for_check' ] ;    then \rm 'mhm_for_check' ;    fi

    case_dirs=$(\ls -d case_*)
    for case_dir in ${case_dirs} ; do
	cd ${case_dir}
	if [ -f 'cdo.error' ] ;        then \rm 'cdo.error' ;        fi
	if [ -f 'cdo.out' ] ;          then \rm 'cdo.out' ;          fi
	if [ -f 'cdo.terminal' ] ;     then \rm 'cdo.terminal' ;     fi
        if [ -f 'dds_results.out' ] ;  then \rm 'dds_results.out' ;  fi
        if [ -f 'mhm_terminal.out' ] ; then \rm 'mhm_terminal.out' ; fi
	if [ -f 'mhm_for_check' ] ;    then \rm 'mhm_for_check' ;    fi
	\rm -r output_b1
	mkdir output_b1
	\rm -r output_b2
	mkdir output_b2
	cd ..
    done

    cd ..
    make clean
}

# global variables
mhm_exe='mhm_for_check'

# switches
compilers=''
system='eve'
release='debug'
while getopts "c:hr:s:" Option ; do
    case ${Option} in
        h) usage 1>&2; exit 0;;
	c) compilers="${OPTARG}";;
	r) release="${OPTARG}";;
	s) system="${OPTARG}";;
        *) printf "Error ${pprog}: unimplemented option.\n\n";  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

if [[ ${compilers} == '' ]] ; then
    compilers=$(\ls ../make.config/${system}.* | grep -v 'alias' | grep -v 'gnu41' | grep -v 'gnu44' | cut -d '/' -f 3 | cut -d '.' -f 2)
fi

for compiler in ${compilers} ; do

    # Compile mHM
    cd ..
    # echo ''
    # echo '//////////////////////////////////////////////////////////////////////////////////////////////////////////////////'
    # echo 'Setup compilation: '
    # echo 'make system='${system}' compiler='${compiler}' release='${release}' PROGNAME='${mhm_exe}
    # echo '//////////////////////////////////////////////////////////////////////////////////////////////////////////////////'
    # echo ''
    make system=${system} compiler=${compiler} release=${release} PROGNAME=${mhm_exe}

    cd check
    case_dirs=$(\ls -d case_*)
   # case_dirs=$(\ls -d case_0)

    nRecordsDiffer_allCases=0
    for case_dir in ${case_dirs} ; do
        echo ''
        echo '##################################################################################################################'
        echo '#                                          '${case_dir}' will be checked ...                                            #'
        echo '##################################################################################################################'
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
        if [[ ${found} == 'mHM: Finished!' ]] ; then 
            echo '   mhm finished'
            
            # load cdo
            cdo -h diffn 1> cdo.out 2> cdo.error || true
            if [[ $(more cdo.error) != '' ]] ; then 
                echo '   try to load module cdo ...'
                source /etc/profile.d/000-modules.sh  || true
                module load cdo  || true
                cdo -h diffn 1> cdo.out 2> cdo.error || true
                if [[ $(more cdo.error) != '' ]] ; then
                    echo '   cdo not installed, but necessary for checking mhm outputs'
                    echo '   Abort'
                    exit 1
                else
                    echo '   cdo is loaded'
                fi
            else
                echo '   cdo is installed'
            fi
    	    echo ''

            # compare outputs in output_b1/ with reference ones in output_save/
    	    nRecordsDiffer=0
    	    cd output_save
    	    nc_files=$(\ls *.nc)
    	    cd ..
    	    for nc_file in ${nc_files} ; do
    	    # cdo -b 16 :: precision=double=16bit
		
    		echo '   *  cdo diff output_save/'${nc_file}' output_b1/'${nc_file}
    		cdo -s diffn output_save/${nc_file} output_b1/${nc_file} > cdo.terminal
    		iRecordsDiffer=$( grep -i 'differ' cdo.terminal | tail -1 | cut -f 1 -d 'o')
    		nRecordsDiffer=$((nRecordsDiffer+iRecordsDiffer))
		nRecordsDiffer_allCases=$((nRecordsDiffer_allCases+iRecordsDiffer))
    		echo '      Number of records that differ more than 0.001: '${iRecordsDiffer}
    		if ((iRecordsDiffer > 0)) ; then
    		    echo '            -------------------'
    		    echo '            Details '
    		    echo '            -------------------'
    		    cdo diffn output_save/${nc_file} output_b1/${nc_file}
    		    echo '            -------------------'
    		fi
    		
    		echo ''
    	    done
    	    echo '   Total Number of records that differ more than 0.001: '${nRecordsDiffer}
    	    
    	    if (( ${nRecordsDiffer} == 0 )) ; then
    		echo ''
		echo '   ###########################'
		echo '   #     '${case_dir}'     o.k.     #'
    		echo '   ###########################'
		echo '   using'
		echo '         System:   '${system} 
		echo '         Compiler: '${compiler}
		echo '         Release:  '${release}
    		echo ''
    	    else
    		echo ''
		echo '   ###########################'
		echo '   #     '${case_dir}'    failed    #'
    		echo '   ###########################'
		echo '   using'
		echo '         System:   '${system} 
		echo '         Compiler: '${compiler}
		echo '         Release:  '${release}
    		echo ''
    	    fi
            cd ..
        else
            echo '   mhm was aborted'
            echo ''
            echo '   ###########################'
            echo '   #     '${case_dir}'    failed    #'
    	    echo '   ###########################'
	    echo '   using'
	    echo '         System:   '${system} 
	    echo '         Compiler: '${compiler}
	    echo '         Release:  '${release}
    	    echo ''
        fi
	
    done # case loop

    if ((nRecordsDiffer_allCases > 0)) ; then
	echo ''
	echo '##################################################################################################################'
        echo '#     At least one case failed in this comiler setup                                                             #'
	echo '#     check.sh is aborted without cleaning up the output in the case directories,                                #'
	echo '#     such that it can be used to identify the bug                                                               #'
    	echo '##################################################################################################################'
	exit 1
    fi
	
done # compiler loop

# Clean up mhm directory
cleanup
echo ''
echo '##################################################################################################################'
echo '#     All cases with the specified compilers were o.k.                                                          #'
echo '#     Output is cleaned up.                                                                                      #'
echo '##################################################################################################################'
echo ''

exit 0
