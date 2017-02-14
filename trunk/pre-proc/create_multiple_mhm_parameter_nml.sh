#!/bin/bash
#
# Creates multiple mhm_parameter.nml files. The parameter sets will be read from a file and 
# arranged in the namelist format feasible for mHM.
#
# Copyright 2013 Juliane Mai, Matthias Cuntz
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

# You should have received a copy of the GNU Lesser General Public License
# along with The UFZ bash library.  If not, see <http://www.gnu.org/licenses/>.
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
function usage () {
    printf "${pprog} [-hq] [-n nhead] parameter_file\n"
    printf "Creates multiple namelists mhm_parameter.nml-<ID> with parameter sets specified in an <infile>.\n"
    printf "\n"
    printf "Input\n"
    printf "    parameter_file   File with parameter sets.\n"
    printf "\n"
    printf "Options\n"
    printf "    -h               Prints this help screen.\n"
    printf "    -n nhead    	 # of header lines to skip at the beginning of parameter file.\n"
    printf "                	 Default: 0.\n"
    printf "    -p nskip    	 # of parameter sets to skip at the beginning of parameter file.\n"
    printf "                	 Default: 0.\n"
    printf "    -q          	 Proceed quietly.\n"
    printf "    -z          	 Fill ID with zeros.\n"
    printf "\n"
    printf "Example\n"
    printf "    ${pprog} -n 1 -p 10000 params\n"
}
#
# cleanup at end wnd when interupted
function cleanup ()
{
  \rm -f *.${pid}
}
# ---------------------------------------------------------------------------------------------------------------------
# Start
#
# switches
skip=0
pskip=0
quiet=0
zeros=0
while getopts "hn:p:qz" Option ; do
    case ${Option} in
        h) usage 1>&2; exit 0;;
        n) skip="${OPTARG}";;
        p) pskip="${OPTARG}";;
        q) quiet=1;;
        z) zeros=1;;
        *) printf "Error ${pprog}: unimplemented option.\n\n";  usage 1>&2; exit 1;;
    esac 
done
shift $((${OPTIND} - 1))
#
# Check args
NO_ARGS=1
if [[ $# -lt ${NO_ARGS} ]] ; then
    printf "Error ${pprog}: input files missing.\n\n"
    usage 1>&2
    exit 1
fi
infile="$@"
#
# Check input file exist
if [ ! -f ${infile} ] ; then
    printf "Error ${pprog}: Input file not found: %s\n" ${infile}
    exit 1
fi
#
# determine # of namelist files
anz=$(wc -l  ${infile} | awk '{print $1}')
anz=$(( ${anz}-${skip} ))
if [[ ${quiet} -eq 0 ]] ; then
    echo ' '
    echo ${anz} 'namelists will be written.'
    echo ' '
fi

# ----------------------------------------------------------------------------------------------------
# Here comes the CORE:
# ----------------------------------------------------------------------------------------------------
#   1. replace space by semicolon
#   2. fill in semicolons after parameters given by ${list}
printf "%${anz}s" | tr ' ' '\n' >> semifile.${pid}
delim=' '
if [[ ${skip} -gt 0 ]] ; then
    # deletes trailing and leading spaces and skips some lines
    sed -e 's/^ *//g' -e 's/ *$//g' -e "1,${skip}d" ${infile} | tr -s "${delim}" | tr "${delim}" ';' > sinfile.${pid}
else
    sed -e 's/^ *//g' -e 's/ *$//g' ${infile} | tr -s "${delim}" | tr "${delim}" ';' > sinfile.${pid}
fi

# number of parameters and betas
npara=$(tail -1 sinfile.${pid} | tr ';' '\n' | wc -l)
nbeta=$((npara-43))

# Not all lines in the final namelist contain parameter values, 
# e.g. between parameter 38 and 39 there are 4 lines which do not contain a parameter information.
# <list> contains the information how many lines have to be skipped after which parameter,
# e.g. 4 lines between 38 and 39 --> '38 38 38 38'
# Must be sorted reversely.
list='43 43 43 43 43 38 38 38 38 35 35 35 35 35 30 30 30 30 27 27 27 27 26 26 26 26 9 9 9 9 1 1 1 1'
list="${npara} ${npara} ${list}"

for k in ${list} ; do
    cut -f 1-${k} -d ';' sinfile.${pid} > links.${pid}
    cut -f $((k+1))- -d ';' sinfile.${pid} > rechts.${pid}
    paste -d ';' links.${pid} semifile.${pid} rechts.${pid} > tmptmp.${pid}
    mv tmptmp.${pid} sinfile.${pid}
done
# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# Write the namelist in three files that will be put together later.
#   1. The first few lines of the namelist containing comments, etc.
#   2. The left half of the namelist and
#   3. the right half of the namelist
# The final namelist is then:
#  2 and 3 will be pasted together with the parameter values in between.
#  The whole is appended to 1.
# ----------------------------------------------------------------------------------------------------

#
# The first lines of the namelist --> header.${pid}
str='! This file was automatically generated from file containing multiple parameter sets:''\n'
str=${str}'!     '${PWD}'/'${infile}'\n'
str=${str}'! This is set #xxx''\n' # Put placeholder xxx for later replacement
str=${str}'!''\n'
str=${str}'!global_parameters''\n'
str=${str}'!PARAMETER                       lower_bound  upper_bound          value   FLAG  SCALING''\n'
str=${str}'!interception''\n'
str=${str}'&interception1''\n'
printf "${str}" > header.${pid}

#
# The values will be inserted later between the first few and the last few columns.
# The first columns with variable names and ranges of parameter                   The last columns with Flag and scaling law
# -->  front.${pid}                                                               -->  back.${pid}
# Parameter 1
str='canopyInterceptionFactor           =  0.1500,      0.4000,     ''\n'       ; str1=',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! snow''\n'                                                          ; str1=${str1}'\n'                      
str=${str}'&snow1''\n'                                                          ; str1=${str1}'\n'                     
# Parameter 2
str=${str}'snowTreshholdTemperature           = -2.0000,      2.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 3
str=${str}'degreeDayFactor_forest             =  0.0001,      4.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 4
str=${str}'degreeDayFactor_impervious         =  0.0000,      1.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 5
str=${str}'degreeDayFactor_pervious           =  0.0000,      2.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 6
str=${str}'increaseDegreeDayFactorByPrecip    =  0.1000,      0.9000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 7
str=${str}'maxDegreeDayFactor_forest          =  0.0000,      8.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 8
str=${str}'maxDegreeDayFactor_impervious      =  0.0000,      8.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 9
str=${str}'maxDegreeDayFactor_pervious        =  0.0000,      8.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! soilmoisture''\n'                                                  ; str1=${str1}'\n'                      
str=${str}'&soilmoisture1''\n'                                                  ; str1=${str1}'\n'       
# Parameter 10
str=${str}'orgMatterContent_forest            =  0.0000,      20.000,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 11
str=${str}'orgMatterContent_impervious        =  0.0000,      1.0000,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 12
str=${str}'orgMatterContent_pervious          =  0.0000,      4.0000,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 13
str=${str}'PTF_lower66_5_constant             =  0.6462,      0.9506,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 14
str=${str}'PTF_lower66_5_clay                 =  0.0001,      0.0029,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 15
str=${str}'PTF_lower66_5_Db                   = -0.3727,     -0.1871,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 16
str=${str}'PTF_higher66_5_constant            =  0.5358,      1.1232,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 17
str=${str}'PTF_higher66_5_clay                = -0.0055,      0.0049,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 18
str=${str}'PTF_higher66_5_Db                  = -0.5513,     -0.0913,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 19
str=${str}'PTF_Ks_constant                    = -1.2000,     -0.2850,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 20
str=${str}'PTF_Ks_sand                        =  0.0060,      0.0260,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 21
str=${str}'PTF_Ks_clay                        =  0.0030,      0.0130,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 22
str=${str}'PTF_Ks_curveSlope                  =  1.0000,     150.000,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 23
str=${str}'rootFractionCoefficient_forest     =  0.9000,      0.9990,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 24
str=${str}'rootFractionCoefficient_impervious =  0.9000,      0.9500,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 25
str=${str}'rootFractionCoefficient_pervious   =  0.0010,      0.0900,    ''\n'  ; str1=${str1}',     1,       1''\n'
# Parameter 26
str=${str}'infiltrationShapeFactor            =  1.0000,      4.0000,    ''\n'  ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! directSealedAreaRunoff''\n'                                        ; str1=${str1}'\n'                      
str=${str}'&directRunoff1''\n'                                                  ; str1=${str1}'\n'                      
# Parameter 27
str=${str}'imperviousStorageCapacity          =  0.0000,      5.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! meteoCorrectionFactor''\n'                                         ; str1=${str1}'\n'                      
str=${str}'&actualET1''\n'                                                      ; str1=${str1}'\n'                      
# Parameter 28
str=${str}'minCorrectionFactorPET             =  0.7000,      1.3000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 29
str=${str}'maxCorrectionFactorPET             =  0.0000,      0.2000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 30
str=${str}'aspectTresholdPET                  =  160.00,      200.00,     ''\n' ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! interflow''\n'                                                     ; str1=${str1}'\n'                      
str=${str}'&interflow1''\n'                                                     ; str1=${str1}'\n'                      
# Parameter 31
str=${str}'interflowStorageCapacityFactor     =  75.000,      200.00,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 32
str=${str}'interflowRecession_slope           =  0.0000,      10.000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 33
str=${str}'fastInterflowRecession_forest      =  1.0000,      3.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 34
str=${str}'slowInterflowRecession_Ks          =  1.0000,      30.000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 35
str=${str}'exponentSlowInterflow              =  0.0500,      0.3000,     ''\n' ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! percolation''\n'                                                   ; str1=${str1}'\n'                      
str=${str}'&percolation1''\n'                                                   ; str1=${str1}'\n'                      
# Parameter 36
str=${str}'rechargeCoefficient                =  0.0000,      50.000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 37
str=${str}'rechargeFactor_karstic             = -5.0000,      5.0000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 38
str=${str}'gain_loss_GWreservoir_karstic      =  1.0000,      1.0000,     ''\n' ; str1=${str1}',     0,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'                      
str=${str}'''\n'                                                                ; str1=${str1}'\n'                      
str=${str}'! routing''\n'                                                       ; str1=${str1}'\n'                      
str=${str}'&routing1''\n'                                                       ; str1=${str1}'\n'                      
# Parameter 39
str=${str}'muskingumTravelTime_constant       =  0.3100,      0.3500,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 40
str=${str}'muskingumTravelTime_riverLength    =  0.0700,      0.0800,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 41
str=${str}'muskingumTravelTime_riverSlope     =  1.9500,      2.1000,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 42
str=${str}'muskingumTravelTime_impervious     =  0.0900,      0.1100,     ''\n' ; str1=${str1}',     1,       1''\n'
# Parameter 43
str=${str}'muskingumAttenuation_riverSlope    =  0.0100,      0.5000,     ''\n' ; str1=${str1}',     1,       1''\n'
str=${str}'/''\n'                                                               ; str1=${str1}'\n'
str=${str}'''\n'                                                                ; str1=${str1}'\n'
str=${str}'! geological parameters (ordering according to file "geology_classdefinition.txt")''\n' ; str1=${str1}'\n'                      
str=${str}'! this parameters are NOT REGIONALIZED yet, i.e. these are <beta> and not <gamma>''\n'  ; str1=${str1}'\n'                      
str=${str}'&geoparameter''\n'                                                                      ; str1=${str1}'\n'  
# Parameter 44 ... npara = 43+nbeta
for ((i=1 ; i<=${nbeta} ; i++)) ; do      
    str=${str}'GeoParam('${i}',:)                      =  1.0000,      1000.00,     ''\n' ; str1=${str1}',     1,       1''\n';
done
str=${str}'/''\n'                                                                ; str1=${str1}'\n'
str=${str}'''\n'                                                                 ; str1=${str1}'\n'
printf "${str}" > front.${pid}                                                   ; printf "${str1}" > back.${pid}

# ----------------------------------------------------------------------------------------------------
# The result parameter file of the CORE is tinfile...
# ----------------------------------------------------------------------------------------------------

# Now put everything together...
i=0
while read line ; do
    i=$((i+1))
    if [[ ${i} -le ${pskip} ]] ; then continue ; fi
    # New file name
    if [[ ${zeros} -eq 1 ]] ; then
	newfile=$(printf "mhm_parameter.nml-%0${nanz}i" ${i})
    else
	newfile="mhm_parameter.nml-${i}"
    fi
    # If namelist already exists
    if [[ ${quiet} -eq 0 ]] ; then
	if [ -f ${newfile} ] ; then
            printf "Deleting and newly generating: ${newfile}\n"
            rm ${newfile}
	else
            printf "Generating: ${newfile}\n"
	fi
    else
	if [ -f ${newfile} ] ; then rm ${newfile} ; fi
    fi

    # Replace id in header
    sed -e "s/xxx/${i}/" header.${pid} > header1.${pid}
    # (1) Cutting the ith column from prepared parameter set file "tinfile" and 
    # (2) paste it between front.${pid} and back.${pid} and
    # (3) append that to the header.${pid} and
    # (4) store everything in <newfile>
    echo ${line} | tr ';' '\n' | paste front.${pid} - back.${pid} | cat header1.${pid} - > ${newfile}
    \rm header1.${pid}
done < sinfile.${pid}

# Perform cleanup, i.e. all temporary files *.${pid} are deleted
cleanup

exit 0
