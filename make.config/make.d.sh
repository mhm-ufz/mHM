#!/bin/bash

#
# Produces makefile dependencies of Fortran files
#
# Copyright 2013-2015 Matthias Cuntz - mc (at) macu.de
#
# License
# This file is part of the makefile library.
#
# The UFZ bash library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The UFZ bash library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with The UFZ bash library. If not, see <http://www.gnu.org/licenses/>.

set -e
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

# Perform a cleanup if script is interupted
trap cleanup 1 2 3 6

# --------------------------------------------------------------------------------------------------
# functions
#
function usage () {
    printf "${pprog} [-h] FortranFile FortranFileName Source2ObjectPath AllSrcFiles\n"
    printf "\n"
    printf "Produces makefile dependency files with file ending .d of Fortran files.\n"
    printf "\n"
    printf "Look for USE statements in FortranFile "
    printf "and search in AllSrcFiles for the files that provide the modules.\n"
    printf "Write a dependency file .d in Source2ObjectPath for FortranFileName and not FortranFile. "
    printf "This allows to check preprocessed Fortran files but write the dependency for the original files. "
    printf "Use the filename twice if no preprocessed file.\n"
    printf "\n"
    printf "Input\n"
    printf "    FortranFile         Fortran file for which dependency file will be produced.\n"
    printf "    FortranFileName     Fortran file name used in dependency file.\n"
    printf "    Source2ObjectPath   Script assumes that object path is $(dirname ${FortranFile})/${Source2ObjectPath}).\n"
    printf "    AllSrcFiles         All source files of project that will be scanned for dependencies.\n"
    printf "\n"
    printf "Options\n"
    printf "    -h          Prints this help screen.\n"
    printf "\n"
    printf "Examples\n"
    printf "    ${pprog} src/main.f90     src/main.f90 .gnu.release src/*.f90\n"
    printf "or\n"
    printf "    ${pprog} src/main.f90.pre src/main.f90 .gnu.release src/*.f90\n"
}
#
# cleanup at end wnd when interupted
function cleanup ()
{
  \rm -f *.${pid}
}

# --------------------------------------------------------------------------------------------------
# Arguments
#
# Switches
while getopts "h" Option ; do
    case ${Option} in
	h) usage; exit;;
	*) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

# Check that enough arguments
if [[ $# -lt 3 ]] ; then
    printf "Error ${pprog}: not enough input arguments.\n\n" 1>&2
    usage 1>&2
    exit 1
fi

# infile and objectpath
thisfile=$1
thisfilename=$2
src2obj=$3
shift 3
# all source files
srcfiles=$@

# --------------------------------------------------------------------------------------------------
# Dependencies
#
# All module names and filenames into a dictionary
# Same dictionary for all input directories
firstdir=$(dirname $(echo ${srcfiles} | tr ' ' '\n' | sort | uniq | sed '/^$/d' | tr '\n' '\t' | cut -f 1))
# dictionary in first input directory
dict="${firstdir}/${src2obj}/${pprog}.dict"
if [[ ! -f ${dict} ]] ; then # new dict only if it does not exist in directory yet
    if [[ ! -d $(dirname ${dict}) ]] ; then mkdir -p $(dirname ${dict}) ; fi
    for i in ${srcfiles} ; do # all files in all input dirs
	# 0. cariage returns to newlines,
	# 1. all blanks to one space, 2. rm f90 comments, 3. rm f77 comments,
	# 4. rm leading blank, 5. rm trailing blank, 6. rm blank lines,
	# 7. lowercase, 8. squeeze blanks (redundant),
	# 9. lines with word 'module', 10. lines with only 'module name', i.e. not module procedure,
	# 11. rm word 'module'
	ismod=$(echo "${i}:$(tr '\r' '\n' < ${i} | sed -e 's/[[:blank:]]\{1,\}/ /g' -e 's/ *\!.*//' -e '/^[Cc]/d' -e 's/^ //' -e 's/ $//' -e 's/^$//' | tr [A-Z] [a-z] | sed -n -e '/^module /p' | sed -n -e '/^module [[:alnum:]_]\{1,\}$/p' | sed -e 's/^module //')")
	if [[ "${ismod}z" != "${i}:z" ]] ; then echo ${ismod} >> ${dict} ; fi
    done
fi

# Modules used in the input file
molist=$(sed -e 's/\!.*//' -e '/^[Cc]/d' ${thisfile} | tr [A-Z] [a-z] | tr -s ' ' | grep -E '^[[:blank:]]*use[[:blank:]]+' | sed 's/,.*//' | sed 's/.*use //' | sort | uniq)
is=$(echo ${molist} | tr ' ' '|')

# Query dictionary for filenames of modules used in input file
# Remove own file name for circular dependencies if more than one module in input file
if [[ "${is}" != "" ]] ; then
    olist=$(cut -f 1 -d ':' ${dict} | sed -n $(echo $(grep -nEw "${is}" ${dict} | cut -f 1 -d ':') | sed -e 's/\([0-9]*\)/-e \1p/g') | tr '\n' ' ' | sed "s|${thisfilename}||")
fi

# Write output .d file
s2ofile="$(dirname ${thisfilename})/${src2obj}/$(basename ${thisfilename})"
tmpfile=${s2ofile}.${pid}
printf "${s2ofile/\.[fF]*/.d} : ${thisfilename}\n" > ${tmpfile}
printf "${s2ofile/\.[fF]*/.o} : ${s2ofile/\.[fF]*/.d}" >> ${tmpfile}
for i in ${olist} ; do
    is2ofile="$(dirname ${i})/${src2obj}/$(basename ${i})"
    printf " ${is2ofile/\.[fF]*/.o}" >> ${tmpfile}
done
printf "\n" >> ${tmpfile}
outfile=${s2ofile/\.[fF]*/.d}
mv ${tmpfile} ${outfile}

# replace .d file if changed
# outfile=${s2ofile/\.[fF]*/.d}
# if [[ -f ${outfile} ]] ; then
#     set +e
#     tt=$(sdiff -s ${tmpfile} ${outfile})
#     set -e
#     if [[ "${tt}" != "" ]] ; then
# 	mv ${tmpfile} ${outfile}
#     else
# 	rm ${tmpfile}
#     fi
# else
#     mv ${tmpfile} ${outfile}
# fi

# Perform cleanup, i.e. all temporary files *.${pid} are deleted
cleanup

exit 0
