#!/bin/bash
version="'"$(cat version.txt)"'"
version_date="'"$(cat version_date.txt)"'"
if [[ -n $FC && -n $(command -v $FC) ]]; then
  if [[ $FC == *nagfor ]]; then
    $FC -fpp -F -DMHMVERSION="${version}" -DMHMDATE="${version_date}" -o tempfile-doxygen.fpp $1 &>/dev/null
    cat tempfile-doxygen.fpp
    rm -f tempfile-doxygen.fpp
  else
    $FC -cpp -E -DMHMVERSION="${version}" -DMHMDATE="${version_date}" $1
  fi
else
  if [[ -n $(command -v gfortran) ]]; then
    gfortran -cpp -E -DMHMVERSION="${version}" -DMHMDATE="${version_date}" $1
  fi
fi
