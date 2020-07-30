#!/bin/bash
for pdfile in *.pdf ; do
  convert -verbose -density 300 -resize '800' "${pdfile}" "${pdfile%.*}".png
done
