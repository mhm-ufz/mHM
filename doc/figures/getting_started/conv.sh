#!/bin/bash
for pdfile in *.pdf ; do
  convert -verbose "${pdfile}" "${pdfile%.*}".png
done
