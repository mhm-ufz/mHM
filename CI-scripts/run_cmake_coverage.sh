#!/bin/bash
for path in ../check/case_[0-9][0-9]; do
    echo $path
    out="$path/output_b1"
    mkdir -p $out
    ./mhm $path
done
