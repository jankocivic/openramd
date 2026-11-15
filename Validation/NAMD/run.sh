#!/bin/bash

for i in $(seq 1 100); do
    ~/Software/namd_3.0.2/namd3 +p16 +setcpuaffinity +devices 0 namd.input > run.log
    tail run.log > run${i}.log
    tail -n 1 ramd_run1.log > ramd${i}.log
    rm run*.log
    rm ramd_run1.log
    find . -maxdepth 1 -type f ! -name namd.input ! -name 'ramd-5*.tcl' ! -name 'run.sh*' -delete
done
