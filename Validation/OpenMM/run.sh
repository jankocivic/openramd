#!/bin/bash

for i in $(seq 1 100); do
    python3 run.py > run${i}.log 2>&1
    tail -n 1 ramd.log > ramd${i}.log
    rm ramd.log
    rm log.txt
    rm run${i}.log
    find . -maxdepth 1 -type f ! -name '*.log' ! -name '*.txt' ! -name '*.py' ! -name 'run.sh*' -delete
done
