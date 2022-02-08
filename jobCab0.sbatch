#!/bin/bash
source ~/.bashrc

We="20.0"
Ohd="1e-5"
Ohs="1e-5"
Bo="0.25"
Lo="4.0"
Level="10"
tmax="10.0"

qcc -fopenmp -Wall -O2 bounce.c -o bounce -lm
export OMP_NUM_THREADS=64
./bounce $Level $tmax $We $Ohd $Ohs $Bo $Lo
