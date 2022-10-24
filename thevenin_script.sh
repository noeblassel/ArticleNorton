#!/bin/bash

Nmin=$1
Nmax=$2
for N in `seq $Nmin 2 $Nmax`
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=16 thevenin_shear.jl 0.8 0.7 1e-3 1.0 0.5 SINUSOIDAL 100.0 100 $N BAOAB 2.5 2.0 0.5 > nohup.out &
done
