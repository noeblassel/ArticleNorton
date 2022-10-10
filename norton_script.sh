#!/bin/bash
Nmin=$1
Nmax=$2
for N in `seq $Nmin $Nmax`
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 norton_shear.jl 0.8 0.7 1e-3 1.0 0.1 SINUSOIDAL 100.0 1000 $N 2.5 > nohup.out &
done