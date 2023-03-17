#!/bin/bash

Nmin=$1
Nmax=$2
method=$3

for N in `seq $Nmin 2 $Nmax`
do

if [ $method == "norton" ]

then
    nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=16 norton_cd.jl 1.25 0.6 5e-3 1.0 1.21 10.0 1000 $N BAOAB > nohup.out &
else
    nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=16 thevenin_cd.jl 1.25 0.6 5e-3 1.0 10.0  10.0 1000 $N BAOAB > nohup.out &
fi
done