#!/bin/bash

for N in `seq 15`
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 thevenin_shear.jl 0.8 0.7 1e-3 1.0 0.3 SINUSOIDAL 100.0 1000 $N OBABO 2.5 > nohup.out &
done