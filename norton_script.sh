#!/bin/bash
for N in `seq 6 15`
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 test_norton_splitting.jl 0.8 0.7 1e-3 1.0 0.1 SINUSOIDAL 100.0 10000 $N 2.5 > nohup.out &
done