for eta in 0.1 0.2 0.5 1.0 2.0 5.0 10.0 20.0
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 thevenin_cd.jl 2.5 1.0 1e-3 1.0 $eta 10.0 1000 400 BAOAB > nohup.out &
done