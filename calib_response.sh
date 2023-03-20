for eta in 100.0 200.0 300.0 400.0 500.0 600.0 700.0 
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 thevenin_cd.jl 2.5 1.0 1e-3 1.0 $eta 10.0 1000 400 BAOAB > nohup.out &
done