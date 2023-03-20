for eta in 10.0 20.0 30.0 40.0 50.0 60.0 70.0 
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 thevenin_cd.jl 2.5 1.0 1e-3 1.0 $eta 10.0 1000 400 BAOAB > nohup.out &
done