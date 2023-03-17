using Molly, LinearAlgebra

include("utils.jl")
include("norton_integrators.jl")

println("Usage: T ρ dt γ r t_equilibration n_iter_sim N scheme")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
r=parse(Float64,ARGS[5])

t_eq=parse(Float64,ARGS[6])
n_iter_sim=parse(Int64,ARGS[7])

N=parse(Int64,ARGS[8])
splitting=ARGS[9]


L=cbrt(N/ρ)
box_size = CubicBoundary(L,L,L)

simulator=NortonSplittingColorDrift(dt,r,T,γ,N,splitting)

coords=place_atoms(N,box_size;min_dist=0.8)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

inter=LennardJones(force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=floor(Int64,t_eq/dt)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),boundary=box_size,force_units=NoUnits,energy_units=NoUnits,k=1.0)
_=simulate!(sys,simulator,n_steps_eq)

sys= System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),boundary=box_size,force_units=NoUnits,energy_units=NoUnits,k=1.0)

for i=1:n_iter_sim
    println("iteration $i")
    force = simulate!(sys,simulator,n_steps_eq)
    f=open("thermo_results/norton_forcing_colordrift_$(r)_$(N).out","a")
    write(f,force)
    close(f)
end