using Molly, LinearAlgebra

include("utils.jl")
println("Usage: T ρ dt γ η t_equilibration n_iter_sim N scheme")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
η=parse(Float64,ARGS[5])

t_eq=parse(Float64,ARGS[6])
n_iter_sim=parse(Int64,ARGS[7])

N=parse(Int64,ARGS[8])
splitting=ARGS[9]


L=cbrt(N/ρ)
box_size = CubicBoundary(L,L,L)

color_drift = ones(N)
color_drift[1:2:end] .= -1
color_drift /= sqrt(N)

struct ColorDriftForcing
    η::Float64
end

function Molly.forces(inter::ColorDriftForcing,s::System,neighbors=nothing)
    f=zero(s.velocities)
    f_x=view(reinterpret(reshape,Float64,f),1,:)
    f_x .= color_drift
    return forcing.η*f
end

function colordrift_response(s::System,args...;kwargs...)
    p_x=view(reinterpret(reshape,Float64,s.velocities),1,:)
    return dot(p_x,color_drift)
end

forcing=ColorDriftForcing(η)

atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms(N,box_size;min_dist=0.8)
velocities=[velocity(1.0,T,1.0) for i=1:N]
inter=LennardJones(force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))

sim=LangevinSplitting(dt=dt,friction=γ,temperature=T,splitting=splitting;remove_CM_motion=false)
sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,force_units=NoUnits,energy_units=NoUnits,k=1.0)

simulate!(sys,sim,n_steps_eq)
sys=System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,force_units=NoUnits,energy_units=NoUnits,k=1.0,loggers=(mobility=GeneralObservableLogger(colordrift_response,Float64,1),))

for i=1:n_iter_sim
    simulate!(sys,sim,n_steps_eq)
    f=open("thermo_results/nemd_response_colordrift_$(η)_$(N).out","a")
    write(f,values(sys.loggers.mobility))
    close(f)
    empty!(sys.loggers.mobility.history)
end