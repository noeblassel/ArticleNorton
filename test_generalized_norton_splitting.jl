using Molly, Plots


(!isdir("gal_norton_plots")) && mkdir("gal_norton_plots")

include("utils.jl")
include("animate.jl")
include("norton_integrators.jl")

γ_r = 1.0
r_bar = 10.0

φ(r,dt)= r_bar
ξ(dt) = 0.0

F_const(y) = (y < L / 2) ? 1.0 : -1.0
F_sin(y) = sin(2π * y / L)
G_imag(y) = sin(2π * y / L) / N

T=0.8
γ=1.0
ρ=0.7
dt=2e-3
r_c =2.5
Npd=10

N=Npd^3
L=(N/ρ)^(1//3)
box_size=CubicBoundary(L,L,L)

max_speed=10.0*sqrt(T)
n_steps_neighbors=floor(Int64,0.2*r_c/(dt*max_speed))

splitting="BAOAB"
simulator=GeneralizedNortonSplitting(dt,T,γ,splitting, F_sin, G_imag, φ, ξ)

nf = (3.6r_c < L) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff=1.2r_c)

coords = place_atoms_on_3D_lattice(Npd, box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq= 2000

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0,loggers=(coords=GeneralObservableLogger((sys,args...;kwargs...)->deepcopy(sys.coords),typeof(coords),1),))

λ_hist,λ_mart_hist,r_hist=simulate!(sys,simulator,n_steps_eq)
T_range=dt*(1:n_steps_eq)
println("done simulating")
plot(T_range,λ_hist,xlabel="t",ylabel="λ",label="finite variation part")
savefig("gal_norton_plots/lambda_bar_sin_eq.pdf")

plot(T_range,λ_mart_hist,xlabel="t",ylabel="λ",label="martingale part")
savefig("gal_norton_plots/lambda_tilde_sin_eq.pdf")

plot(T_range,r_hist,xlabel="t",ylabel="r",label="response",legend=:bottom)
savefig("gal_norton_plots/response_sin_eq.pdf")

animate_system(sys,"gal_norton_plots/anim_sin_eq.mp4",F_sin)