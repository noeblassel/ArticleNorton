using Molly

include("utils.jl")
include("norton_integrators.jl")

println("Usage: T ρ dt γ r forcing_type=SINUSOIDAL|LINEAR|CONSTANT t_equilibration n_iter_sim Nx splitting cutoff_radius y_ratio z_ratio")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
r=parse(Float64,ARGS[5])

forcing_type=ARGS[6]
t_eq=parse(Float64,ARGS[7])
n_iter_sim=parse(Int64,ARGS[8])

Nx=parse(Int64,ARGS[9])
splitting=ARGS[10]
r_c=parse(Float64,ARGS[11])
y_ratio = parse(Float64,ARGS[12])
z_ratio = parse(Float64,ARGS[13])

Ny=round(Int64,Nx*y_ratio)
Nz=round(Int64,Nx*z_ratio)

N=Nx*Ny*Nz

Lx=Nx*cbrt(ρ)
Ly=Lx*y_ratio
Lz=Lx*z_ratio

box_size=CubicBoundary(Lx,Ly,Lz)

max_speed=10.0*sqrt(T)
n_steps_neighbors=floor(Int64,0.2*r_c/(dt*max_speed))


F_sin(y) = sin(2π * y / Ly)
F_const(y) = (y < Ly / 2) ? 1.0 : -1.0
F_lin(y) = (y < Ly / 2) ? 4 * (y - Ly / 4) / L : 4 * (3Ly / 4 - y) / Ly

G_imag(y) = sin(2π * y / Ly) / N
G_real(y) = cos(2π * y/ Ly) / N

F=F_sin
G=G_imag

if forcing_type=="LINEAR"
    F=F_lin
    G=G_real
elseif forcing_type=="CONSTANT"
    F=F_const
end

simulator=NortonSplitting(dt,r,T,γ,splitting,F,G)

nf = (3.6r_c < min(Lx,Ly,Lz)) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff=1.2r_c)
coords=place_atoms_on_3D_lattice(Nx,Ny,Nz,box_size)
atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=floor(Int64,t_eq/dt)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)
_=simulate!(sys,simulator,n_steps_eq)

loggers=(temp=TemperatureLogger(Float64,1),)
sys= System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0, loggers=loggers)

for i=1:n_iter_sim
    force = simulate!(sys,simulator,n_steps_eq)
    f=open("thermo_results/norton_forcing_$(forcing_type)_$(r)_$(N).out","a")
    write(f,force)
    close(f)

    temps=values(sys.loggers.temp)
    f=open("thermo_results/norton_temp_$(forcing_type)_$(r)_$(N).out","a")
    write(f,temps)
    close(f)
    #println(sum(values(sys.loggers.temp))/length(values(sys.loggers.temp)))
    empty!(sys.loggers.temp.history)
end