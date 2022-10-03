using Random, LinearAlgebra,Molly

struct NortonSVIntegrator{TF,TG}
    dt::Float64
    η::Float64
    T::Float64
    γ::Float64
    F::TF # forcing profile
    G::TG # response profile
end
NortonSVIntegrator(dt::Float64, η::Float64, T::Float64, γ::Float64, F::Function, G::Function) = NortonSVIntegrator{typeof(F),typeof(G)}(dt, η, T, γ, F, G)

function Molly.simulate!(sys::System, sim::NortonSVIntegrator, n_steps; parallel::Bool=true, rng=Random.GLOBAL_RNG)
    force_hist=Float64[]

    neighbors = find_neighbors(sys, sys.neighbor_finder; parallel=parallel)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt(1 - α^2)

    accels = accelerations(sys, neighbors; parallel=parallel)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot products
    FdotG = dot(F_y, G_y)

    #initialize state on constant response manifold
    λ = (sim.η-dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y
    
    run_loggers!(sys, neighbors, 0; parallel=parallel)
    for step_n = 1:n_steps
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_13 = (sim.η-dot(G_y,v_x))/FdotG#analytic expression for Lagrange multiplier
        v_x .+= λ_13 * F_y 

        #A step
        sys.coords .+= sys.velocities * sim.dt
        sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

        accels .= accelerations(sys, neighbors; parallel=parallel)

        F_y .= sim.F.(q_y)
        G_y .= sim.G.(q_y)

        FdotG = dot(F_y, G_y)

        (isnan(FdotG)) && return force_hist #abort if system NaNs out
        λ_12 = (sim.η-dot(G_y,v_x))/FdotG #correction term to reproject momenta on manifold
        v_x .+=λ_12 * F_y
        #λ_12/Δt is a good approximation to the term ∇R_q(q_t,p_t)⋅p_t/F(q_t)⋅G(q_t) forcing the dynamics to remain on the cotangent bundle
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_23 = (sim.η-dot(G_y,v_x))/FdotG
        v_x .+= λ_23 * F_y 
        
        #O_step
        velocities .= α*sys.velocities + σ * random_velocities(sys, sim.T; rng=rng)  #equilibrium fd solution
        λ_fd=(sim.η - dot(G_y,v_x))/FdotG #analytic expression for Lagrange multiplier
        v_x .+= λ_fd * F_y

        F_ham=(λ_13+λ_23)/sim.dt
        F_ou=sim.γ*sim.η/FdotG
        F_corr= λ_12/sim.dt 

        push!(force_hist, F_ham+F_corr+F_ou)

        neighbors = find_neighbors(sys, sys.neighbor_finder,neighbors,step_n; parallel=parallel)
        run_loggers!(sys, neighbors, step_n; parallel=parallel)
    end
    return force_hist
end

struct NortonSplitting{S,E,K,F,W,TF,TG}
    dt::S
    r::E
    T::K
    γ::F
    splitting::W
    F::TF # forcing profile
    G::TG # response profile
end

function Norton_O_step!(s,α_eff,σ_eff, γ, rng,temperature,r,v_x, F_y, G_y, FdotG)
    s.velocities .= α_eff * s.velocities + σ_eff * random_velocities(s,temperature; rng=rng)
    λ = (r -dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y

    return γ*r / FdotG # no need to keep martingale part of the forcing
end

function Norton_A_step!(s,dt_eff, neighbors,accelerations,n_threads,r,q_y,v_x,F_y,G_y)
    s.coords .+= s.velocities * dt_eff
    s.coords .= Molly.wrap_coords.(s.coords, (s.boundary,))

    accelerations .= accelerations(s,neighbors, n_threads=n_threads)

    F_y .= F.(q_y)
    G_y .= G.(q_y)
    FdotG = dot(F_y,G_y)
    λ = (r - dot(G_y, v_x))/ FdotG #reprojection in p onto the constant response manifold
    v_x .+= λ * F_y

    return (λ , FdotG) # update FdotG to avoid recomputation
end

function Norton_B_step!(s,dt_eff, accelerations, r, v_x,F_y,G_y,FdotG)
    s.velocities .+= accelerations * dt_eff
    λ = (r- dot(v_x, G_y)) /FdotG
    v_x .+= λ * F_y

    return λ
end


function Molly.simulate!(sys::System, sim::NortonSplitting, n_steps; n_threads::Integer=Threads.nthreads(), rng=Random.GLOBAL_RNG)
    α_eff = exp(-sim.γ * sim.dt / count('O', sim.splitting))
    σ_eff = sqrt(sim.T * (1 - α_eff^2)) #noise on velocities, not momenta

    effective_dts=[sim.dt / count(c,sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    sys.coords = wrap_coords.(sys.coords, (sys.boundary,))
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)

    accels = accelerations(sys, neighbors; n_threads=n_threads)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot product
    FdotG = dot(F_y, G_y)

    #initialize state on constant response manifold
    λ = (sim.r-dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y


    λ_hist=Float64[]

    for step_n=1:n_steps
        λ_A=0.0
        λ_B=0.0
        λ_O=0.0

        for (i,op)=enumerate(sim.splitting)
            if op=='A'
                sys.coords .+= sys.velocities * effective_dts[i]
                sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

                accels .= accelerations(sys,neighbors, n_threads=n_threads)

                F_y .= F.(q_y)
                G_y .= G.(q_y)
                FdotG = dot(F_y,G_y)
                λ = (r - dot(G_y, v_x))/ FdotG #reprojection in p onto the constant response manifold
                v_x .+= λ * F_y

                λ_A += λ

            elseif op=='B'
                ( force_computation_steps[i] ) &&  ( accels .= accelerations(sys,neighbors, n_threads=n_threads) )
                sys.velocities .+= accels * effective_dts[i]
                λ = (r- dot(v_x, G_y)) /FdotG
                v_x .+= λ * F_y
                λ_B += λ

            elseif op=='O'
                sys.velocities .= α_eff * sys.velocities + σ_eff * random_velocities(sys,sim.T; rng=rng)
                λ = (r -dot(G_y,v_x))/FdotG
                v_x .+= λ * F_y

                λ_O += sim.γ*r / FdotG # no need to keep martingale part of the forcing
            end
        end

        λ_est= (λ_A + λ_B) / sim.dt + λ_O / count('O',sim.splitting)

        push!(λ_hist, λ_est)

        run_loggers!(sys,neighbors,step_n)

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n ; n_threads=n_threads)
        end
    end
    
    return λ_hist
end