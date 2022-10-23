using LinearAlgebra,Base.Threads

N=1000
β=5.0
η=15.0


function F!(q,F_vect;num_threads=Threads.nthreads())
    N=first(size(q))
    chunk_size = Int(ceil(N / num_threads))
    ix_ranges = [i:min(i + chunk_size - 1, N) for i in 1:chunk_size:N]

    @threads for ixs in ix_ranges
        for i in ixs
            x,y=q[i,:]
            r=sqrt(x^2+y^2)
            F_vect[i,:]=[y/r,-x/r]  
        end
    end

end

function simulate_euler_maruyama!(q,Δt,n_iterations,β,η; record_hist::Bool=false, num_threads=Threads.nthreads())

    N=first(size(q))
    shape=size(q)
    σ=sqrt(2Δt/β)

    F_vect=zero(q)

    if record_hist
        #q_record=typeof(q)[]
        r_record=Float64[]
    end

    for i=1:n_iterations
        F!(q,F_vect;num_threads=num_threads)

        if record_hist
             #push!(q_record,copy(q))
             push!(r_record,dot(q,q)/2N)
        end

        q.+= -Δt*q + σ*randn(shape) +η*F_vect*Δt # Euler-Maruyama for overdamped Langevin
    end

    if record_hist
        return r_record .- 1/β
    else
        return
    end

end


#= 

n_eq_steps=5000
n_sim_steps=5000000

q_history,v_history=simulate_euler_maruyama!(q,5e-3,n_steps,β,η; record_hist=true)

xlims=(-2,2)
ylims=(-2,2) =#

#= anim=@animate for i=1:n_steps
    println("$i/$n_steps")
    q_point=q_history[i]
    scatter(q_point[:,1],q_point[:,2],label="",xlims=xlims,ylims=ylims,showaxis=false,aspect_ratio=1,markersize=1)
end
mp4(anim,"thevenin_dynamics.mp4")
plot(v_history,label="")
savefig("thevenin_v.pdf")
 =#

q=1 .- 2rand(N,2)
η_range=0.0:1.0:10.0
responses=Float64[]

n_eq_steps=5000
n_sim_steps=5000000

for η in η_range
    simulate_euler_maruyama!(q,1e-3,n_eq_steps,β,η;record_hist=false)
    println(η)
    r_hist=simulate_euler_maruyama!(q,1e-3,n_sim_steps,β,η;record_hist=true)
    push!(responses,sum(r_hist)/length(r_hist))
    println(sum(r_hist)/length(r_hist))
    flush(STDOUT)
end
f=open("toy_thevenin_profile.txt")
println(f,join(η_range," "))
println(f,join(responses," "))