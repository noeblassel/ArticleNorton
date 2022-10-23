using Plots,LinearAlgebra


N=1000
β=5.0


function F(p)
    x,y=p
    r=sqrt(dot(p,p))
    return [y/r,-x/r]
end
function simulate_euler_maruyama!(q,Δt,n_iterations,β; record_hist::Bool=false)

    N=first(size(q))
    shape=size(q)
    σ=sqrt(2Δt/β)

    if record_hist
        #q_record=typeof(q)[]
        r_record=Float64[]
    end

    for i=1:n_iterations
        if record_hist
            # push!(q_record,copy(q))
             push!(r_record,dot(q,q)/2N)
        end
        q.+= -Δt*q + σ*randn(shape) # Euler-Maruyama for overdamped Langevin
    end

    if record_hist
        return q_record, (r_record .- 1/β)
    else
        return
    end

end

q=1 .- 2rand(N,2)
n_steps=100000

q_history,r_history=simulate_euler_maruyama!(q,1e-3,n_steps,β; record_hist=true)

#= xlims=(-2,2)
ylims=(-2,2)

anim=@animate for i=1:n_steps
    println("$i/$n_steps")
    q_point=q_history[i]
    scatter(q_point[:,1],q_point[:,2],label="",xlims=xlims,ylims=ylims,showaxis=false,aspect_ratio=1,markersize=1)
end
mp4(anim,"equilibrium_dynamics.mp4") =#
plot(r_history,label="")
savefig("equilibrium_r.pdf")