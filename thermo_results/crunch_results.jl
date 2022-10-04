using Statistics

function asymptotic_var(v)
    data=copy(v)
    avg=mean(v)
    data .-= avg
    N_steps=floor(Int64,log2(length(data)))
    data=data[1:2^N_steps]#crop series to largest possible power of two
    L=1
    N=length(data)
    K=N
    max_var=0.0
    while K>1000
        max_var=max(max_var,varm(data,0.0)*N*inv(K))
        K >>=1
            data=(data[1:2:end]+data[2:2:end])/2
        end
    return max_var
end

output_file="norton_forcing_SINUSOIDAL_0.1.dat"
f_output=open(output_file,"w")
println(f_output,"N r λ T N_samps AV_λ AV_T")
println(n)

for Npd in 1:15
    λ_file=open("norton_forcing_SINUSOIDAL_0.1_$Npd.dat","r")
    λ_series=reinterpret(Float64,read(λ_file))
    close(λ_file)
    T_file=open("norton_temp_SINUSOIDAL_0.1_$Npd")
    T_series=reinterpret(Float64,read(T_file))
    close(T_file)
    av_λ=asymptotic_variance(λ_series)
    av_T=asymptotic_variance(T_series)
    N_samps=length(λ_series)

    join(f_output,[Npd^3, 0.1, mean(λ_series),mean(T_series),N_samps,av_λ,av_T]," ","\n")
end

close(f_output)
#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)