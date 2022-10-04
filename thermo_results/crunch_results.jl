using Statistics

function asymptotic_variance(v)
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

output_file="thevenin_results.txt"
f_output=open(output_file,"w")
println(f_output,"N η R T N_samps AV_R AV_T")

for Npd in 1:15
    println(Npd)
    R_file=open("../thevenin_response_SINUSOIDAL_0.3_$Npd.out","r")
    R_series=reinterpret(Float64,read(R_file))
    close(R_file)
    T_file=open("../thevenin_temp_SINUSOIDAL_0.3_$Npd.out","r")
    T_series=reinterpret(Float64,read(T_file))
    close(T_file)
    av_R=asymptotic_variance(R_series)
    av_T=asymptotic_variance(T_series)
    N_samps=length(R_series)
    println(f_output,"$(Npd^3) $(0.3) $(mean(R_series)) $(mean(T_series)) $(N_samps) $(av_R) $(av_T)")
end

close(f_output)
#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)