using Statistics

method = ARGS[1]

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

output_file= (method == "THEVENIN") ? "thevenin_results.txt" : "norton_results.txt"
header = (method =="THEVENIN") ? "N η R T N_samps AV_R AV_T" : "N λ r T N_samps AV_λ AV_T"
f_output=open(output_file,"w")

println(f_output,header)

for Npd in 6:20
    println(Npd)
    S_file= (method == "THEVENIN") ? open("thevenin_response_SINUSOIDAL_0.3_$Npd.out","r") : open("./norton_forcing_SINUSOIDAL_0.1_$Npd.out","r")
    S_series=reinterpret(Float64,read(S_file))
    close(S_file)
    T_file=(method == "THEVENIN") ? open("thevenin_temp_SINUSOIDAL_0.3_$Npd.out","r") : open("./norton_temp_SINUSOIDAL_0.1_$Npd.out","r")
    T_series=reinterpret(Float64,read(T_file))
    close(T_file)
    av_S=asymptotic_variance(S_series)
    av_T=asymptotic_variance(T_series)
    N_samps=length(S_series)
    output_str= (method =="THEVENIN") ? "$(Npd^3) 0.3 $(mean(S_series)) $(mean(T_series)) $(N_samps) $(av_S) $(av_T)" : "$(Npd^3) $(mean(S_series)) 0.1 $(mean(T_series)) $(N_samps) $(av_S) $(av_T)"
    println(f_output,output_str)
end

close(f_output)
#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)