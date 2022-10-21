using Statistics

println(ARGS)
n_bins,output_file,input_files... = ARGS
n_bins=parse(Int64,n_bins)

if !isdir("histograms")
    mkdir("histograms")
end

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

header = "FILENAME MEAN N_SAMPLES ASYMPTOTIC_VARIANCE HISTOGRAM_PATH"
f_output=open(output_file,"w")

println(f_output,header)

for f in input_files
    series=reinterpret(Float64,read(f))
    av=asymptotic_variance(series)
    avg=mean(series)
    N_samps=length(series)
    println(f_output,"$f $avg $N_samps $av $(ENV["PWD"])/histograms/histogram_$f")

    # compute histogram
    m,M=minimum(series),maximum(series)
    dx=(M-m)/n_bins
    hist=zeros(n_bins)
    ts = (M .- series) / (M-m)
    is = ceil.(Int64,ts*n_bins)
    clamp!(is,1,n_bins)
    println(m," ", M, " ",minimum(is)," ",maximum(is))
    map(i-> hist[i]+=1,is)

    hist_file=open("histograms/histogram_$(n_bins)_$f")
    println(hist_file,"$m $M $n_bins")
    print(hist_file,join(hist,", "))
    close(hist_file)
end

close(f_output)
#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)