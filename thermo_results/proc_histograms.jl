#!/libre/blasseln/julia-1.8.2/bin/julia

n_bins, files...=ARGS

if !isdir("histograms")
    mkdir("histograms")
end

n_bins=parse(Int64,n_bins)

for f in files
    println(f)
    fh=open(f,"r")
    series=reinterpret(Float64,read(fh))
    close(fh)
    m,M=minimum(series),maximum(series)
    dx=(M-m)/n_bins
    hist=zeros(n_bins)

    ts = (M .- series) / (M-m)
    is =floor.(Int64,ts*n_bins)
    map(i-> hist[i]+=1,is)

    output_f=open("histograms/$f","w")
    println(output_f,"$m $M $n_bins")
    print(output_f,join(hist,", "))
    close(output_f)
end
