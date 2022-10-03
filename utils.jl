function place_atoms_on_3D_lattice(N_per_dim::Integer,box_size)
    (Lx,Ly,Lz)=box_size.side_lengths 
    reshape([SVector(i * Lx/N_per_dim, j * Ly/N_per_dim, k * Lz/N_per_dim) for i = 0:N_per_dim-1, j = 0:N_per_dim-1, k = 0:N_per_dim-1],N_per_dim^3)
end


function animate_system(sys,filename)
    l,l,l=sys.box_size
    map_color(y)=RGB(1-(y+1)/2,0.0,(y+1)/2)#color particles based on y coordinate
    P=sys.loggers[:coords].coords
    n_steps=length(P)
    N=length(first(P))
    anim=@animate for i=1:n_steps
        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end
        X=[P[i][j][1] for j=1:N]
        Y=[P[i][j][2] for j=1:N]
        Z=[P[i][j][3] for j=1:N]
        
        plot(X,Y,seriestype=:scatter,color=map_color.(f_dict[method].(Y)),label="",xlims=(0,l),ylims=(0,l),zlims=(0,l),showaxis=true,ticks=false,camera=(0,90),msw=0,markersize=3)
    end
    mp4(anim,filename,fps=30)

end