function animate_system(sys,filename,F)
    l,l,l=sys.boundary.side_lengths
    map_color(y)=RGB(1-(y+1)/2,0.0,(y+1)/2)#color particles based on y coordinate
    P=values(sys.loggers.coords)
    n_steps=length(P)
    N=length(first(P))
    anim=@animate for i=1:n_steps
        if i%100==0
            println("Frame $(i)/$(n_steps)")
        end
        X=[P[i][j][1] for j=1:N]
        Y=[P[i][j][2] for j=1:N]
        Z=[P[i][j][3] for j=1:N]
        
        plot(X,Y,seriestype=:scatter,color=map_color.(F.(Y)),label="",xlims=(0,l),ylims=(0,l),zlims=(0,l),showaxis=true,ticks=false,camera=(0,90),msw=0,markersize=3)
    end
    mp4(anim,filename,fps=30)

end