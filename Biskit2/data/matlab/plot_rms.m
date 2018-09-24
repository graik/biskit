function plot_rms( f, color )
    echo off
    
    dat = load( f )
    '..plot_rms..'
    
    n = length( dat ) / 121
    colors = colormap( hsv(n))
    
    hold on
    
    for i=1:n
        
        sub = dat( (i-1)*121+1 : i*121, : )
        plot_single_rms( sub, colors(i, :) )
    end
    
    echo on