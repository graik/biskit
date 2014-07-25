function plot_single_rms( dat, color )
%3D plot of receptor and ligand rms versus docking performance

    'single_rms'
    p = plot3( dat(:,4), dat(:,9), dat(:,11) )

    set( p, 'Marker', '.')
    set( p, 'LineStyle', 'none')
    set( p, 'Color', color )

    fig = get( p, 'Parent')
    set( fig, 'YGrid', 'on')
    set( fig, 'XGrid', 'on')
    set( fig, 'ZGrid', 'on')

    label = get( fig, 'ZLabel')
    set( label, 'String', 'docking performance' )

    label = get( fig, 'XLabel')
    set( label, 'String', 'interface rms rec f->b' )

    label = get( fig, 'YLabel')
    set( label, 'String', 'interface rms lig f->b' )

