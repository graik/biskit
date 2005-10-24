##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##
## last $Author$
## last $Date$
## $Revision$
"""bar-plotting for Biggles"""

import biggles as B

##############################
## low-level tools

def multibar_curve( values, x0=None, xwidth=0.25, xsep=0.15 ):
    """
    -> ([float],[float]), x,y values that draw a bar for each y value
    not used at the moment because it draws connecting line at the bottom
    """
    x  = x0 or xsep
    vy = []
    vx = []
    
    for y in values:
        vy += [ 0, y,y, 0]
        vx += [ x, x, x+xwidth, x+xwidth]
        x += xwidth + xsep

    return vx, vy

def bar_curve( height, x0, xwidth=0.25, y0=0. ):
    """ -> ([ float ], [ float ]) - x,y values that draw a single bar."""
    y = height - y0
    return [ x0, x0, x0+xwidth, x0+xwidth ], [y0,y,y,y0]
    

def boxed_diagonal( x0, y0, x1, y1, vy ):
    """
    Get a diagonal connecting left (x0) to right (x1) edge but which is only
    visible within a rectangular box.
    x0,y0, x1,y1 - rectangular window of visibility
    vy - virtual starting point (y-coordinate, may be outside of the box)
    -> xa, ya, xb, yb - start and end point of a diagonal line originating at
       (x0, vy) but beeing only visible within the rectangle (x0,y0,x1,y1)
    """
    dx = x1 - x0
    dy = y1 - y0

    ## line completely within rectangle
    xa = x0
    ya = vy
    xb = x1
    yb = vy + dx

    ## cut off low end
    if vy < y0:
        xa = x0 + (y0-vy)
        ya = y0 ## mhm, doesn't work as it is supposed to

    ## cut of high end
    if vy > y1 - dx:
        xb = x0 + dy - vy
        yb = y1 

    return xa,ya, xb,yb


##########################
## Fill patterns

def line_fill( x0, y0, x1, y1, sep=0.1, size=0.06, color='black', **kw ):
    """
    -> [ Biggles.Curve ], horizontal lines filling the given rectangle
    """
    if not 'color' in kw:
        kw['color'] = color
    if not 'width' in kw:
        kw['width'] = size
        
    r = []
    y = y1
    while y >= y0:

        r += [ B.Curve( [x0,x1], [y,y], **kw ) ]
        y -= sep
        
    return r

def bar_fill( x0, y0, x1, y1, sep=0.1, size=0.06, color='grey', **kw ):
    """
    -> [ Biggles.FillBetween ], horizontal bars filling the given rectangle
    """
    kw.update( {'color':color} )
    r = []
    y = y1
    while y >= y0:

        y_low = y - size
        if y_low < y0:
            y_low = y0

        r += [ B.FillBetween( [x0,x1], [y,y], [x0,x0,x1,x1], [y,y_low,y_low,y],
                              **kw ) ]
        y -= sep
        
    return r

def diagonal_fill( x0, y0, x1, y1, sep=0.2, size=0.05, invert=0, **kw ):
    """
    x0,y0,x1,y1 - lowerleft and upper right corner of rectangle to be filled
    sep    - float, offset between (low edges of) diagonal lines [0.2]
    size   - float, width of diagonal lines [0.05]
    invert - 1|0, mirror diagonals [0]
    -> [ Biggles.FillBetween ], set of diagonal bars filling rectangle
    """
    r = [] ## will hold result

    vy = y0 - (x1 - x0) - sep ## virtual y coordinate

    ## shift mirrored lines horizontally to get better criss-cross effect
    if invert: vy = vy - sep/2

    while vy < (y1 - y0):

        ## lower diagonal
        xa,ya, xb,yb = boxed_diagonal( x0,y0, x1,y1, vy )

        ## upper diagonal
        if vy + size < y1 - y0:
            xc,yc, xd,yd = boxed_diagonal( x0,y0, x1,y1, vy+size )
        else: ## never start outside box
            xc, yc = xa, y1-y0
            xd, yd = xb, y1-y0

        ## mirror line
        if invert:
            xa = x1 - (xa - x0)
            xb = x0 - (xb - x1)
            xc = x1 - (xc - x0)
            xd = x0 - (xd - x1)
        
        r += [ B.FillBetween( [xc,xa,xb,xd],[yc,ya,yb,yd], [xc,xd],[yc,yd],
                              **kw)]
        vy += sep

    return r


def diagonal_line_fill( x0, y0, x1, y1, sep=0.2, size=0.05, invert=0, **kw ):
    """
    x0,y0,x1,y1 - lowerleft and upper right corner of rectangle to be filled
    sep    - float, offset between (low edges of) diagonal lines [0.2]
    size   - float, width of diagonal lines [0.05]
    invert - 1|0, mirror diagonals [0]
    -> [ Biggles.Curve ], set of diagonal lines filling rectangle
    """
    r = [] ## will hold result

    vy = y0 - (x1 - x0) - sep ## virtual y coordinate

    kw[ 'width' ] = kw.get('width',None) or size
    
    ## shift mirrored lines horizontally to get better criss-cross effect
    if invert: vy = vy - sep/2

    while vy <= (y1 - y0):

        ## lower diagonal
        xa,ya, xb,yb = boxed_diagonal( x0,y0, x1,y1, vy )

        ## mirror line
        if invert:
            xa = x1 - (xa - x0)
            xb = x0 - (xb - x1)
        
        r += [ B.Curve( [xa,xb],[ya,yb], **kw) ]
        vy += sep

    return r

def solid_fill( x0, y0, x1, y1, **kw ):
    """ -> [ Biggles.FillBetween ], list with single filled area """
    return [ B.FillBetween( [x0,x0,x1,x1], [y1,y0,y0,y1], [x0,x1], [y1,y1],
                          filltype=1, **kw) ]


#####################
## Bar-plotting

def box_curve( x0, y0, x1, y1, **kw ):
    """-> biggles.Curve, describing rectangle"""
    return B.Curve( [x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0], **kw )

def box_fill( x0, y0, x1, y1, **kw ):
    """-> biggles.FillBetween"""
    return B.FillBetween( [x0,x1,x1],[y0,y0,y1],[x1,x0,x0],[y1,y1,y0], **kw)
    
def add_box( p, x0, y0, x1, y1, fillfunc=solid_fill, **kw ):
    """
    Add single filled box to plot.
    p        - biggles.FramedPlot
    x0, y0   - float, lower left corner
    x1, y1   - float, upper right corner
    fillfunc - function
    key=value arguments for fillfunc and (starting with 'l') for Curve
    """

    for f in fillfunc( x0,y0, x1,y1, **kw ):
        p.add( f )

    ## extract attributes starting with 'l'
    attr = {}
    for k,v in kw.items():
        if k[0] == 'l':
            attr[ k[1:] ] = v

    p.add( box_curve( x0,y0, x1,y1, **attr ) )

    
def fill_bars( p, values, x0=None, xwidth=0.25, xoffset=1, fillfunc=solid_fill,
               margin=0.01, **kw ):
    """
    Add filling to bars described y values, left-most x, width and offset. 
    p        - Biggles.FramedPlot, fill objects are added directly to the plot
    values   - [ float ], bar heights
    x0       - float, left edge of first bar [xwidth]
    xwidth   - float, width of each bar [0.25]
    offset   - float, distance between bars (measured between left edges) [1]
    fillfunc - function, solid_fill | line_fill | diagonal_fill | bar_fill
    **kw     - key=value, attributes for fill function, e.g. color, size, sep
    """
    vx, vy = multibar_curve( values, x0=x0 or xwidth, xwidth=xwidth,
                             xsep=xoffset-xwidth )

    for i in range( 0, len(vy), 4 ):

        fill_objects = fillfunc( vx[i]+margin,vy[i]+margin,
                                 vx[i+2]-margin,vy[i+2]-margin, **kw )

        for l in fill_objects: p.add( l )


def add_bars( p, values, x0=None, xwidth=0.25, xoffset=1, fillfunc=solid_fill,
              **kw ):
    """
    Add bars to plot, described by y values, left-most x, width and offset. 
    p        - Biggles.FramedPlot, fill objects are added directly to the plot
    values   - [ float ], bar heights
    x0       - float, left edge of first bar [xwidth]
    xwidth   - float, width of each bar      [0.25]
    offset   - float, distance between bars (measured between left edges) [1]
    fillfunc - func, solid_fill | line_fill | diagonal_fill | bar_fill | None
    **kw     - key=value, attributes for fill function (e.g. color, size, sep)
               and Curve (the latter ones have to start with "l", e.g. lcolor)
    """
    if x0 is None:
        x0 = xwidth

    ## draw filling
    if fillfunc is not None:
        fill_bars( p, values,x0=x0, xwidth=xwidth, xoffset=xoffset,
                   fillfunc=fillfunc, **kw)

    ## extract attributes starting with 'l'
    attr = {}
    for k,v in kw.items():
        if k[0] == 'l':
            attr[ k[1:] ] = v

    ## draw bar outline
    for i in range( len( values ) ):

        x = x0 + i * xoffset
        p.add( B.Curve( *bar_curve( values[i], x, xwidth=xwidth  ), **attr) )


######################
## other stuff

def prepare_plot( xlabel='', ylabel='', yrange=None, xrange=None,
                  width=500, height=350 ):

    B.configure( 'screen', 'height', height )
    B.configure( 'screen', 'width', width )

    p = B.FramedPlot()

    p.xlabel = xlabel
    p.ylabel = ylabel
    if yrange: p.yrange = yrange
    if xrange: p.xrange = xrange

    return p

#####################
## Testing / Example

def plot_test():

    p = prepare_plot(xlabel='', ylabel='flex $\langle{x}\rangle$',
                     xrange=(0,5), yrange=(0,4) )
    
    add_bars( p, [ 1, 2.5, 1.25, 0.3 ], fillfunc=diagonal_fill,
              color='grey', size=0.1, invert=1,
              lcolor='black', lwidth=1  )

    add_bars( p, [0.4, 1.5, 2., 0.6], x0=0.5,
              fillfunc=solid_fill, lcolor='blue', lwidth=2 )

    add_box( p, 3.5, 3, 4, 3.5, fillfunc=line_fill,
             color='grey', size=5 )


    return p

if __name__ == '__main__':

    p = plot_test()
    p.write_eps(T.absfile('~/test.eps'), width='10cm', height='8.7cm')

    p.show()
    
