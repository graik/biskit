##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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
"""
bar-plotting for Biggles
"""

try:
    import biggles as B
except:
    B = 0


##############################
## low-level tools


def multibar_curve( values, x0=None, xwidth=0.25, xsep=0.15 ):
    """
    Get x,y values to draw many bars.
    
    @param values: values to bar-plot
    @type  values: [float]
    @param x0: start of first bar
    @type  x0: float
    @param xwidth: width of bars (default: 0.25)
    @type  xwidth: float
    @param xsep: space between bars (default: 0.15)
    @type  xsep: float
    
    @return: x,y values that draw a bar for each y value
    @rtype: [float],[float]
    
    @note: not used at the moment because it draws connecting
           line at the bottom
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
    """
    Get x,y values for a single bar.
    
    @param height: height of bar
    @type  height: float
    @param x0: start of first bar
    @type  x0: float
    @param xwidth: width of bars (default: 0.25)
    @type  xwidth: float
    @param y0: start of first bar
    @type  y0: float
    
    @return: x,y values that draw a single bar.
    @rtype: [ float ], [ float ]
    """
    y = height - y0
    return [ x0, x0, x0+xwidth, x0+xwidth ], [y0,y,y,y0]


def boxed_diagonal( x0, y0, x1, y1, vy ):
    """
    Get a diagonal connecting left (x0) to right (x1) edge but which is only
    visible within a rectangular box.

    @param x0: rectangular window of visibility  
    @type  x0: float
    @param y0: rectangular window of visibility
    @type  y0: float
    @param x1: rectangular window of visibility
    @type  x1: float
    @param y1: rectangular window of visibility
    @type  y1: float
    
    @param vy: virtual starting point (y-coordinate, may be outside of the box)
    @type  vy: float
    
    @return: xa, ya, xb, yb - start and end point of a diagonal
             line originating at (x0, vy) but beeing only visible
             within the rectangle (x0,y0,x1,y1)
    @rtype: float, float, float, float
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
    Fill the rectangle described by x0, y0, x1, y1 with horizontal lines.
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param sep: separation between lines (default: 0.1)
    @type  sep: float
    @param size: line thickness (default: 0.06)
    @type  size: float
    @param color: color name (default: black)
    @type  color: str
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: list of biggles plot objects
    @rtype: [Biggles.Curve]
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
    Fill the rectangle described by x0, y0, x1, y1 with horizontal bars.
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param sep: separation between lines (default: 0.1)
    @type  sep: float
    @param size: line thickness (default: 0.06)
    @type  size: float
    @param color: color name (default: grey)
    @type  color: str
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: list of biggles plot objects
    @rtype: [Biggles.FillBetween]
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
    Fill the rectangle described by x0,y0, x1,y1 (lower left and
    upper right corner) with diagonal bars.
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param sep: offset between (low edges of) diagonal lines (default: 0.2)
    @type  sep: float
    @param size: width of diagonal line (default: 0.05)
    @type  size: float
    @param invert: mirror diagonals (default: 0)
    @type  invert: 1|0
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: list of biggles plot objects
    @rtype: [Biggles.FillBetween]
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
    Fill the rectangle described by x0,y0, x1,y1 (lower left and
    upper right corner) with diagonal lines.
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param sep: offset between (low edges of) diagonal lines (default: 0.2)
    @type  sep: float
    @param size: width of diagonal line (default: 0.05)
    @type  size: float
    @param invert: mirror diagonals (default: 0)
    @type  invert: 1|0
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: list of biggles plot objects
    @rtype: [Biggles.Curve]
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
    """
    Fill the rectangle described by x0,y0, x1,y1 (lower left and
    upper right corner).
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: list of biggles plot objects
    @rtype: [Biggles.FillBetween]
    """
    return [ B.FillBetween( [x0,x0,x1,x1], [y1,y0,y0,y1], [x0,x1], [y1,y1],
                          filltype=1, **kw) ]


#####################
## Bar-plotting


def box_curve( x0, y0, x1, y1, **kw ):
    """
    A rectangle described by x0,y0, x1,y1 (lower left and
    upper right corner).
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: biggles plot object
    @rtype: biggles.Curv
    """
    return B.Curve( [x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0], **kw )


def box_fill( x0, y0, x1, y1, **kw ):
    """
    Fill for a rectangle described by x0,y0, x1,y1 (lower left and
    upper right corner).
    
    @param x0: rectangular coorinates
    @type  x0: float
    @param y0: rectangular coorinates
    @type  y0: float
    @param x1: rectangular coorinates
    @type  x1: float
    @param y1: rectangular coorinates
    @type  y1: float
    @param kw: additional key-value pairs
    @type  kw: key=value
    
    @return: biggles plot object
    @rtype: biggles.FillBetween
    """
    return B.FillBetween( [x0,x1,x1],[y0,y0,y1],[x1,x0,x0],[y1,y1,y0], **kw)


def add_box( p, x0, y0, x1, y1, fillfunc=solid_fill, **kw ):
    """
    Add single filled box to plot.

    @param p: plot object which to add rectangle to
    @type  p: biggles.FramedPlot
    @param x0: rectangular coorinates, lower left corner
    @type  x0: float
    @param y0: rectangular coorinates, lower left corner
    @type  y0: float
    @param x1: rectangular coorinates, upper right corner
    @type  x1: float
    @param y1: rectangular coorinates, upper right corner
    @type  y1: float
    @param fillfunc: function name (default: solid_fill)
    @type  fillfunc  str
    @param kw: arguments for fillfunc and (starting with 'l') for Curve
    @type  kw: key=value
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
    
    @param p: fill objects are added directly to the plot
    @type  p: Biggles.FramedPlot
    @param values: bar heights
    @type  values: [float]
    @param x0: left edge of first bar (default: None, xwidth)
    @type  x0: float
    @param xwidth: width of each bar (default: 0.25)
    @type  xwidth: float
    @param xoffset: distance between bars (measured between left edges)
                   (default: 1)
    @type  xoffset: float
    @param fillfunc: function solid_fill | line_fill | diagonal_fill | bar_fill
    @type  fillfunc: str
    @param kw: attributes for fill function, e.g. color, size, sep
    @type  kw: key=value
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
    
    @param p: fill objects are added directly to the plot
    @type  p: Biggles.FramedPlot
    @param values: bar heights
    @type  values: [float]
    @param x0: left edge of first bar (default: None, xwidth)
    @type  x0: float
    @param xwidth: width of each bar  (default 0.25)
    @type  xwidth: float
    @param xoffset: distance between bars (measured between left edges)
                    (default 1)
    @type  xoffset: float
    @param fillfunc: function solid_fill | line_fill | diagonal_fill
                              | bar_fill | None
    @type  fillfunc: str
    @param kw: attributes for fill function (e.g. color, size, sep)
               and Curve (the latter ones have to start with "l", e.g. lcolor)
    @type  kw: key=value
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
    """
    Initiate a biggles.FramedPlot object.

    @param xlabel: label for x-axis 
    @type  xlabel: str
    @param ylabel: label for y-axis 
    @type  ylabel: str
    @param yrange: range of y-axis
    @type  yrange: (float,float)
    @param xrange: range of x-axis
    @type  xrange: (float,float)
    @param width: hard plot width (in pixels or cm)
    @type  width: int
    @param height: hard plot height
    @type  height: int

    @return: biggles plot object
    @rtype: biggles.FillBetween    
    """
    if not B:
        raise ImportError('biggles is required for this method')
    
    B.configure( 'screen', 'height', height )
    B.configure( 'screen', 'width', width )

    p = B.FramedPlot()

    p.xlabel = xlabel
    p.ylabel = ylabel
    if yrange: p.yrange = yrange
    if xrange: p.xrange = xrange

    return p



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    TAGS = [BT.OLD, BT.BIGGLES]

    def test_plotUtils(self):
        """plotUtils test"""
        self.p = prepare_plot(xlabel='', ylabel='flex $\\langle{x}\\rangle$',
                              xrange=(0,5), yrange=(0,4) )

        add_bars( self.p, [ 1, 2.5, 1.25, 0.3 ], fillfunc=diagonal_fill,
                  color='grey', size=0.1, invert=1,
                  lcolor='black', lwidth=1  )

        add_bars( self.p, [0.4, 1.5, 2., 0.6], x0=0.5,
                  fillfunc=solid_fill, lcolor='blue', lwidth=2 )

        add_box( self.p, 3.5, 3, 4, 3.5, fillfunc=line_fill,
                 color='grey', size=5 )

##         ## example of how to write an eps plot to disc
##         p.write_eps(T.absfile('~/test.eps'), width='10cm', height='8.7cm')
        
        if self.local or self.VERBOSITY > 2:
            self.p.show()
        
       

if __name__ == '__main__':

    BT.localTest()


