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

## last $Author$
## last $Date$
## $Revision$
##
## Contributions:
## mostly copied over from the MatrixPlot class of Wolfgang Rieping
"""
Create color scales.
"""

import Numeric as N

from Errors import BiskitError

class ColorError( BiskitError ):
    pass

class ColorSpectrum:
    """
    Translate a range of numeric values into a range of color codes.

    Example::
     >>> p = ColorSpectrum( 'grey', 1, 500 )
     >>> single_color= p.color( 250 )
     >>> color_range = p.colors( range(25,250), resetLimits=0 )

    Available palettes are:
    * grey
    * plasma
    * plasma2 (default)
    * sausage (seems to have a discontinuity at 50%, see example below)
    """
    
    MAX_COL = {'grey': 3 * 255,
               'plasma': 3 * 255,
               'plasma2': 3 * 255,
               'sausage': 2 * 255}


    def __init__(self, palette="plasma2", vmin=0., vmax=1. ):
        """
        Create a new palette of given type and value range.
        
        @param palette: palette type (grey, sausage, plasma, plasma2)
                        (default: plasma2)
        @type  palette: str
        @param vmin: smallest value covered by the color range (default: 0.)
        @type  vmin: float
        @param vmax: largest value covered by the color range (default: 1.)
        @type  vmax: float

        @raise ColorError: if palette unknown
        """

        try:
            self.col_func = eval( 'self._ColorSpectrum__' + palette )

            self.vmin = vmin * 1.0
            self.vmax = vmax * 1.0
            self.col_max = self.MAX_COL[ palette ]

        except AttributeError:
            raise ColorError, 'Unknown palette: ' + str(palette)

        except IndexError, why:
            raise ColorError, 'Undefined palette: ' + str(why)


    def __make_col(self, red, green, blue):
        """
        Create color.
        
        @param red: rgb color, 0-255
        @type  red: int
        @param green: rgb color, 0-255
        @type  green: int
        @param blue: rgb color, 0-255
        @type  blue: int

        @return: color
        @rtype: int
        """
        return ((red << 16) + (green << 8) + blue)        


    def __normalize( self, value ):
        """
        Normalize values

        @param value: normalization value
        @type  value: float

        @return: normalized color
        @rtype: int 
        """
        return (value - self.vmin) / ( self.vmax - self.vmin ) * self.col_max


    def color(self, value ):
        """
        Translate a single value into a color.
        
        @param value: value to be translated into color
        @type  value: float
        @return: color code for value
        @rtype: int
        """
        r = self.__make_col( *self.col_func(self.__normalize(value)) )

        return r


    def colors( self, values, resetLimits=1 ):
        """
        Translate a list of values into a list of colors.
        
        @param values: values to be translated into colors
        @type  values: [float]
        @param resetLimits: re-define color range on max and min of values
                            (default: 1)
        @type  resetLimits: 1|0
        
        @return: color codes
        @rtype: [int]
        """
        if resetLimits:
            self.vmax = max( values )
            self.vmin = min( values )

        return [ self.color(v) for v in values ]


    def color_array( self, a, resetLimits=1 ):
        """
        @param a: array of float
        @type  a: array of float
        @param resetLimits: re-define color range on max and min of values
                            (default: 1)
        @type  resetLimits: 1|0
        
        @return: matrix of color codes with same dimensions as a
        @rtype: array of float
        """
        s = N.shape( a )
        v = N.ravel( a )

        r = self.colors( v, resetLimits=resetLimits )

        r = N.reshape( r, s )

        return r


    def legend( self ):
        """
        @return: color mapping for each color
        @rtype: [ (float,int) ], value
        """
        r = []
        step = (self.vmax - self.vmin) / self.col_max

        for i in range( self.col_max ):

            v = i*step + self.vmin
            c = self.color( v )

            r.append( (v,c) )

        return r  

    def __map(self, x):
        return int(min(255, round(x)))

    ##
    ## Available color palettes
    ##

    def __grey(self, x):
        x = self.__map(255 * x / self.col_max)
        return x, x, x

    def __plasma2(self, x):

        blue_range = 150
        red_range = 255
        green_range = 255

        if x <= blue_range:
            red = 0
            green = self.__map(x)
            blue = self.__map(blue_range - x)

        elif x <= 255:
            red = 0
            blue = 0
            green = self.__map(x)

        elif x > 255 + green_range:
            x -= 255 + green_range
            blue = 0
            red = 255
            green = self.__map(255 - x)
        else:
            x -= 255
            blue = 0
            red = self.__map(x)
            green = 255

        return red, green, blue

    def __plasma(self, x):

        blue_range = 255
        red_range = 255
        green_range = 255

        if x <= blue_range:
            red = 0
            green = self.__map(x)
            blue = self.__map(blue_range - x)

        elif x > 255 + green_range:
            x -= 255 + green_range
            blue = 0
            red = 255
            green = self.__map(255 - x)
        else:
            x -= 255
            blue = 0
            red = self.__map(x)
            green = 255

        return red, green, blue

    def __sausage(self, x):

        first_half = 255

        if x <= first_half:
            red = self.__map(x)
            green = 0
            blue = self.__map(first_half - x)

        else:
            x -= 255
            red = 255
            green = self.__map(x)
            blue = int(0.3 * 255)

        return red, green, blue

#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0

        @return: grey color spectrum
        @rtype:  [ int ]
        """
        import biggles as B
        
        c_grey    = ColorSpectrum( 'grey', 0, 100 )
        c_sausage = ColorSpectrum( 'sausage', 0, 100 )
        c_plasma  = ColorSpectrum( 'plasma', 0, 100 )
        c_plasma2 = ColorSpectrum( 'plasma2', 0, 100 )

        p = B.FramedPlot()
        
        result = []
        for i in range( 100 ):

            x = (i, i+1 )
            p.add( B.FillBelow( x, (1., 1.),
                                color = c_grey.color( i ) ) )
            result += [ c_grey.color( i ) ]
            
            p.add( B.FillBelow( x, (0.75, 0.75),
                                color = c_sausage.color( i ) ) )
            p.add( B.FillBelow( x, (0.5, 0.5),
                                color = c_plasma.color( i ) ) )
            p.add( B.FillBelow( x, (0.25, 0.25),
                                color = c_plasma2.color( i ) ) )

        p.add( B.Curve( (0,100), (1.,1.)) )
        p.add( B.Curve( (0,100), (.75,.75)) )
        p.add( B.Curve( (0,100), (.5,.5) ))
        p.add( B.Curve( (0,100), (0.25, 0.25)) )
        p.add( B.Curve( (0,100), (0.0, 0.0)) )

        p.add( B.PlotLabel(  0.5 ,0.9, 'grey') )
        p.add( B.PlotLabel(  0.5 ,0.65, 'sausage') )
        p.add( B.PlotLabel(  0.5 ,0.4, 'plasma') )
        p.add( B.PlotLabel(  0.5 ,0.15, 'plasma2') )

        if local:
            globals().update( locals() )
            p.show()
            
        return result


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: grey color spectrum
        @rtype:  [ int ]
        """
        return [0, 197379, 328965, 526344, 657930, 855309, 986895, 1184274, 1315860, 1513239, 1710618, 1842204, 2039583, 2171169, 2368548, 2500134, 2697513, 2829099, 3026478, 3158064, 3355443, 3552822, 3684408, 3881787, 4013373, 4210752, 4342338, 4539717, 4671303, 4868682, 5066061, 5197647, 5395026, 5526612, 5723991, 5855577, 6052956, 6184542, 6381921, 6513507, 6710886, 6908265, 7039851, 7237230, 7368816, 7566195, 7697781, 7895160, 8026746, 8224125, 8421504, 8553090, 8750469, 8882055, 9079434, 9211020, 9408399, 9539985, 9737364, 9868950, 10066329, 10263708, 10395294, 10592673, 10724259, 10921638, 11053224, 11250603, 11382189, 11579568, 11776947, 11908533, 12105912, 12237498, 12434877, 12566463, 12763842, 12895428, 13092807, 13224393, 13421772, 13619151, 13750737, 13948116, 14079702, 14277081, 14408667, 14606046, 14737632, 14935011, 15132390, 15263976, 15461355, 15592941, 15790320, 15921906, 16119285, 16250871, 16448250, 16579836]
        
        

if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()
    
    

