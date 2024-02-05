## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

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
## Contributions:
## mostly copied over from the MatrixPlot class of Wolfgang Rieping
"""
Create color scales.
"""

import biskit.core.oldnumeric as N0

from biskit import BiskitError


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


    def __init__(self, palette="plasma2", vmin=0., vmax=1., default=0xffffff ):
        """
        Create a new palette of given type and value range.
        
        :param palette: palette type (grey, sausage, plasma, plasma2)
                        (default: plasma2)
        :type  palette: str
        :param vmin: smallest value covered by the color range (default: 0.)
        :type  vmin: float
        :param vmax: largest value covered by the color range (default: 1.)
        :type  vmax: float
        :param default: default color for values below the palette range
        :type  default: int (color code, default=0xffff)

        :raise ColorError: if palette unknown
        """

        try:
            self.col_func = eval( 'self._ColorSpectrum__' + palette )

            self.vmin = vmin * 1.0
            self.vmax = vmax * 1.0
            self.col_max = self.MAX_COL[ palette ]
            self.default_color = default

        except AttributeError:
            raise ColorError('Unknown palette: ' + str(palette))

        except IndexError as why:
            raise ColorError('Undefined palette: ' + str(why))


    def __make_col(self, red, green, blue):
        """
        Create color.
        
        :param red: rgb color, 0-255
        :type  red: int
        :param green: rgb color, 0-255
        :type  green: int
        :param blue: rgb color, 0-255
        :type  blue: int

        :return: color
        :rtype: int
        """
        return ((red << 16) + (green << 8) + blue)        


    def __normalize( self, value ):
        """
        Normalize values

        :param value: normalization value
        :type  value: float

        :return: normalized color
        :rtype: int 
        """
        if self.vmax == self.vmin:
            return self.col_max
        return (value - self.vmin) / ( self.vmax - self.vmin ) * self.col_max


    def color(self, value ):
        """
        Translate a single value into a color.
        
        :param value: value to be translated into color
        :type  value: float
        :return: color code for value
        :rtype: int
        """
        if value < self.vmin:
            return self.default_color
        
        r = self.__make_col( *self.col_func(self.__normalize(value)) )

        return r


    def colors( self, values, resetLimits=1 ):
        """
        Translate a list of values into a list of colors.
        
        :param values: values to be translated into colors
        :type  values: [float]
        :param resetLimits: re-define color range on max and min of values
                            (default: 1)
        :type  resetLimits: 1|0
        
        :return: color codes
        :rtype: [int]
        """
        if resetLimits:
            self.vmax = max( values ) * 1.
            self.vmin = min( values ) * 1.

        return [ self.color(v) for v in values ]


    def color_array( self, a, resetLimits=1 ):
        """
        :param a: array of float
        :type  a: array of float
        :param resetLimits: re-define color range on max and min of values
                            (default: 1)
        :type  resetLimits: 1|0
        
        :return: matrix of color codes with same dimensions as a
        :rtype: array of float
        """
        s = N0.shape( a )
        v = N0.ravel( a )

        r = self.colors( v, resetLimits=resetLimits )

        r = N0.reshape( r, s )

        return r


    def legend( self ):
        """
        :return: color mapping for each color
        :rtype: [ (float,int) ], value
        """
        r = []
        step = (self.vmax - self.vmin) // self.col_max

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


def colorRange( nColors, palette='plasma2' ):
    """Quick access to a range of colors.

    :param nColors: number of colors needed
    :type  nColors: int
    :param palette: type of color spectrum
    :type  palette: str
    :return: a range of color values
    :rtype: [ int ]
    """
    c = ColorSpectrum( palette=palette, vmin=0, vmax=1. )

    r = 1. * N0.arange( 0, nColors ) / nColors

    return c.colors( r )



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """ColorSpectrum test"""
    
    def test_ColorSpectrum( self ):
        """ColorSpectrum test"""
        try:
            import biskit.tools as T
            import biggles as B
        except:
            B = 0
        
        c_grey    = ColorSpectrum( 'grey', 0, 100 )
        c_sausage = ColorSpectrum( 'sausage', 0, 100 )
        c_plasma  = ColorSpectrum( 'plasma', 0, 100 )
        c_plasma2 = ColorSpectrum( 'plasma2', 0, 100 )

        if B:
            self.p = B.FramedPlot()

##        old_spectrum = T.colorSpectrum( 100 )
        
        self.result = []
        for i in range( -1, 100 ):

            x = (i, i+1 )

            if B:
                self.result += [ c_grey.color( i ) ]
    
                self.p.add( B.FillBelow( x, (1., 1.),
                                         color = c_grey.color( i ) ) )
                
                self.p.add( B.FillBelow( x, (0.75, 0.75),
                                         color = c_sausage.color( i ) ) )
                self.p.add( B.FillBelow( x, (0.5, 0.5),
                                         color = c_plasma.color( i ) ) )
                self.p.add( B.FillBelow( x, (0.25, 0.25),
                                         color = c_plasma2.color( i ) ) )

##                self.p.add( B.FillBelow( x, (0., 0.),
##                                      color = old_spectrum[i] ))

        if B:
            self.p.add( B.Curve( (0,100), (1.,1.)) )
            self.p.add( B.Curve( (0,100), (.75,.75)) )
            self.p.add( B.Curve( (0,100), (.5,.5) ))
            self.p.add( B.Curve( (0,100), (0.25, 0.25)) )
            self.p.add( B.Curve( (0,100), (0.0, 0.0)) )
    
            self.p.add( B.PlotLabel(  0.5 ,0.9, 'grey') )
            self.p.add( B.PlotLabel(  0.5 ,0.65, 'sausage') )
            self.p.add( B.PlotLabel(  0.5 ,0.4, 'plasma') )
            self.p.add( B.PlotLabel(  0.5 ,0.15, 'plasma2') )

        if (self.local or self.VERBOSITY > 2) and B:
            self.p.show()
            
        ##self.assertEqual(self.result, self.EXPECTED)
        ## tolerate two differences to account for Python 3 result
        if B:
            a = N0.array(self.result)
            b = N0.array(self.EXPECTED)
            self.assert_(N0.count_nonzero(a-b)<3)

    EXPECTED = [16777215, 0, 197379, 328965, 526344, 657930, 855309, 986895, 1184274, 1315860, 1513239, 1710618, 1842204, 2039583, 2171169, 2368548, 2500134, 2697513, 2829099, 3026478, 3158064, 3355443, 3552822, 3684408, 3881787, 4013373, 4210752, 4342338, 4539717, 4671303, 4868682, 5066061, 5197647, 5395026, 5526612, 5723991, 5855577, 6052956, 6184542, 6381921, 6513507, 6710886, 6908265, 7039851, 7237230, 7368816, 7566195, 7697781, 7895160, 8026746, 8224125, 8421504, 8553090, 8750469, 8882055, 9079434, 9211020, 9408399, 9539985, 9737364, 9868950, 10066329, 10263708, 10395294, 10592673, 10724259, 10921638, 11053224, 11250603, 11382189, 11579568, 11776947, 11908533, 12105912, 12237498, 12434877, 12566463, 12763842, 12895428, 13092807, 13224393, 13421772, 13619151, 13750737, 13948116, 14079702, 14277081, 14408667, 14606046, 14737632, 14935011, 15132390, 15263976, 15461355, 15592941, 15790320, 15921906, 16119285, 16250871, 16448250, 16579836]
        
if __name__ == '__main__':

    BT.localTest()
    

