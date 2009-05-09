## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Wolfgang Rieping
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

## last $Author$
## last $Date$
## $Revision$
"""
Plot a 2D matrix (up to 100 x 100)
"""

import numpy.oldnumeric as N

from Biskit import ColorSpectrum, EHandler
import Biskit.tools as T

T.tryImport( 'biggles', 'FramedPlot', namespace=globals() )
try:
    import biggles
except:
    biggles = False

class Legend(FramedPlot):
    """
    Class to create a legend to use with a Matrix plot.
    """

    def __init__(self, values):
        """
        @param values: color mapping for each color used
        @type  values: [(float, int)]
        """
        if not biggles:
            raise ImportError, 'biggles module could not be imported.'

        FramedPlot.__init__(self)

        values = N.array(values)

        self.frame.draw_spine = 1

        n_values = 4 ## number of labeled ticks in legend
        step = len(values) / (n_values - 1) + 1

        indices = range(0, len(values), step)
        indices.append(len(values) - 1)

        labels = ['%.1f' % values[i, 0] for i in indices]

        self.y.ticks = len(labels)
        self.y1.ticklabels = labels
        self.y2.draw_ticks = 0
        self.x.draw_ticks = 0
        self.x.ticklabels = []

        i = 2
        x = (2, 3)

        for value, color in values:

            y1 = (i, i)
            y2 = (i + 1, i + 1)

            cell = biggles.FillBetween(x, y1, x, y2, color = int(color))
            self.add(cell)

            i += 1


class MatrixPlot(FramedPlot):
    """
    Class to plot the values of a matix, the rows and the columns
    will be plotted along the x- and y-axis, respectively. The value
    of each cell will be illutrated using the selected color range.
    """

    def __init__(self, matrix, mesh=0, palette="plasma", legend=0, step=1):
        """
        @param matrix: the 2-D array to plot
        @type  matrix: array
        @param mesh: create a plot with a dotted mesh
        @type  mesh: 1|0
        @param palette: color palette name see L{Biskit.ColorSpectrum}
        @type  palette: str
        @param legend: create a legend (scale) showing the walues of the
                       different colors in the plot.  
        @type  legend: 1|0
        @param step: reduce matrix -- take only each step position in x and y
        @type  step: int

        @return: biggles plot object, view with biggles.FramedPlot.show() or
                 save with biggles.FramedPlot.write_eps(file_name).
        @rtype: biggles.FramedPlot
        """
        if not biggles:
            raise ImportError, 'biggles module could not be imported.'

        FramedPlot.__init__(self)

        self.palette = ColorSpectrum( palette )

        if step != 1:
            matrix = self.__thinarray( matrix, step )

        self.matrix = self.palette.color_array( matrix )
        s = N.shape( self.matrix )

        for i in range(s[0]):
            for j in range(s[1]):

                col = self.matrix[i,j]

                x1 = (j, j + 1)
                y1 = (i, i)
                y2 = (i + 1, i + 1)

                cell = biggles.FillBetween(x1, y1, x1, y2, color = col)

                self.add(cell)

        if mesh:

            for i in range(s[0] + 1):
                self.add(biggles.LineY(i, linetype='dotted'))

            for i in range(s[1] + 1):
                self.add(biggles.LineX(i, linetype='dotted'))

        if legend:

            legend = self.__make_legend()

            self.add(legend)

            self.add(biggles.PlotBox((-0.17, -0.1), (1.25, 1.1)))

        self.aspect_ratio = 1.0


    def __make_legend(self):
        """
        Create and position the legend.

        @return: biggles legend object
        @rtype: biggles.Inset
        """
        l = self.palette.legend()

        legend = Legend( l )

        inset = biggles.Inset((1.1, 0.60), (1.2, .97), legend)

        return inset


    def __thinarray( self, a, step ):
        """
        @param a: input array
        @type  a: N.array
        @param step: stepping in both dimensions
        @type  step: int
        @return: smaller array
        @rtype: N.array
        """
        r = N.take( a, range( 0, len(a), step ), axis=0 )
        r = N.take( r, range( 0, len(r[0]), step ), axis=1 )
        return r


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    def test_MatrixPlot( self ):
        """MatrixPlot test"""
        n = 30

        z = N.zeros((n,n), N.Float)

        for i in range(N.shape(z)[0]):
            for j in range(N.shape(z)[1]):
                z[i,j] = N.exp(-0.01*((i-n/2)**2+(j-n/2)**2))

        self.p = MatrixPlot(z, palette='sausage', legend=1)

        if self.local or self.VERBOSITY > 2:
            self.p.show()

        self.assert_( self.p is not None )



if __name__ == '__main__':

    BT.localTest()

