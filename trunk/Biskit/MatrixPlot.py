##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Wolfgang Rieping
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
"""Plot a 2D matrix (up to 100 x 100)"""

import biggles
import Numeric as N
import ColorSpectrum as C 

class Legend(biggles.FramedPlot):

    def __init__(self, values):

        biggles.FramedPlot.__init__(self)

        values = N.array(values)
        
        self.frame.draw_spine = 1

        n_values = 4
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

class MatrixPlot(biggles.FramedPlot):

    def __init__(self, matrix, mesh = 0, palette = "plasma", legend = 0):

        biggles.FramedPlot.__init__(self)

        self.palette = C.ColorSpectrum( palette )

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

        l = self.palette.legend()

        legend = Legend( l )

        inset = biggles.Inset((1.1, 0.60), (1.2, .97), legend)

        return inset

    
if __name__ == '__main__':

    from Numeric import *

    n = 30

    z = zeros((n,n), Float)

    for i in range(shape(z)[0]):
        for j in range(shape(z)[1]):
            z[i,j] = exp(-0.01*((i-n/2)**2+(j-n/2)**2))

    p = MatrixPlot(z, palette='sausage', legend=1)
    p.show()
