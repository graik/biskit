## Plot a 2D matrix (up to 100 x 100)
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Wolfgang Rieping; All rights reserved
##
## create a histigram from data
## 
## Author: Wolfgang Rieping
##
## last $Author$
## last $Date$
## $Revision$


import biggles
import Numeric as N

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

    max_col = {'grey': 3 * 255,
               'plasma': 3 * 255,
               'plasma2': 3 * 255,
               'sausage': 2 * 255}

    def __init__(self, matrix, mesh = 0, palette = "plasma", legend = 0):

        biggles.FramedPlot.__init__(self)

        f_pal = {'grey': self.__grey,
                 'plasma': self.__plasma,
                 'plasma2': self.__plasma2,
                 'sausage': self.__sausage}

        self.palette = palette
        self.f_pal = f = f_pal[palette]

        self.min, self.max, matrix = self.__normalize(matrix)
        self.m = matrix
        s = N.shape(matrix)

        for i in range(s[0]):
            for j in range(s[1]):

                x = int(round(matrix[i,j]))

                col = self.__make_col(*f(x))

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

        max_col = self.max_col[self.palette]

        step = (self.max - self.min) / max_col
        values = N.arange(self.min, self.max + step, step)

        l = []

        for i in range(max_col):
            col = self.__make_col(*self.f_pal(i))
            l.append((values[i], col))

        legend = Legend(l)

        inset = biggles.Inset((1.1, 0.60), (1.2, .97), legend)

        return inset

    def __normalize(self, matrix):

        _max = max(N.ravel(matrix))
        _min = min(N.ravel(matrix))
        
        s = N.shape(matrix)

        max_col = self.max_col[self.palette]

        matrix = (matrix - _min) * max_col / (_max - _min)

        return _min, _max, matrix

    def __map(self, x):
        return int(min(255, round(x)))

    def __make_col(self, red, green, blue):
        return ((red << 16) + (green << 8) + blue)        

    def __grey(self, x):
        x = self.__map(255 * x / self.max_col['grey'])
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
    
if __name__ == '__main__':

    from Numeric import *

    n = 30

    z = zeros((n,n), Float)

    for i in range(shape(z)[0]):
        for j in range(shape(z)[1]):
            z[i,j] = exp(-0.01*((i-n/2)**2+(j-n/2)**2))

    p = MatrixPlot(z, palette='sausage', legend = 1)
    p.show()
