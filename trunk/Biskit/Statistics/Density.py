## 
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (c) 2004 Wolfgang Rieping & Michael Habeck
## Copyright (c) 2005 Raik Gruenberg
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

## Contributions: Wolfgang Rieping, Michael Habeck, Raik Gruenberg
## last $Author$
## last $Date$
## $Revision$

"""
Analyze a density distribution of values.
"""

import numpy as N
import math
import Biskit.hist as H

class Density:
    """
    Analyze a density distribution of values.
    Can be created from a list of samples or from a discrete distribution::

      Density( values = [ float ] )
        or
      Density( p = array(2xN,'f') )
    """

    def __init__(self, p = None, values = None, bins = 20):
        """
        @param p: discrete distribution, array (2 x len(data) ) with
                  start of bin and witdh of bin (default: None)
        @type  p: array
        @param values: list of samples  (default: None)
        @type  values: [float]
        @param bins: number of bins  (default: 20)
        @type  bins: int
        """
        if p is None and values is None:
            raise ValueError, 'Either a discrete distribution or ' + \
                  'a list of samples must be specified'

        if values is not None:

            p = H.density(values, bins, steps = 0)

        self.store(p)

        self.verbose = 0


    def store(self, val):
        """
        Analyze distribution data.

        @param val: array (2 x len(data) ) with start of bin and witdh
                    of bin (default: None)
        @type  val: array        
        """
        self.val = N.array(val, N.float32)
        self.x = self.val[:,0]
        self.p = self.val[:,1]

        self.delta_x = abs(self.x[0] - self.x[1])

        Z = N.sum(self.p) * self.delta_x

        self.p /= Z


    def get(self):
        return self.val


    def __getslice__(self, *args, **kw):
        from operator import getslice
        return getslice(self.val, *args, **kw)


    def __getitem__(self, *args, **kw):
        from operator import getitem
        return getitem(self.val, *args, **kw)


    def confidenceInterval(self, level):
        """
        confidenceInterval(self, level)

        @param level: confidence level (e.g. 0.68 for stdev interval)
        @type  level: float

        @return: start and end of the confidence interval
                 containing |level|*100 % of the probability
        @rtype: float, float
        """          
        order = N.argsort(self.p).tolist()
        cumulative = N.cumsum(N.take(self.p, order)) * self.delta_x

        ind = N.nonzero(N.greater_equal(cumulative, 1. - level))[0]

        sub_set = order[ind[0]:]

        intervals = self.__find_intervals(sub_set)

        boundaries = [(self.x[i[0]], self.x[i[-1]]) for i in intervals]

        return tuple(boundaries)


    def findConfidenceInterval(self, x):
        """
        findConfidenceInterval(self, x)
        Find the smallest possible density interval that still includes x.

        @param x: value
        @type  x: float

        @return: convidence level, interval start and end
        @rtype: float, (float,float)
        """
        closest = N.argmin(abs(self.x - x))

        ind = N.nonzero(N.greater_equal(self.p, self.p[closest]))[0].tolist()

        intervals = self.__find_intervals(ind)

##        lens = N.array([len(i) for i in intervals])
        levels = [N.sum(N.take(self.p, i)) for i in intervals]
        level = N.sum(levels) * self.delta_x

        boundaries = [(self.x[i[0]], self.x[i[-1]]) for i in intervals]

        return level, tuple(boundaries)


    def median(self):
        """
        Median of distribution.
        """
        cum = N.cumsum(self.p) * self.delta_x
        index = N.argmin(abs(cum - 0.5))

        return self.x[index]


    def average(self):
        """
        Average of distribution.
        """
        return self.delta_x * N.sum(self.p * self.x)


    def max(self):
        """
        Max height of distribution.
        """        
        index = N.argmax(self.p)
        return self.x[index]


    def __find_intervals(self, l):
        l = N.array(l)
        l = N.take(l, N.argsort(l))

        globals().update( locals() )

        break_points = N.nonzero(N.greater(l[1:] - l[:-1], 1))[0]

        start = 0
        intervals = []

        for i in range(len(break_points)):
            index = break_points[i]
            intervals.append(tuple(N.take(l, range(start, index + 1))))
            start = index + 1

        intervals.append(tuple(l[start:]))

        return intervals


def p_lognormal(x, alpha, beta):
    """
    Get the probability of x from a lognormal density distribution that
    is described by two parameters alpha and beta. Alpha and beta are not
    the usual mean and standard dev of the lognormal distribution itself but
    are the mean and stdev of the distribution after log-transformation.
    The two parameters can hence be calculated from n sample values v::

      alpha = 1/n N.sum( ln(vi) )
      beta = N.sqrt( 1/(n-1) N.sum( ln(vi) - alpha )^2  )

    @param x: value
    @type  x: float
    @param alpha: mean of the log-transformed random variable
    @type  alpha: float
    @param beta: stdev of the log-transformed random variable
    @type  beta: float

    @return: probability of x
    @rtype: float
    """
    return 1. / math.sqrt(2. * math.pi) / beta / x * \
           math.exp(-0.5 / beta ** 2 * (math.log(x) - alpha) ** 2)


def logConfidence( x, R, clip=1e-32 ):
    """
    Estimate the probability of x NOT beeing a random observation from a
    lognormal distribution that is described by a set of random values.
    The exact solution to this problem is in L{Biskit.Statistics.lognormal}.

    @param x: observed value
    @type  x: float
    @param R: sample of random values; 0 -> don't clip (default: 1e-32)
    @type  R: [float]
    @param clip: clip zeros at this value
    @type  clip: float

    @return:  confidence that x is not random, mean of random distrib.
    @rtype: (float, float)
    """
    if clip and 0 in R:
        R = N.clip( R, clip, max( R ) )
    ## get mean and stdv of log-transformed random sample
    mean = N.average( N.log( R ) )

    n = len( R )

    stdv = N.sqrt(N.sum(N.power(N.log( R ) - mean, 2)) / (n - 1.))

    ## create dense lognormal distribution representing the random sample
    stop = max( R ) * 50.0
    step = stop / 100000
    start = step / 10.0

    X = [(v, p_lognormal(v, mean, stdv) ) for v in N.arange(start, stop, step)]

    ## analyse distribution
    d = Density( X )

    return d.findConfidenceInterval( x * 1.0 )[0], d.average()


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """

    def test_Density(self):
        """Statistics.Density test"""
        import random

        ## a lognormal density distribution the log of which has mean 1.0
        ## and stdev 0.5
        self.X = [ (x, p_lognormal(x, 1.0, 0.5))
                   for x in N.arange(0.00001, 50, 0.001)]

        alpha = 2.
        beta = 0.6

        self.R = [ random.lognormvariate( alpha, beta )
                   for i in range( 10000 )]

        p = logConfidence( 6.0, self.R )[0]#, area(6.0, alpha, beta)


if __name__ == '__main__':

    BT.localTest()
