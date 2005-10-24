## 
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (c) 2004 Wolfgang Rieping & Michael Habeck
## Copyright (c) 2005 Raik Gruenberg
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

## Contributions: Wolfgang Rieping, Michael Habeck, Raik Gruenberg
## last $Author$
## last $Date$
## $Revision$

import Numeric as N
import math
import Biskit.hist as H

class Density:
    """
    Analyze a density distribution of values.
    Can be created from a list of samples or from a discrete distribution:
    Density( values = [ float ] )
    or
    Density( p = array(2xN,'f') )
    by W.+M.
    """

    def __init__(self, p = None, values = None, bins = 20):

        if p is None and values is None:
            raise ValueError, 'Either a discrete distribution or ' + \
                  'a list of samples must be specified'

        if values is not None:

            p = H.density(values, bins, steps = 0)

        self.store(p)

        self.verbose = 0

    def store(self, val):
        
        self.val = N.array(val, N.Float)
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
        level - float, confidence level (e.g. 0.68 for stdev interval)
        -> ( float, float ), start and end of the confidence interval
                             containing |level|*100 % of the probability 
        """
        order = N.argsort(self.p).tolist()
        cumulative = N.add.accumulate(N.take(self.p, order)) * self.delta_x

        ind = N.nonzero(N.greater_equal(cumulative, 1. - level))
        
        sub_set = order[ind[0]:]

        intervals = self.__find_intervals(sub_set)

        boundaries = [(self.x[i[0]], self.x[i[-1]]) for i in intervals]

        return tuple(boundaries)

    def findConfidenceInterval(self, x):
        """
        findConfidenceInterval(self, x)
        Find the smallest possible density interval that still includes x.
        x - float, value
        -> ( float, (float, float) ), convidence level, interval start and end 
        """
        closest = N.argmin(abs(self.x - x))

        ind = N.nonzero(N.greater_equal(self.p, self.p[closest])).tolist()

        intervals = self.__find_intervals(ind)

        lens = N.array([len(i) for i in intervals])
        levels = [N.sum(N.take(self.p, i)) for i in intervals]
        level = N.sum(levels) * self.delta_x

        boundaries = [(self.x[i[0]], self.x[i[-1]]) for i in intervals]

        return level, tuple(boundaries)

    def median(self):
        cum = N.add.accumulate(self.p) * self.delta_x
        index = N.argmin(abs(cum - 0.5))

        return self.x[index]

    def average(self):
        return self.delta_x * N.sum(self.p * self.x)

    def max(self):
        index = N.argmax(self.p)
        return self.x[index]


    def __find_intervals(self, l):

        l = N.array(l)
        l = N.take(l, N.argsort(l))

        break_points = N.nonzero(N.greater(l[1:] - l[:-1], 1))

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
    The two parameters can hence be calculated from n sample values v:
    alpha = 1/n N.sum( ln(vi) )
    beta = N.sqrt( 1/(n-1) N.sum( ln(vi) - alpha )^2  )
    
    x     - float, value
    alpha - float, mean of the log-transformed random variable
    beta  - float, stdev of the log-transformed random variable
    -> float, probability of x
    """
    return 1. / math.sqrt(2. * math.pi) / beta / x * \
           math.exp(-0.5 / beta ** 2 * (math.log(x) - alpha) ** 2)


def logConfidence( x, R, clip=1e-32 ):
    """
    Estimate the probability of x NOT beeing a random observation from a
    lognormal distribution that is described by a set of random values.
    The exact solution to this problem is in Biskit/Statistics/lognormal.py.
    x    - float, observed value
    R    - [ float ], sample of random values 0->don't clip (1e-32)
    clip - float, clip zeros at this value
    -> (float, float) confidence that x is not random, mean of random distrib.
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


#########
## TEST

if __name__ == '__main__':

    import random

    ## a proper lognormal density distribution the log of which has mean 1.0
    ## and stdev 0.5
##  X = [ (x, p_lognormal(x, 1.0, 0.5)) for x in N.arange(0.00001, 50, 0.001)]

    for i in range( 10 ):
        ## Some random values drawn from the same lognormal distribution 

        alpha = 2.
        beta = 0.6

        R = [ random.lognormvariate( alpha, beta ) for i in range( 10000 ) ]
        
        print logConfidence( 6.0, R )[0], area(6.0, alpha, beta) 
        
