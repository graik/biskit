## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004 Michael Habeck
## Copyright (C) 2005 Raik Gruenberg
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

## Contributions: Michael Habeck, Raik Gruenberg
## $Revision$
## last $Date$
## last $Author$

"""
lognormal distribution
"""

import numpy as N


def rand_log_normal(alpha, beta, shape):
    return N.exp(N.random.normal(alpha, beta, shape))


def ln(r, alpha, beta):
    return N.exp(-0.5/beta**2 * (N.log(r) - alpha)**2 \
                 - 0.5*N.log(2*N.pi)-N.log(beta*r))


def erf(x):
    """
    Approximation to the erf-function with fractional error
    everywhere less than 1.2e-7

    @param x: value
    @type  x: float

    @return: value
    @rtype: float
    """
    if x > 10.: return 1.
    if x < -10.: return -1.

    z = abs(x)
    t = 1. / (1. + 0.5 * z)

    r = t * N.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + \
                                                               t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * \
                                                                                                         (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * \
                                                                                                                                               0.17087277)))))))))

    if x >= 0.:
        return 1. - r
    else:
        return r - 1.


def logArea(x, alpha, beta):
    """
    Area of the smallest interval of a lognormal distribution that still
    includes x.

    @param x: border value
    @type  x: float
    @param alpha: mean of log-transformed distribution
    @type  alpha: float
    @param beta: standarddev of log-transformed distribution
    @type  beta: float

    @return: probability that x is NOT drawn from the given distribution
    @rtype: float
    """
    r_max = N.exp(alpha - beta**2)

    if x < r_max: x = r_max**2 / x

    upper = (N.log(x) - alpha) / beta 

    return 0.5 * (erf(upper / N.sqrt(2)) - erf(-(upper + 2*beta) / N.sqrt(2)))


def logMean( alpha, beta ):
    """
    @param alpha: mean of log-transformed distribution
    @type  alpha: float
    @param beta: standarddev of log-transformed distribution
    @type  beta: float

    @return: mean of the original lognormal distribution
    @rtype: float
    """
    return N.exp( alpha + (beta**2)/2. )


def logSigma( alpha, beta ):
    """
    @param alpha: mean of log-transformed distribution
    @type  alpha: float
    @param beta: standarddev of log-transformed distribution
    @type  beta: float

    @return: 'standard deviation' of the original lognormal distribution
    @rtype: float
    """
    return logMean( alpha, beta ) * N.sqrt( N.exp(beta**2) - 1.)


def logMedian( alpha, beta=None ):
    """
    @param alpha: mean of log-transformed distribution
    @type  alpha: float
    @param beta: not needed
    @type  beta: float

    @return: median of the original lognormal distribution
    @rtype: float
    """
    return N.exp( alpha )


def logConfidence( x, R, clip=0 ):
    """
    Estimate the probability of x NOT beeing a random observation from a
    lognormal distribution that is described by a set of random values.

    @param x: observed value
    @type  x: float
    @param R: sample of random values
    @type  R: [float]
    @param clip: clip zeros at this value  0->don't clip (default: 0)
    @type  clip: float

    @return: confidence that x is not random, median of random distr.
    @rtype: (float, float)
    """
    if clip and 0 in R:
        R = N.clip( R, clip, max( R ) )
    if clip and x == 0:
        x = clip

    ## remove 0 instead of clipping
    R = N.compress( R, R )
    if x == 0:
        return 0, 0

    ## get mean and stdv of log-transformed random sample
    alpha = N.average( N.log( R ) )

    n = len( R )

    beta = N.sqrt(N.sum(N.power(N.log( R ) - alpha, 2)) / (n - 1.))

    return logArea( x, alpha, beta ), logMedian( alpha )



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """

    def test_lognormal(self):
        """Statistics.lognormal test"""
        import random
        import Biskit.gnuplot as gnuplot
        import Biskit.hist as H

        cr = []
        for i in range( 10000 ):
            ## Some random values drawn from the same lognormal distribution 

            alpha = 1.5
            beta = .7
            x = 10.

            R = [ random.lognormvariate( alpha, beta ) for j in range( 10 ) ]

            cr += [ logConfidence( x, R )[0] ]


        ca = logArea( x, alpha, beta )

        if self.local:
            gnuplot.plot( H.density( N.array(cr) - ca, 100 ) )

            globals().update( locals() )

        self.assertAlmostEqual( ca,  0.86877651432955771, 7)


if __name__ == '__main__':

    BT.localTest()
