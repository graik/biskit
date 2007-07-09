##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (c) 2004 Michael Habeck
## Copyright (c) 2005 Raik Gruenberg
## All rights reserved
##
## Contribution: Michael Habeck, Raik Gruenberg
##
## $Revision$
## last $Author$

import Numeric as N
import RandomArray as R


def rand_log_normal(alpha, beta, shape):

    return N.exp(R.normal(alpha, beta, shape))


def ln(r, alpha, beta):

    return N.exp(-0.5/beta**2 * (N.log(r) - alpha)**2 \
                 - 0.5*N.log(2*pi)-N.log(beta*r))


def erf(x):
    """
    Approximation to the erf-function with fractional error
    everywhere less than 1.2e-7
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
    x     - float, border value
    alpha - float, mean of log-transformed distribution
    beta  - float, standarddev of log-transformed distribution
    -> float, probability that x is NOT drawn from the given distribution
    """
    r_max = N.exp(alpha - beta**2)

    if x < r_max: x = r_max**2 / x

    upper = (N.log(x) - alpha) / beta 

    return 0.5 * (erf(upper / N.sqrt(2)) - erf(-(upper + 2*beta) / N.sqrt(2)))


def logMean( alpha, beta ):
    """
    alpha - float, mean of log-transformed distribution
    beta  - float, standarddev of log-transformed distribution
    -> float, mean of the original lognormal distribution
    """
    return N.exp( alpha + (beta**2)/2. )
    

def logSigma( alpha, beta ):
    """
    alpha - float, mean of log-transformed distribution
    beta  - float, standarddev of log-transformed distribution
    -> float, 'standard deviation' of the original lognormal distribution
    """
    return logMean( alpha, beta ) * N.sqrt( N.exp(beta**2) - 1.)


def logMedian( alpha, beta=None ):
    """
    alpha - float, mean of log-transformed distribution
    beta  - float, not needed
    -> float, median of the original lognormal distribution
    """
    return N.exp( alpha )


def logConfidence( x, R, clip=0 ):
    """
    Estimate the probability of x NOT beeing a random observation from a
    lognormal distribution that is described by a set of random values.
    x    - float, observed value
    R    - [ float ], sample of random values
    clip - float, clip zeros at this value  0->don't clip [0]
    -> (float, float) confidence that x is not random, median of random distr.
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


if __name__ == '__main__':

    import random

    cr = []
    for i in range( 10000 ):
        ## Some random values drawn from the same lognormal distribution 

        alpha = 1.5
        beta = .7
        x = 10.

        R = [ random.lognormvariate( alpha, beta ) for i in range( 10 ) ]

        cr += [ logConfidence( x, R )[0] ]


    ca = logArea( x, alpha, beta )

    plot( density( N.array(cr) - ca, 100 ) )