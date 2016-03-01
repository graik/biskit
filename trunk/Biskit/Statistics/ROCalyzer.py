## 
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (c) 2007 Raik Gruenberg
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

## Contributions: Benjamin Bouvier
## last $Author$
## last $Date$
## $Revision$
"""Generate and analyze ROC curves of specificity versus sensitivity."""
import numpy as N
import copy, random

import lognormal as L
try:
    import scipy.stats as stats
    USING_SCIPY = True
except:
    # import outdated stand-alone internal copy with some numeric issues
    import stats
    USING_SCIPY = False

class ROCError( Exception ):
    pass

class ROCalyzer( object ):
    """
    Generate and analyze ROC (reciever operating characteristic) curves
    which describe how successful a predictor score is at identifying a list
    of elements (positives) within a larger sequence (of negatives).

    A ROC curve plots the sensitivity of the score (fraction of true positives
    captured) against its specificity (fraction of true negatives identified
    as such). An ideal score would hit the point 1 : 1 for some combination of
    specificity and sensitivity -- that means it would discover all true
    positives without picking any false positives. A random score simply
    follows the diagonal in the plot -- that means it captures true positives
    at the same rate as false positives. 

    We stick to the definition used here: [quote Lichtarge JMB 2007]

    Usage::

      a = ROCalyzer( mask_labelling_positives )

      curve = a.roccurve( score )
      Biskit.gnuplot.plot( curve )     ## visualize the ROC curve

      success = a.rocarea( score )     ## the area under the roccurve - 0.5
      justrandom = a.isnoise( score )  ## probability that score is just random
    """

    def __init__( self, positives ):
        """
        @param positives: mask identifying all true positives
        @type  positives: [ 1|0 ]
        """
        self.positives = positives


    def random_roccurves(self, score, n=100 ):
        """
        Generate roc curves for randomized target masks.
        @param score: sequence of score values for target sequence
        @type  score: [ int ] or [ float ]
        @param n: number of curves to generate (100)
        @type  n: int
        @return: a set of sensitivity/specifity curves for the given score
        applied to random 'positive' items.
        @rtype: [ (sens, spec), ]
        """
        N.random.seed()
        r = []
        for i in range(n):
            ref = copy.copy( self.positives )
            N.random.shuffle( ref )

            r += [ self.roccurve( score, ref=ref ) ]

        return r


    def roccurve( self, score, ref=None ):
        """
        Calculate the ROC curve of the given score.
        @todo: doesn't handle target containing only true positives or
        only true negatives

        @param score: sequence of score values for target sequence
        @type  score: [ int ] or [ float ]
        @param ref  : alternative mask of positives (overrides self.positives)
        @type  ref  : [ 1|0 ]
        @return: curve describing 1-specificity (1st column) versus sensitivity
        @rtype: [ (specifity, false_positive_rate) ]
        """
        if ref is None:
            ref = self.positives

        if N.sum( ref ) == len( ref ) or N.sum(ref) == 0:
            raise ROCError,\
                  'Cannot compute Roc curves for all positive or all '+\
                  'negative target'
        order = N.argsort( score ).tolist()
        order.reverse()

        ref = N.take( ref, order, 0 )

        #: number of true positives identified with decreasing score
        n_pos = N.add.accumulate( ref )

        # number of false positives picked up with decreasing score
        neg = N.logical_not( ref )
        n_neg = N.add.accumulate( neg  )

        sensitivity = 1. * n_pos / n_pos[-1] ## true positive rate
        specificity = 1. * n_neg / n_neg[-1] ## FP_rate = 1-specificity

        return zip( specificity, sensitivity )


    def area( self, curve, start=0.0, stop=1.0 ):
        """
        Numerically add up the area under the given curve.
        The curve is a 2-D array or list of tupples as returned by roccurve().
        The x-axis is the first column of this array (curve[:,0]).

        @param curve: a list of x,y coordinates
        @type  curve: [ (y,x), ] or N.array
        @param start: lower boundary (in x) (default: 0.0)
        @type  start: float
        @param stop: upper boundary (in x) (default: 1.0)
        @type  stop: float
        @return: the area underneath the curve between start and stop.
        @rtype: float
        """
        c = N.array( curve )
        assert len( N.shape( c ) ) == 2

        ## apply boundaries  ## here we have a problem with flat curves
        mask = N.greater_equal( c[:,0], start )
        mask *= N.less( c[:,0], stop+1e-15 )
        c = N.compress( mask, c, axis=0 )

        ## fill to boundaries -- not absolutely accurate: we actually should
        ## interpolate to the neighboring points instead
        if c[0,0] > start:
            c = N.concatenate( (N.array([[start, c[0,1]],]), c) )
        if c[-1,0] < stop:
            c = N.concatenate( (c, N.array([[stop, c[-1,1]],])) )

        x = c[:,0]
        y = c[:,1]

        dx = x[1:] - x[:-1] # distance on x between points 
        dy = y[1:] - y[:-1] # distance on y between points

        areas1 = y[:-1] * dx  # the rectangles between all points
        areas2 = dx * dy / 2.0 # the triangles between all points

        return N.sum(areas1) + N.sum(areas2)


    def rocarea( self, score ):
        """
        Calculate how much the score outperforms a perfect random score.
        The measure used is the area underneath the ROC curve. A perfect
        random score should, on average, form a diagonal roc curve with an
        area of 0.5. Anything above that *may* indicate better-than random
        performance. See also isnoise() for quantifying this.
        @param score: the score predicted for each item
        @type  score: [ float ]
        @return: a.area( a.roccurve( score) )
        @rtype: float
        """
        return self.area( self.roccurve( score ) )


    def isnoise( self, score, n_samples=1000 ):
        """
        Test sample how a given score performs at predicting items in the
        positive list compared to its 'performance' at  predicting random
        elements. The result corresponds to a two-tailed P value.
        See L{utest} for the analytical solution.
        @param score: the score predicted for each item
        @type  score: [ float ]
        @param n_samples: number of random samples
        @type  n_samples: int

        @return: probability P that the prediction success of score is just
        a random effect (1.0 means it's just perfectly random).
        """
        from Biskit import EHandler

        ## list of random deviations from diagonal area 0.5
        a_rand = [ self.area(c)-0.5
                   for c in self.random_roccurves(score,n_samples) ]

        sd_rand = N.std( a_rand )
        av_rand = N.mean(a_rand )

        if round(av_rand,2) != 0.0:
            EHandler.warning( 'random sampling is skewed by %f'% (av_rand-0.0))

        a = self.rocarea( score )
        z = a / sd_rand

        ## probability that a sample falls *within* z stdevs from the mean
        p = L.erf( z / N.sqrt(2) )

        ## probability that the score hits just at random
        return 1.0 - p


    def utest( self, score ):
        """
        Gives the Mann-Withney U test probability that the score is
        random.  See:

        Mason & Graham (2002) Areas beneath the relative operating
        characteristics (ROC) and relative operating levels (ROL)
        curves: Statistical significance and interpretation

        Note (1): P-values below ~1e-16 are reported as 0.0.
        See zprob() in Biskit.Statistics.stats!

        Note (2): the P-value does not distinguish between positive
        and negative deviations from random -- a ROC area of 0.1 will
        get the same P-value as a ROC area of 0.9.

        @param score: the score predicted for each item
        @type  score: [ float ]

        @return: 1-tailed P-value
        @rtype: float
        """
        sample1 = N.compress( self.positives, score )
        sample1 = sample1[-1::-1]  # invert order

        sample2 = N.compress( N.logical_not( self.positives ), score )
        sample2 = sample2[-1::-1]  # invert order

        sample1 = sample1.tolist()
        sample2 = sample2.tolist()

        p = stats.mannwhitneyu( sample1, sample2 )
        return p[1]



class ROCThreshold(object):

    def __init__( self, target ):
        """
        @param target: the target profile to be predicted by a score
        @type  target: [ float ]
        """
        self.target = target

    def target2mask( self, n=1 ):
        """
        Select n highest points from target profile as True Positives.
        @param n: number of true positives to be selected
        @type  n: int
        """
        order = N.argsort( self.target ).tolist()
        order.reverse()

        r = N.zeros( len(self.target), int )
        N.put( r, order[:n], 1 )
        return r

    def threshold_curve( self, score ):
        """
        @param target: the target profile to be predicted
        @type  target: [ float ]
        @param score: the score predicted for each item
        @type  score: [ float ]
        """
        import Biskit.tools as T
        ## order from highest to lowest
        ordered = N.take( self.target, N.argsort( self.target ), 0 )
        ordered = ordered[::-1]

        r = N.zeros( (len(self.target), 2) )
        r[:,0] = ordered

        for n in range(len(self.target)-1):
            roc = ROCalyzer( self.target2mask( n+1 ) )

            r[n,0] = ordered[n]
            r[n,1] = roc.rocarea( score )

        return r

def pfisher( pvalues ):
    """
    Combine independent P-values into one according to

    Fisher, R. A. (1948) Combining independent tests of significance.
    American Statistician, vol. 2, issue 5, page 30.

    ('Fisher method' or 'inverse ChiSquare method') See also book:
    Walter W. Piegorsch, A. John Bailer: Analyzing Environmental Data.
    Wiley 2005

    @param pvalues: list of independent P-values
    @type  pvalues: [ float ]
    @return: P-value 
    @rtype: float
    """
    ## stats.mannwhitneyu minimal P ~ stats.zprob( 8.2 );
    ## all below becomes 0. which is not handled by the fisher test 
    clipped = N.clip( pvalues, 1.0e-16, 1.0 )

    x2 = -2 * N.sum( N.log( clipped ) )

    if not USING_SCIPY:
        x2 = float( x2 )

    return stats.chisqprob( x2, 2*len(pvalues) )


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """
    def prepare( self ):
        import Biskit.tools as T

        self.cl = T.load( T.testRoot()+'/dock/hex/complexes.cl')

        self.score = self.cl.valuesOf('hex_eshape')
## ## convert hex energies into positive score
##  self.score = N.array(self.cl.valuesOf('hex_etotal')) * -1

        ## define complexes with less the 6 A rmsd from reference as positives
        self.hits = N.less( self.cl.valuesOf('rms'), 6 )

    def test_roccurve(self):
        """Statistics.ROCalyzer test"""
        from Biskit.gnuplot import plot

        self.a = ROCalyzer( self.hits )
        self.roc = self.a.roccurve( self.score )

        if self.local:
            plot( self.roc )

    def test_noise(self):
        """Statistics.ROCalyzer.isnoise/utest test"""
        a = ROCalyzer( self.hits )

        rand = a.random_roccurves( self.score, 500 )
        self.avg = N.average( N.array( rand ), axis=0 )

        self.assertAlmostEqual( a.area(self.avg), 0.5, 1 )

        p1 = a.isnoise( self.score, n_samples=500 )
        p2 = a.utest( self.score )

        if self.local:
            print "\nsimulation  P: ", p1
            print "mannwithney P: ", p2

        self.assertAlmostEqual( p1, 0.0, 3 )
        self.assertAlmostEqual( p2, 0.0, 3 )

##         r =  p1 / p2
##         self.assert_( (r<3.25) and (r>0.75),
##                       'isnoise P should be about twice that of utest'+\
##                       ' isnoise : utest = %f : %f' % (p1, p2))


    def test_area(self):
        """Statistics.ROCalyzer.area test"""
        a = ROCalyzer( self.hits )

        perfect_y = N.arange(1.0, 0,   -1.0/len(self.score) )
        perfect_x = N.arange(0.0, 1.0, +1.0/len(self.score) )
        self.perfect = zip( perfect_x, perfect_y )

        self.assertAlmostEqual( a.area( self.perfect ), 0.5, 5 )

    def test_threshold(self):
        """Statistics.ROCThreshold test"""
        from Biskit.gnuplot import plot

        target = 1./ N.array( self.cl.valuesOf('rms') )

        rt = ROCThreshold( target )
        positives = rt.target2mask( N.sum( self.hits ) )

        self.t_curve = rt.threshold_curve( self.score )

        if self.local:
            plot( self.t_curve )

        self.assert_( max(self.t_curve[:,1]) > 0.4 )

    def test_pfisher( self ):
        """Statistics.ROCalyzer.pfisher test"""
        p = pfisher( [0.046, .072, .0005, .999, 0.198, .944, .0015, 0.0075,
                      .0625, .003, .0014, .0005, .238, .988])
        self.assertAlmostEqual( p, 3.2675687582319262e-10, 5 )


if __name__ == '__main__':

    BT.localTest()

