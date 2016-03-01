## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2002-2005 Michael Habeck
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

## Contributions: Michael Habeck, Raik Gruenberg, Johan Leckner
## last $Author$
## last $Date$
## $Revision$
"""
superimpose 2 structures iteratively
"""

import Biskit.mathUtils as MU
import numpy as N
from numpy.linalg import svd

## def average(x):
##     return N.sum(x) / len(x)

## def variance(x):
##     return N.average(N.power(x - N.average(x), 2))

## def standardDev(x):
##     return N.sqrt(variance(x))


def findTransformation(x, y):
    """
    Match two arrays by rotation and translation. Returns the
    rotation matrix and the translation vector.

    @param x: first set of coordinates
    @type  x: array('f')
    @param y: second set of coordinates
    @type  y: array('f')

    @return: rotation matrix (3x3) and translation vector (1x3)
    @rtype:  array, array
    """
    ## center configurations
    x_av = N.average(x)
    y_av = N.average(y)

    x = x - x_av
    y = y - y_av


    ## svd of correlation matrix
    v, l, u = svd(N.dot(N.transpose(x), y))

    ## build rotation matrix and translation vector
    r = N.dot(v, u)

    t = x_av - N.dot(r, y_av)

    return r, t


def match(x, y, n_iterations=1, z=2, eps_rmsd=0.5, eps_stdv=0.05):
    """
    Matches two arrays onto each other, while iteratively removing outliers.
    Superimposed array y would be C{ N.dot(y, N.transpose(r)) + t }.

    @param n_iterations: number of calculations::
                           1 .. no iteration 
                           0 .. until convergence
    @type  n_iterations: 1|0
    @param z: number of standard deviations for outlier definition (default: 2)
    @type  z: float
    @param eps_rmsd: tolerance in rmsd (default: 0.5)
    @type  eps_rmsd: float
    @param eps_stdv: tolerance in standard deviations (default: 0.05)
    @type  eps_stdv: float

    @return: (r,t), [ [percent_considered, rmsd_for_it, outliers] ]
    @rtype: (array, array), [float, float, int]
    """
    iter_trace = []

    rmsd_old = 0
    stdv_old = 0

    n = 0
    converged = 0

    mask = N.ones(len(y), N.int32 )

    while not converged:

        ## find transformation for best match
        r, t = findTransformation(N.compress(mask, x, 0),
                                  N.compress(mask, y, 0))

        ## transform coordinates
        xt = N.dot(y, N.transpose(r)) + t

        ## calculate row distances
        d = N.sqrt(N.sum(N.power(x - xt, 2), 1)) * mask

        ## calculate rmsd and stdv
        rmsd = N.sqrt(N.average(N.compress(mask, d)**2))
        stdv = MU.SD(N.compress(mask, d))

        ## check conditions for convergence
        d_rmsd = abs(rmsd - rmsd_old)
        d_stdv = abs(1 - stdv_old / stdv)

        if d_rmsd < eps_rmsd and d_stdv < eps_stdv:
            converged = 1
        else:
            rmsd_old = rmsd
            stdv_old = stdv

        ## store result
        perc = round(float(N.sum(mask)) / float(len(mask)), 2)

        ## throw out non-matching rows
        mask = N.logical_and(mask, N.less(d, rmsd + z * stdv))
        outliers = N.nonzero( N.logical_not( mask ) )[0]
        iter_trace.append([perc, round(rmsd, 3), outliers])

        n += 1

        if n_iterations and n >= n_iterations:
            break

    return (r, t), iter_trace


def rowDistances( x, y ):
    """
    Calculate the distances between the items of two arrays (of same shape)
    after least-squares superpositioning.

    @param x: first set of coordinates
    @type  x: array('f')
    @param y: second set of coordinates
    @type  y: array('f')  

    @return: array( len(x), 'f' ), distance between x[i] and y[i] for all i
    @rtype: array
    """
    ## find transformation for best match
    r, t = findTransformation(x, y)

    ## transform coordinates
    z = N.dot(y, N.transpose(r)) + t

    ## calculate row distances
    return N.sqrt(N.sum(N.power(x - z, 2), 1)) 



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_rmsFit( self ):
        """rmsFit test"""
        import Biskit.tools as T

        self.traj = T.load( T.testRoot() + '/lig_pcr_00/traj.dat' )

        rt, rmsdLst = match( self.traj.ref.xyz, self.traj[-1].xyz)

        if self.local:
            print 'RMSD: %.2f' % rmsdLst[0][1]

        # return rotation matrix
        r = abs( N.sum( N.ravel( rt[0] )))
        e = abs( N.sum( N.ravel( self.EXPECT )))

        self.assertAlmostEqual(r, e, 6)

    EXPECT = N.array( [[ 0.9999011,   0.01311352,  0.00508244,],
                       [-0.01310219,  0.99991162, -0.00225578,],
                       [-0.00511157,  0.00218896,  0.99998454 ]] )


if __name__ == '__main__':

    BT.localTest()





