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
"""
superimpose 2 structures iteratively
"""

from . import mathUtils as MU
from biskit.core import oldnumeric as N0
from numpy.linalg import svd

## def average(x):
##     return N0.sum(x) / len(x)

## def variance(x):
##     return N0.average(N0.power(x - N0.average(x), 2))

## def standardDev(x):
##     return N0.sqrt(variance(x))


def findTransformation(x, y):
    """
    Match two arrays by rotation and translation. Returns the
    rotation matrix and the translation vector.

    :param x: first set of coordinates
    :type  x: array('f')
    :param y: second set of coordinates
    :type  y: array('f')

    :return: rotation matrix (3x3) and translation vector (1x3)
    :rtype:  array, array
    """
    ## center configurations
    x_av = N0.average(x)
    y_av = N0.average(y)

    x = x - x_av
    y = y - y_av


    ## svd of correlation matrix
    v, l, u = svd(N0.dot(N0.transpose(x), y))

    ## build rotation matrix and translation vector
    r = N0.dot(v, u)

    t = x_av - N0.dot(r, y_av)

    return r, t


def match(x, y, n_iterations=1, z=2, eps_rmsd=0.5, eps_stdv=0.05):
    """
    Matches two arrays onto each other, while iteratively removing outliers.
    Superimposed array y would be C{ N0.dot(y, N0.transpose(r)) + t }.

    :param n_iterations: number of calculations::
                           1 .. no iteration 
                           0 .. until convergence
    :type  n_iterations: 1|0
    :param z: number of standard deviations for outlier definition (default: 2)
    :type  z: float
    :param eps_rmsd: tolerance in rmsd (default: 0.5)
    :type  eps_rmsd: float
    :param eps_stdv: tolerance in standard deviations (default: 0.05)
    :type  eps_stdv: float

    :return: (r,t), [ [percent_considered, rmsd_for_it, outliers] ]
    :rtype: (array, array), [float, float, int]
    """
    iter_trace = []

    rmsd_old = 0
    stdv_old = 0

    n = 0
    converged = 0

    mask = N0.ones(len(y), N0.Int32 )

    while not converged:

        ## find transformation for best match
        r, t = findTransformation(N0.compress(mask, x, 0),
                                  N0.compress(mask, y, 0))

        ## transform coordinates
        xt = N0.dot(y, N0.transpose(r)) + t

        ## calculate row distances
        d = N0.sqrt(N0.sum(N0.power(x - xt, 2), 1)) * mask

        ## calculate rmsd and stdv
        rmsd = N0.sqrt(N0.average(N0.compress(mask, d)**2))
        stdv = MU.SD(N0.compress(mask, d))

        ## check conditions for convergence
        d_rmsd = abs(rmsd - rmsd_old)
        d_stdv = abs(1 - stdv_old / stdv)

        if d_rmsd < eps_rmsd and d_stdv < eps_stdv:
            converged = 1
        else:
            rmsd_old = rmsd
            stdv_old = stdv

        ## store result
        perc = round(float(N0.sum(mask)) / float(len(mask)), 2)

        ## throw out non-matching rows
        mask = N0.logical_and(mask, N0.less(d, rmsd + z * stdv))
        outliers = N0.nonzero( N0.logical_not( mask ) )
        iter_trace.append([perc, round(rmsd, 3), outliers])

        n += 1

        if n_iterations and n >= n_iterations:
            break

    return (r, t), iter_trace


def rowDistances( x, y ):
    """
    Calculate the distances between the items of two arrays (of same shape)
    after least-squares superpositioning.

    :param x: first set of coordinates
    :type  x: array('f')
    :param y: second set of coordinates
    :type  y: array('f')  

    :return: array( len(x), 'f' ), distance between x[i] and y[i] for all i
    :rtype: array
    """
    ## find transformation for best match
    r, t = findTransformation(x, y)

    ## transform coordinates
    z = N0.dot(y, N0.transpose(r)) + t

    ## calculate row distances
    return N0.sqrt(N0.sum(N0.power(x - z, 2), 1)) 



#############
##  TESTING        
#############
from . import test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_rmsFit( self ):
        """rmsFit test"""
        from . import tools as T

        self.traj = T.load( T.testRoot('lig_pcr_00/traj.dat') )

        rt, rmsdLst = match( self.traj.ref.xyz, self.traj[-1].xyz)

        if self.local:
            print('RMSD: %.2f' % rmsdLst[0][1])

        # return rotation matrix
        r = abs( N0.sum( N0.ravel( rt[0] )))
        e = abs( N0.sum( N0.ravel( self.EXPECT )))

        self.assertAlmostEqual(r, e, 6)

    EXPECT = N0.array( [[ 0.9999011,   0.01311352,  0.00508244,],
                       [-0.01310219,  0.99991162, -0.00225578,],
                       [-0.00511157,  0.00218896,  0.99998454 ]] )


if __name__ == '__main__':

    BT.localTest()





