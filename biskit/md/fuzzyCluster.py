## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2002-2004; Wolfgang Rieping
## Copyright (C) 2005; Raik Gruenberg & Johan Leckner
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


"""
Implementation of the fuzzy c-means algorithm

Author: Wolfgang Rieping 1998, 2002

Reference::
  Heather L. Gordon, Rajmund L. Somorjai
  Fuzzy Cluster Analysis of Molecular Dynamics Trajectories
  PROTEINS: Structure, Function and Genetics 14:249-264 1992
"""

import biskit.core.oldnumeric as N0
import biskit.mathUtils as MU
import biskit.tools as tools

import numpy.random.mtrand as R # seed, random / converted from oldnumeric/random_array

## def average(x):
##     return N0.sum(N0.array(x)) / len(x)

## def variance(x, avg = None):
##     if avg is None:
##         avg = N0.average(x)

##     return N0.sum(N0.power(N0.array(x) - avg, 2)) / (len(x) - 1.)

## def standardDeviation(x, avg = None):
##     return N0.sqrt(variance(x, avg))

def squared_distance_matrix(x, y):

    d1 = N0.diagonal(N0.dot(x, N0.transpose(x)))
    d2 = N0.diagonal(N0.dot(y, N0.transpose(y)))

    a1 = N0.add.outer(d1,d2)
    a2 = N0.dot(x, N0.transpose(y))

    return a1 - 2 * a2


def distance_matrix(x, y):
    return N0.sqrt(squared_distance_matrix(x, y))


class FuzzyCluster:
    def __init__(self, data, n_cluster, weight, seedx = 0, seedy = 0):
        """
        @param data: cluster this
        @type  data: [float] OR array
        @param n_cluster: number of clusters
        @type  n_cluster: int
        @param weight: fuzziness weigth
        @type  weight: float
        @param seedx: random seed value for RandomArray.seed (default: 0)
        @type  seedx: int OR 0
        @param seedy: random seed value for RandomArray.seed
        (default: 0, set seed from clock)
        @type  seedy: int OR 0
        """
        self.data = N0.array(data, N0.Float)
        self.w = weight
        self.n_cluster = n_cluster
        self.npoints, self.dimension = N0.shape(data)
        self.seedx = seedx
        self.seedy = seedy


    def calc_membership_matrix(self, d2):
        ## remove 0s (if a cluster center is exactly on one item)
        d2 = N0.clip( d2, N0.power(1e200, 1-self.w), 1e300 )
        q = N0.power(d2, 1. / (1. - self.w))
        return q / N0.sum(q)


    def calc_cluster_center(self, msm):
        p = N0.power(msm, self.w)
        ccenter = N0.transpose(N0.dot(p, self.data))
        return N0.transpose(ccenter / N0.sum(p, 1))


    def updateDistanceMatrix(self):
        return squared_distance_matrix(self.cluster_center, self.data)


    def iterate(self, centers):
        """
        @param centers: array with cluster centers
        @type  centers: array('f')

        @return: distance to the centers, membership matrix, array of cenetrs
        @rtype: array, array, array
        """
        d2 = squared_distance_matrix(centers, self.data)
        msm = self.calc_membership_matrix(d2)
        centers = self.calc_cluster_center(msm)

        return d2, msm, centers


    def error(self, msm, d2):
        """
        @param msm: membership matrix
        @type  msm: array('f')
        @param d2: distance from data to the centers
        @type  d2: array('f')

        @return: weighted error 
        @rtype: float
        """
        p = N0.power(msm, self.w)
        product = N0.dot(p, N0.transpose(d2))
        return N0.trace(product)


    def create_membership_matrix(self):
        """
        Create a random membership matrix.

        @return: random array of shape length of data to
                 cluster times number of clusters
        @rtype: array('f')
        """
        ## default signature has changed oldnumeric->numpy
        if (self.seedx==0 or self.seedy==0):  
            R.seed()
        else:
            R.seed((self.seedx, self.seedy))

        r = R.random_sample((self.npoints, self.n_cluster))
        return N0.transpose(r / N0.sum(r))


    def go(self, errorthreshold, n_iterations=1e10, nstep=10, verbose=1):
        """
        Start the cluestering. Run until the error is below the error
        treshold or the max number of iterations have been run.

        @param errorthreshold: treshold value for error 
        @type  errorthreshold: float
        @param n_iterations: treshold value for number of iterations
                             (default: 1e10)
        @type  n_iterations: int
        @param nstep: print information for every n'th step in the iteration
        @type  nstep: int

        @return: array with cluster centers
        @rtype: array('f')
        """
        iteration = 0
        rel_err = 1e10
        error = 1.e10

        msm = self.create_membership_matrix()
        centers = self.calc_cluster_center(msm)

        while rel_err > errorthreshold and iteration < n_iterations:
            d2, msm, centers = self.iterate(centers)
            old_error = error
            error = self.error(msm, d2)
            rel_err = abs(1. - error/old_error)
            iteration = iteration+1
            if not iteration % nstep and verbose:
                tools.errWrite( "%i %f\n" % (iteration, error) )

        self.centers = centers
        self.msm = msm
        self.d2 = d2

        return centers


    def clusterEntropy(self):
        centropy = N0.diagonal(N0.dot(self.msm,
                                    N0.transpose(N0.log(self.msm))))
        return -1/float(self.npoints)*centropy


    def entropy(self):
        return N0.sum(self.clusterEntropy())


    def nonFuzzyIndex(self):
        p = N0.power(self.msm, self.w)
        return (self.n_cluster*N0.sum(N0.sum(p))-
                self.npoints)/(self.npoints*(self.n_cluster-1))

    def clusterPartitionCoefficient(self):
        return N0.sum(N0.power(self.msm, self.w), 1)/self.npoints


    def partitionCoefficient(self):
        return N0.sum(self.clusterPartitionCoefficient())


    def getMembershipMatrix(self):
        return self.msm


    def getClusterCenter(self):
        return self.cluster_center


    def entropySD(self):
        centropy = N0.sum(-N0.log(self.msm)*\
                         self.msm)/float(self.n_cluster)
        return MU.SD(centropy)


    def standardDeviation(self):
        sd = MU.SD(self.msm)
        return sd


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """FuzzyCluster test"""    

    def test_FuzzyCluster( self):
        """FuzzyCluster test"""
        import biskit.gnuplot as G

        x1 = R.random_sample((500,2))
        x2 = R.random_sample((500,2)) + 1
        x3 = R.random_sample((500,2)) + 2

        self.x = N0.concatenate((x1, x2, x3))

        self.fuzzy = FuzzyCluster(self.x, n_cluster=5, weight=1.5)

        self.centers = self.fuzzy.go(1.e-30, n_iterations=50, nstep=10,
                                     verbose=self.local)

        if self.local:
            print("cluster centers are displayed in green")
            G.scatter( self.x, self.centers )

        self.assertEqual( N0.shape(self.centers), (5, 2) )

if __name__ == '__main__':

    BT.localTest()
