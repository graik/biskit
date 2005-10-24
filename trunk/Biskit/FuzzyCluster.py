##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2002-2004; Wolfgang Rieping
## Copyright (C) 2005; Raik Gruenberg & Johan Leckner
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

## $Revision$
## last $Date$
## last $Author$
"""
Implementation of the fuzzy c-means algorithm
Author: Wolfgang Rieping 1998, 2002

reference:	Heather L. Gordon, Rajmund L. Somorjai
		Fuzzy Cluster Analysis of Molecular Dynamics Trajectories
		PROTEINS: Structure, Function and Genetics 14:249-264 1992
"""

import Numeric as N
import tools

def average(x):
    return N.sum(N.array(x)) / len(x)

def variance(x, avg = None):
    if avg is None:
        avg = N.average(x)

    return N.sum(N.power(N.array(x) - avg, 2)) / (len(x) - 1.)

def standardDeviation(x, avg = None):
    return N.sqrt(variance(x, avg))

def squared_distance_matrix(x, y):

    d1 = N.diagonal(N.dot(x, N.transpose(x)))
    d2 = N.diagonal(N.dot(y, N.transpose(y)))

    a1 = N.add.outer(d1,d2)
    a2 = N.dot(x, N.transpose(y))

    return a1 - 2 * a2

def distance_matrix(x, y):
	return N.sqrt(squared_distance_matrix(x, y))

class FuzzyCluster:
	def __init__(self, data, n_cluster, weight, seedx = 0, seedy = 0):
		self.data = N.array(data, N.Float)
		self.w = weight
		self.n_cluster = n_cluster
		self.npoints, self.dimension = N.shape(data)
		self.seedx = seedx
		self.seedy = seedy

	def calc_membership_matrix(self, d2):
            ## remove 0s (if a cluster center is exactly on one item)
            d2 = N.clip( d2, N.power(1e200, 1-self.w), 1e300 )
            q = N.power(d2, 1. / (1. - self.w))
            return q / N.sum(q)

        def calc_cluster_center(self, msm):
		p = N.power(msm, self.w)
		ccenter = N.transpose(N.dot(p, self.data))
		return N.transpose(ccenter / N.sum(p, 1))

	def updateDistanceMatrix(self):
		return squared_distance_matrix(self.cluster_center,
					       self.data)

	def iterate(self, centers):
		d2 = squared_distance_matrix(centers, self.data)
		msm = self.calc_membership_matrix(d2)
		centers = self.calc_cluster_center(msm)

		return d2, msm, centers

	def error(self, msm, d2):
		p = N.power(msm, self.w)
		product = N.dot(p, N.transpose(d2))
		return N.trace(product)

	def create_membership_matrix(self):
		from RandomArray import seed, random
		
		seed(self.seedx, self.seedy)
		
		r = random((self.npoints, self.n_cluster))
		return N.transpose(r / N.sum(r))

	def go(self, errorthreshold, n_iterations = 1e10, nstep = 10):
		iteration = 0
		rel_err = 1e10
		error = 1.e10

		msm = self.create_membership_matrix()
		centers = self.calc_cluster_center(msm)
		
		while rel_err > errorthreshold and iteration<n_iterations:
			d2, msm, centers = self.iterate(centers)
			old_error = error
			error = self.error(msm, d2)
			rel_err = abs(1. - error/old_error)
			iteration = iteration+1
			if not iteration % nstep:
				tools.errWrite( "%i %f\n" % (iteration, error) )

		self.centers = centers
		self.msm = msm
		self.d2 = d2
		
		return centers

	def clusterEntropy(self):
		centropy = N.diagonal(N.dot(self.msm,
				       N.transpose(N.log(self.msm))))
		return -1/float(self.npoints)*centropy

	def entropy(self):
		return N.sum(self.clusterEntropy())
		
	def nonFuzzyIndex(self):
		p = N.power(self.msm, self.w)
		return (self.n_cluster*N.sum(N.sum(p))-
			self.npoints)/(self.npoints*(self.n_cluster-1))

	def clusterPartitionCoefficient(self):
		return N.sum(N.power(self.msm, self.w), 1)/self.npoints

	def partitionCoefficient(self):
		return N.sum(self.clusterPartitionCoefficient())

	def getMembershipMatrix(self):
		return self.msm

	def getClusterCenter(self):
		return self.cluster_center

        def entropySD(self):
            centropy = N.sum(-N.log(self.msm)*\
                           self.msm)/float(self.n_cluster)
            return standardDeviation(centropy)

        def standardDeviation(self):
            sd = standardDeviation(self.msm)
            return sd
    
if __name__ == '__main__':

	from RandomArray import *
	from gnuplot import *
	from time import clock

	x1 = random((500,2))
	x2 = random((500,2)) + 1
	x3 = random((500,2)) + 2

	x = N.concatenate((x1, x2, x3))

	fc = FuzzyCluster(x, n_cluster = 5, weight = 1.5)

	start = clock()
	
	centers = fc.go(1.e-30, n_iterations=50, nstep=10)

	print clock() - start, 's'
	
