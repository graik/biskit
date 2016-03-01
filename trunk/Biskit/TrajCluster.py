## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
##
##
## last $Author$
## last $Date$
## $Revision$

"""
Cluster the members of a trajectory.
"""

import numpy as N
import tools as T

from mathUtils import aboveDiagonal, SD
from FuzzyCluster import FuzzyCluster
import rmsFit
import types

from Errors import BiskitError

class ClusterError( BiskitError ):
    pass

class TrajCluster:

    def __init__( self, traj, verbose=1 ):
        """
        @param traj: Trajectory to cluster (has to be fitted before hand)
        @type  traj: Trajectory
        """
        self.traj = traj
        ## atom mask applied before clustering
        self.aMask = None

        ## last clustering parameters
        self.n_clusters = 0
        self.fcWeight = 1.13
        self.fcConverged = 1e-10

        ## last clustering results
        self.fc = None 
        self.fcCenters = None

        ## by default prints clustering status to stdout
        self.verbose = verbose


    def __raveled( self ):
        """
        Apply current atom mask and return list of raveled frames.
        """
        t = self.traj.compressAtoms( self.aMask )
        return N.array( map( N.ravel, t.frames ) )


    def cluster( self, n_clusters, weight=1.13, converged=1e-11,
                 aMask=None, force=0 ):
        """
        Calculate new clusters.

        @param n_clusters: number of clusters
        @type  n_clusters: int
        @param weight: fuzziness weigth
        @type  weight: float (default: 1.13)
        @param converged: stop iteration if min dist changes less than
                          converged (default: 1e-11)
        @type  converged: float
        @param aMask: atom mask applied before clustering
        @type  aMask: [1|0]
        @param force: re-calculate even if parameters haven't changed
                      (default:0)
        @type  force: 1|0
        """
        if aMask == None:
            aMask = N.ones( self.traj.getRef().lenAtoms() )

        if self.fc == None or force or self.fcWeight != weight \
           or self.n_clusters != n_clusters or N.any( self.aMask != aMask) \
           or self.fcConverged != converged:

            self.n_clusters = n_clusters
            self.fcWeight = weight
            self.aMask = aMask

            self.fc = FuzzyCluster( self.__raveled(), self.n_clusters,
                                    self.fcWeight )

            self.fcCenters = self.fc.go( self.fcConverged,
                                         1000, nstep=10,
                                         verbose=self.verbose )


    def calcClusterNumber( self, min_clst=5, max_clst=30, rmsLimit=1.0,
                           weight=1.13, converged=1e-11, aMask=None, force=0 ):
        """
        Calculate the approximate number of clusters needed to pass
        the average intra-cluster rmsd limit.

        @param min_clst: lower limit for clusters (default: 5)
        @type  min_clst: int
        @param max_clst: upper limit for clusters (default: 30 )
        @type  max_clst: int
        @param rmsLimit: rmsd criteria that the average of all clusters
                         must meet in Angstrom (default: 1.0)
        @type  rmsLimit: float
        @param weight: fuzziness weigth (default: 1.13)
        @type  weight: float
        @param converged: stop iteration if min dist changes less than
                          converged (default: 1e-11)
        @type  converged: float
        @param force: re-calculate even if parameters haven't changed
                      (default: 0)
        @type  force: 1|0

        @return: number of clusters
        @rtype: int

        @raise ClusterError: if can't determining number of clusters
        """
        pos = [ min_clst, max_clst ]

        while 1:
            clst = int( N.average(pos) )
            self.cluster( clst, weight, converged, aMask, force=force )
            rmsLst = [ self.avgRmsd(i, aMask)[0] for i in range(clst)]

            if N.average( rmsLst ) > rmsLimit:
                pos[0] = clst
            else:
                pos[1] = clst

            if pos[1]-pos[0] == 1:
                if self.verbose:
                    T.flushPrint('Converged at %i clusters, current average cluster rmsd %.2f\n'%( clst, N.average( rmsLst ) ))
                return pos[1]

            if pos[1]-pos[0] != 1:
                if self.verbose:
                    T.flushPrint('Current cluster setting %i, current average cluster rmsd %.2f\n'%( clst, N.average( rmsLst ) ))

            if pos[1]-pos[0]<= 0 or pos[0]<min_clst or pos[1]>max_clst:
                raise ClusterError, "Error determining number of clusters"


    def memberships( self ):
        """
        Get degree of membership of each frame to each cluster.

        @return: N.array( n_clusters x n_frames )
        @rtype: array
        """
        return self.fc.getMembershipMatrix()


    def maxMemberships(self):
        """
        Get maximum membership value for each frame.

        @return: list of float
        @rtype: [float]
        """
        msm = self.memberships()
        return map( lambda x: max( msm[:,x]), range(0, self.traj.lenFrames() ))


    def centers( self ):
        """
        Get 'center structure' for each cluster.

        @return: N.array( n_clusters x n_atoms_masked x 3 )
        @rtype: array
        """
        lenAtoms = N.shape( self.fcCenters )[1] / 3
        return N.reshape( self.fcCenters, ( self.n_clusters, lenAtoms, 3))


    def centerFrames( self ):
        """
        Get indices for frame nearest to each cluster center.

        @return: list of cluster center indecies
        @rtype: [int]
        """
        return N.argmax( self.memberships(), 1 )


    def memberFrames( self, threshold=0. ):
        """
        Get indices of all frames belonging to each cluster. Each frame
        is guaranteed to belong, at least, to the cluster for which it has
        its maximum membership. If threshold > 0, it can additionally pop
        up in other clusters.

        @param threshold: minimal cluster membership or 0 to consider
                          only max membership (default: 0)
        @type  threshold: float

        @return: n_cluster, lst of lst of int, frame indices
        @rtype: [[int]]
        """
        ## best cluster for each frame
        msm = self.memberships()
        maxMemb = N.argmax( msm, 0 )

        r = [N.nonzero( N.equal(maxMemb, i) )[0] for i in range(0, self.n_clusters)]
        r = [ x.tolist() for x in r ]

        ## same thing but now taking all above threshold
        ## -> same frame can end up in several clusters
        if threshold > 0.:
            r2 = [ N.nonzero( N.greater( l, threshold) )[0] for l in msm ]

            ## add only additional frames
            for i in range(0, len( r ) ):
                try:
                    frames = r[i].tolist()
                except:
                    frames = r[i]

                r[i] = frames + [ fr for fr in r2[i] if fr not in r[i] ]

        ## sort frames within each cluster by their membership
        r = [ self.membershipSort( r[i], i) for i in range(0, len(r) )]

        return r


    def membershipSort( self, frames, cluster ):
        """
        Sort given list of frame indices by their membership in cluster.

        @param frames: list of frame indecies
        @type  frames: [int]
        @param cluster: cluster number
        @type  cluster: int

        @return: indecies sorted by ther membership to cluster
        @rtype: [int]
        """
        msm = self.memberships()

        pairs = [ (-msm[cluster, fr], fr) for fr in frames ]
        pairs.sort()
        return [ x[1] for x in pairs ]


    def memberTraj( self, cluster, threshold=0. ):
        """
        Get trajectory with all frames belonging to this cluster, sorted
        by their membership-degree (highest first).

        @param cluster: cluster number
        @type  cluster: int
        @param threshold: float 0-1, minimal cluster membership,
                          see L{memberFrames()}
        @type  threshold: float

        @return: Trajectory with all members of a cluster, sorted by membership
        @rtype: Trajectory
        """
        mf = self.memberFrames( threshold )[cluster]
        return self.traj.takeFrames( mf )


    def memberTrajs(self, threshold=0.):
        """
        Get member Trajectories for each cluster. Frames are sorted by their
        membership-degree (highest first).

        @param threshold: float 0-1, minimal cluster membership,
                          see L{memberFrames()}
        @type  threshold:

        @return: lst of Trajectories, with members of each cluster,
                 sorted by membership
        @rtype: [Trajectory]
        """
        mf = self.memberFrames( threshold )
        return [ self.traj.takeFrames(m) for m in mf ]


    def avgRmsd( self, cluster, aMask=None, threshold=0. ):
        """
        Claculate the average pairwise rmsd (in Angstrom) for members
        of a cluter.

        @param cluster: cluster number
        @type  cluster: int       
        @param aMask: atom mask applied before calculation
        @type  aMask: [1|0]
        @param threshold: float 0-1, minimal cluster membership,
                          see L{memberTraj()}
        @type  threshold: float

        @return: average rmsd and the standard deviation 
        @rtype: float, float   
        """
        try:
            rms = self.memberTraj(cluster,threshold).pairwiseRmsd( aMask )
            rms = aboveDiagonal( rms )
        except:
            rms = []

        if len(N.ravel(rms)) == 1:
            ## was: return N.average(rms)[0], 0.0
            return N.average(rms), 0.0
        if len(N.ravel(rms)) == 0:
            return 0.0, 0.0

        return N.average( rms ), SD( rms )


    def avgRmsd2Ref( self, cluster, ref, avg=1 ):
        """
        Claculate the rmsd (or average rmsd) of all frames belonging to a
        cluster to a reference structure (in Angstrom).

        @param cluster: cluster number
        @type  cluster: int
        @param ref: reference structure
        @type  ref: model
        @param avg: return the average rmsd (1) OR a list with all rmsds (0)
                    (default: 1)
        @type  avg: float OR [float]

        """
        eTraj = self.memberTraj(cluster,threshold=0)
        rms = []

        if type(ref) == types.InstanceType:
            ref = ref.xyz

        for frame in eTraj.frames:
            rt, rmsdLst = rmsFit.match( ref, frame)
            rms += [ rmsdLst[0][1] ]
        if avg ==1:
            return  N.average(rms)
        return rms


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test Adaptive clustering"""

    def test_TrajCluster(self):
        """TrajCluster test"""
        from Biskit.EnsembleTraj import traj2ensemble

        traj = T.load( T.testRoot()+'/lig_pcr_00/traj.dat')

        traj = traj2ensemble( traj )

        aMask = traj.ref.mask( lambda a: a['name'] in ['CA','CB','CG'] )

        traj = traj.thin( 1 )

        traj.fit( aMask, verbose=self.local )
        self.tc = TrajCluster( traj, verbose=self.local )

        ## check how many clusters that are needed with the given criteria
        n_clusters = self.tc.calcClusterNumber( min_clst=3, max_clst=15,
                                                rmsLimit=0.7, aMask=aMask )

        ## cluster
        self.tc.cluster( n_clusters, aMask=aMask )

        if self.local:
            member_frames = self.tc.memberFrames()

            print 'There are %i clusters where the members are:'%n_clusters
            for i in range(n_clusters):
                print 'Cluster %i (%i members): %s'%( i+1,
                                                      len(member_frames[i]),
                                                      member_frames[i] )


if __name__ == '__main__':

    BT.localTest()

