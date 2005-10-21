##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$

import Numeric as N
import tools as T

from mathUtils import aboveDiagonal, SD
from FuzzyCluster import FuzzyCluster
import rmsFit
import types

from Errors import BiskitError

class ClusterError( BiskitError ):
    pass

class TrajCluster:

    def __init__( self, traj ):
        """
        traj - Trajectory to cluster (has to be fitted before hand)
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


    def __raveled( self ):
        """
        Apply current atom mask and return list of raveled frames.
        """
        t = self.traj.compressAtoms( self.aMask )
        return N.array( map( N.ravel, t.frames ) )


    def cluster( self, n_clusters, weight=1.13, converged=1e-11, aMask=None,
                 force=0 ):
        """
        Calculate new clusters.
        n_clusters - number of clusters
        weight - fuzziness weigth
        converged - stop iteration if min dist changes less than converged
        force - re-calculate even if parameters haven't changed
        """
        if aMask == None:
            aMask = N.ones( self.traj.getRef().lenAtoms() )
        
        if self.fc == None or force or self.fcWeight != weight \
           or self.n_clusters != n_clusters or self.aMask != aMask \
           or self.fcConverged != converged:

            self.n_clusters = n_clusters
            self.fcWeight = weight
            self.aMask = aMask

            self.fc = FuzzyCluster( self.__raveled(), self.n_clusters,
                                       self.fcWeight )
        
            self.fcCenters = self.fc.go(self.fcConverged,
                                                    1000, nstep=1)


    def calcClusterNumber( self, min_clst=5, max_clst=30, rmsLimit=1.0,
                        weight=1.13, converged=1e-11, aMask=None, force=0 ):
        """
        Calculate the approximate number of clusters needed to pass
        the average intra-cluster rmsd limit.
        
        min_clst - lower limit for clusters
        max_clst - upper limit for clusters
        rmsLimit - rmsd criteria that the average of all clusters must meet
        weight - fuzziness weigth
        converged - stop iteration if min dist changes less than converged
        force - re-calculate even if parameters haven't changed
          -> int, number of clusters 
        """
        pos = [ min_clst, max_clst ]
        
        while 1:
            clst = int( N.average(pos) )
            self.cluster( clst, weight, converged, aMask, force=1 )
            rmsLst = [ self.avgRmsd(i, aMask)[0] for i in range(clst)]

            if N.average( rmsLst ) > rmsLimit:
                pos[0] = clst
            else:
                pos[1] = clst
                
            if pos[1]-pos[0] == 1:
                T.flushPrint('Converged at %i clusters, current average cluster rmsd %.2f\n'%( clst, N.average( rmsLst ) ))
                return pos[1]

            if pos[1]-pos[0] != 1:
                T.flushPrint('Current cluster setting %i, current average cluster rmsd %.2f\n'%( clst, N.average( rmsLst ) ))
            
            if pos[1]-pos[0]<= 0 or pos[0]<min_clst or pos[1]>max_clst:
                raise ClusterError, "Error determining number of clusters"

        
    def memberships( self ):
        """
        Get degree of membership of each frame to each cluster.
        -> N.array( n_clusters x n_frames )
        """
        return self.fc.getMembershipMatrix()
    

    def maxMemberships(self):
        """
        Get maximum membership value for each frame.
        -> list of float
        """
        msm = self.memberships()
        return map( lambda x: max( msm[:,x]), range(0, self.traj.lenFrames() ))


    def centers( self ):
        """
        Get 'center structure' for each cluster.
        -> N.array( n_clusters x n_atoms_masked x 3 )
        """
        lenAtoms = N.shape( self.fcCenters )[1] / 3
        return N.reshape( self.fcCenters, ( self.n_clusters, lenAtoms, 3))


    def centerFrames( self ):
        """
        Get indices for frame nearest to each cluster center.
        """
        return N.argmax( self.memberships(), 1 )
        

    def memberFrames( self, threshold=0. ):
        """
        Get indices of all frames belonging to each cluster. Each frame
        is guaranteed to belong, at least, to the cluster for which it has
        its maximum membership. If threshold > 0, it can additionally pop
        up in other clusters.
        threshold - float, minimal cluster membership or 0 to consider
                    only max membership
        -> lst of lst of int, n_cluster : frame indices
        """
        ## best cluster for each frame
        msm = self.memberships()
        maxMemb = N.argmax( msm, 0 )

        r = [N.nonzero( N.equal(maxMemb, i) ) for i in range(0, self.n_clusters)]
        r = [ x.tolist() for x in r ]

        ## same thing but now taking all above threshold
        ## -> same frame can end up in several clusters
        if threshold > 0.:
            r2 = [ N.nonzero( N.greater( l, threshold) ) for l in msm ]

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
        frames - list of int
        cluster - int
        -> list of int
        """
        msm = self.memberships()
        
        pairs = [ (-msm[cluster, fr], fr) for fr in frames ]
        pairs.sort()
        return [ x[1] for x in pairs ]
    

    def memberTraj( self, cluster, threshold=0. ):
        """
        Get trajectory with all frames belonging to this cluster, sorted
        by their membership-degree (highest first).
        cluster - int, cluster number
        threshold - float 0-1, minimal cluster membership, see memberFrames()
        -> Trajectory
        """
        mf = self.memberFrames( threshold )[cluster]
        return self.traj.takeFrames( mf )


    def memberTrajs(self, threshold=0.):
        """
        Get member Trajectories for each cluster. Frames are sorted by their
        membership-degree (highest first).
        threshold - float 0-1, minimal cluster membership, see memberFrames()
        -> lst of Trajectories
        """
        mf = self.memberFrames( threshold )
        return [ self.traj.takeFrames(m) for m in mf ]
    

    def avgRmsd( self, cluster, aMask=None, threshold=0. ):
        try:
            rms = self.memberTraj(cluster,threshold).pairwiseRmsd( aMask )
            rms = aboveDiagonal( rms )
        except:
            rms = []
            
        if len(N.ravel(rms)) == 1:
            return N.average(rms)[0], 0.0
        if len(N.ravel(rms)) == 0:
            return 0.0, 0.0
        
        return N.average( rms ), SD( rms )

    def avgRmsd2Ref( self, cluster, ref, avg=1 ):
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

    
## TEST

from Biskit.EnsembleTraj import *
import Biskit.tools

if __name__ == '__main__':

    #traj = Biskit.tools.Load( Biskit.tools.testRoot()+'/lig_pc2_00/traj.dat')
    traj = Biskit.tools.Load('/home/Bis/raik/data/tb/interfaces/c11/lig_pcr_00/traj.dat')
    traj = traj2ensemble( traj )
    
    aMask = traj.ref.mask( lambda a: a['name'] in ['CA','CB','CG'] )

    traj = traj.thin( 1 )

    traj.fit( aMask )

    tc = TrajCluster( traj )

    n_clusters = tc.calcClusterNumber( min_clst=3, max_clst=30,
                                       rmsLimit=0.8, aMask=aMask )

    tc.cluster( n_clusters, aMask=aMask )
