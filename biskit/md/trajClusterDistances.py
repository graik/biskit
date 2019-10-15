## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2019 Raik Gruenberg
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

"""
Cluster the members of a trajectory by intra-molecular distances.
"""
import numpy as N
import biskit.tools as T

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.md; __package__ = "biskit.md"

from .trajClusterRmsd import TrajClusterRmsd, ClusterError

class TrajClusterDistances(TrajClusterRmsd):
    """
    Cluster molecular dynamics trajectory by a (random) set of intra-molecular
    distances. This clustering does not require fitting of structure data and
    is more appropriate for very diverse structure ensembles.
    
    Clustering is performed using the c-means fuzzy clustering implemented in
    biskit.md.fuzzyCluster.
    """

    def __init__( self, traj, verbose=1, aMask=None ):
        """
        @param traj: Trajectory to cluster (has to be fitted before hand)
        @type  traj: Trajectory
        @param aMask: atom mask to be applied before clustering (default None)
        @type aMask: numpy.array(bool)
        """
        super(self).__init__(traj, verbose=verbose, aMask=aMask)

    def data4clustering( self ):
        """
        Apply current atom mask and extract a vector of intra-molecular 
        distances for each frame of the trajectory.
        @return: [float] or numpy.array(float)
        """
        pass


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test Adaptive clustering"""

    def test_TrajClusterDistances(self):
        """TrajClusterDistances test"""
        from biskit.md import traj2ensemble

        traj = T.load( T.testRoot()+'/lig_pcr_00/traj.dat')

        traj = traj2ensemble( traj )

        aMask = traj.ref.mask( lambda a: a['name'] in ['CA','CB','CG'] )

        traj = traj.thin( 1 )

        traj.fit( aMask, verbose=self.local )
        self.tc = TrajClusterRmsd( traj, verbose=self.local )

        ## check how many clusters that are needed with the given criteria
        n_clusters = self.tc.calcClusterNumber( min_clst=3, max_clst=15,
                                                rmsLimit=0.7, aMask=aMask )

        ## cluster
        self.tc.cluster( n_clusters, aMask=aMask )

        if self.local:
            member_frames = self.tc.memberFrames()

            print('There are %i clusters where the members are:'%n_clusters)
            for i in range(n_clusters):
                print('Cluster %i (%i members): %s'%( i+1,
                                                      len(member_frames[i]),
                                                      member_frames[i] ))


if __name__ == '__main__':

    BT.localTest()

