## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2019 Raik Gruenberg & Johan Leckner
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
DistanceTrajectory - Reduce trajectory to a set of internal distances 
"""
import copy
import numpy as N
import biskit

class DistanceTrajectory:
    
    def __init__(self, from_atoms=None, to_atoms=None, 
                 n_points=10, refmodel=None, separation=2):
        """
        Prepare object with a given set of (random or manually assigned) 
        atom-atom pairs the distance of which will be tracked.
        
        Args:
            from_atoms ([int]): optional starting atoms for distances [None]
            to_atoms ([int]): optional ending atoms for distances [None]
            n_points (int): number of atom-atom distances to extract [10]
            refmodel (PDBModel): structure to be used for selecting random atoms
            separation (int): minimum residue separation between CA atoms for
                              random distance creation    
        """
        self.from_atoms = from_atoms
        self.to_atoms = to_atoms
        
        if from_atoms is None or to_atoms is None:
            self.from_atoms, self.to_atoms = \
                self.random_atoms(m=refmodel, n=n_points, separation=separation)
            
    
    def random_atoms(self, m, n, separation=2):
        """
        Create a list of random atom pairs.
        
        Args:
            m (biskit.PDBModel): reference PDB
            n (int): number of distance points
            separation (int): minimum separation of residues between start and \
                              end point of each distance
        Returns: 
            tupple: [from_atom_indices], [to_atom_indices]
        """
        
        # position of each CA atom in model
        ca_index = N.compress( m.maskCA(), m.atomRange() ) 

        # First and last residues of reference do not have CA
        # seq = 'XNVTPHHHQHSHIHSYLHLHQQDX'
        # i_res = N.arange(m.lenResidues())
        
        i_res = N.arange(len(ca_index))
        atoms1 = copy.copy(i_res)
        atoms2 = copy.copy(i_res)
        N.random.shuffle(atoms1)
        N.random.shuffle(atoms2)
        
        filtered = N.where(N.abs(atoms1 - atoms2) > 2)[0]
        r1 = N.take(atoms1, filtered[:n])
        r2 = N.take(atoms2, filtered[:n])
                
        ca_1 = N.take(ca_index,r1)
        ca_2 = N.take(ca_index,r2)
        
        return ca_1, ca_2
        
    
    def reduce(self, traj):
        """
        Reduces each frame in a trajectory to a vector with interatomic
        distances

        Args:
            traj (EnsembleTraj): the trajectory that will be reduced

        Returns: 
            numpy.array: N vectors corresponding to the N frames in the
                         trajectory
        """
        print("Reducing traj...")
        t1 = traj.takeAtoms(self.from_atoms)
        t2 = traj.takeAtoms(self.to_atoms)

        distances = N.sqrt(N.sum((t1.frames-t2.frames)**2, axis=2))
        
        return distances, N.array([self.from_atoms, self.to_atoms])


if __name__ == '__main__':

    import biskit.tools as T
    from biskit.md import FuzzyCluster

    ftraj = '~/data/input/traj_step20.dat'
    t = T.load(ftraj) ## Trajectory
    t = t.compressAtoms( t.ref.maskHeavy() )

    d = DistanceTrajectory(n_points=10, refmodel=t.ref)

    v = d.reduce( t )

    fz = FuzzyCluster( v[0:-1:5], n_cluster=10, weight=1.13 )

    centers = fz.go( 1e-10 )

    ## get representative structure for each cluster center:
    msm = fz.getMembershipMatrix()
    i_frames = N.argmax( msm, axis=1 )
    ## models = [ t[i] for i in i_frames ]
    tcenters = t.takeFrames( i_frames )
    tcenters.fit( mask=tcenters.ref.maskCA() )
    tcenters.writePdbs( 'cluster_centers.pdb' )

    ## how many structures per cluster
    frame_membership = N.argmax( msm, axis=0)
    n_members = [ N.sum( frame_membership==i ) for i in range(10) ]

    ## get all frames for each cluster, each into its own trajectory object
    clusters = [ t.takeFrames( N.where(frame_membership==i)[0]) for i in range(10) ]
    
    
