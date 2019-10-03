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
                 n_random=0, refmodel=None, separation=2):
        self.from_atoms = from_atoms
        self.to_atoms = to_atoms
        
        if from_atoms is None or to_atoms is None:
            self.from_atoms, self.to_atoms = \
                self.random_atoms(m=refmodel, n=n_atoms, separation=separation)
            
    
    def random_atoms(self, m, n, separation=2):
        """
        Args:
           -m : biskit.PDBModel, reference PDB
           -n : int, number of distance points
           -separation : int, minimum separation of residues between start and \
                              end point of each distance
        Returns: [from_atom_indices], [to_atom_indices]
        """
        i_res = N.arange(m.lenResidues())
        atoms1 = copy.copy(i_res)
        atoms2 = copy.copy(i_res)
        N.shuffle(atoms1)
        N.shuffle(atoms2)
        
        filtered = N.where(N.abs(atoms1 - atoms2) > 2)[0]
        r1 = N.take(atoms1, filtered[n])
        r2 = N.take(atoms2, filtered[n])
        
        ca_index = N.compress( m.maskCA(), m.atomRange() ) # position of each CA atom in model
        
        ca_1 = N.take(ca_index,r1)
        ca_2 = N.take(ca_index,r2)
        
        return ca_1, ca_2
        
    
    def reduce(self, traj):
        t1 = traj.takeAtoms(self.atoms_from)
        t2 = traj.takeAtoms(self.atoms_to)
        
        
        
        return r
