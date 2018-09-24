## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
Create structures with reduced number of atoms.
"""

from biskit import PDBModel, DictList
from biskit.core import oldnumeric as N0
from biskit import tools as T
from biskit import molUtils as MU


class ReduceCoordinates:
    """
    ReduceCoordinates
    =================
    
    Translate a PDBModel or frames from a trajectory to structure(s) with
    only one backbone and up to 2 side chain atoms per residue.
    The new atoms are the centers of mass of several atoms and carry the
    weight of the pooled atoms in an atom profile called 'mass'.

    Examples
    --------
      >>> ## create with reference PDBModel
      >>> reducer = ReduceCoordinates( m_ref )
      >>> ## creates reduced PDBModel from m_ref
      >>> m_red = reducer.reduceToModel()

    OR:
      >>> m_red_1 = reducer.reduceToModel( m1.getXyz() ) ## reduce many models
      >>> m_red_2 = reducer.reduceToModel( m2.getXyz() ) ## with identical atoms

    OR:
      >>> ## reduce a complete Trajectory
      >>> reducer = ReduceCoordinates( traj.ref )
      >>> red_ref= reducer.reduceToModel()
      >>> frames = reducer.reduceXyz( traj.frames )
      >>> traj_red = Trajectory( ref=red_ref )
      >>> traj_red.frames = frames
    """
    ## modify order of TYR/PHE ring atoms to move centers away from ring axis
    aaAtoms = MU.aaAtoms
    aaAtoms['TYR'] = ['N','CA','C','O','CB','CG','CD1','CE1','CD2',
                      'CE2','CZ','OH', 'OXT']
    aaAtoms['PHE'] = ['N','CA','C','O','CB','CG','CD1','CE1','CD2',
                      'CE2','CZ', 'OXT']

    def __init__( self, model,  maxPerCenter=4 ):
        """
        Prepare reduction of coordinates from a given model.
        
        @param model: reference model defining atom content and order
        @type  model: PDBModel
        @param maxPerCenter: max number of atoms per side chain center atom
                             (default: 4)
        @type  maxPerCenter: int
        """
        self.m = model
        self.__addMassProfile( self.m )

        ## sort atoms within residues into standard order
        def cmpAtoms( a1, a2 ):
            """
            Comparison function for bringing atoms into standard order
            within residues as defined by L{ aaAtoms }.

            @param a1: model
            @type  a1: PDBModel
            @param a2: model
            @type  a2: PDBModel

            @return: int or list of matching positions
            @rtype: [-1|0|1]            
            """
            ## cmp vanished in python 3.x (but still available in past.builtins)
            cmp = lambda x, y: (x > y) - (x < y)

            res = a1['residue_name']
            target = self.aaAtoms[ res ]
            try:
                return cmp(target.index( a1['name'] ), target.index( a2['name'] ))
            except ValueError as why:
                return cmp( a1['name'], a2['name'] )
##                 s = "Unknown atom for %s %i: %s or %s" % \
##                     (res, a1['residue_number'], a1['name'], a2['name'] )
##                 raise PDBError( s )

        self.a_indices = self.m.argsort( cmpAtoms )
        self.m_sorted = self.m.sort( self.a_indices )

        ## remove H from internal model and from list of atom positions
        maskH = self.m_sorted.remove( self.m_sorted.maskH() )
        self.a_indices = N0.compress( maskH, self.a_indices )

        self.makeMap( maxPerCenter )


    def __addMassProfile( self, model ):
        """
        Add a mass profile to the model.

        @param model: model
        @type  model: PDBModel        
        """
        h = MU.atomMasses['H']
        masses = []
        
        elements = model['element']
        resnames = model['residue_name']
        anames   = model['name']
        
        for i in model.atomRange():
            if elements[i] == 'H':
                masses += [ 0 ]
            else:
                masses += [ MU.atomMasses[ elements[i] ] +
                            h * MU.aaAtomsH[ resnames[i] ][ anames[i] ] ]
        
        
        #for a in model.atoms:
            #if a['element'] == 'H':
                #masses += [ 0 ]
            #else:
                #masses += [ MU.atomMasses[ a['element'] ] +
                            #h * MU.aaAtomsH[a['residue_name']][a['name']] ]

        model.atoms.set( 'mass', masses, 
                             comment='mass in D, hydrogen mass added to heavy' )


    def group( self, a_indices, maxPerCenter ):
        """
        Group a bunch of integers (atom indices in PDBModel) so that each
        group has at most maxPerCenter items.
        
        @param a_indices: atom indices
        @type  a_indices: [int]
        @param maxPerCenter: max entries per group
        @type  maxPerCenter: int
        
        @return: list of lists of int
        @rtype: [[int],[int]..]
        """
        ## how many groups are necessary?
        n_centers = len( a_indices ) // maxPerCenter  ## floor division
        if len( a_indices ) % maxPerCenter:
            n_centers += 1

        ## how many items/atoms go into each group?
        nAtoms = N0.ones(n_centers, N0.Int) * int(len( a_indices ) / n_centers)
        i=0
        while N0.sum(nAtoms) != len( a_indices ):
            nAtoms[i] += 1
            i += 1

        ## distribute atom indices into groups
        result = []
        pos = 0
        for n in nAtoms:
            result += [ N0.take( a_indices, N0.arange(n) + pos) ]
            pos += n

        return result


    def nextAtom( self, atom, name ):
        """
        Create an atom dictionary.
        @param atom:
        @type  atom: L{Biskit.CrossView}
        @param name: atom name
        @type  name: str
        
        @return: atom dictionary
        @rtype: dict
        """
        self.currentAtom += 1
        a = atom.toDict()

        a['name'] = name
        a['serial_number'] = self.currentAtom
        return a

    def makeMap( self, maxPerCenter=4 ):
        """
        Calculate mapping between complete and reduced atom list.
        Creates a (list of lists of int, list of atom dictionaries)
        containing groups of atom indices into original model, new center atoms
        
        @param maxPerCenter: max number of atoms per side chain center atom
                             (default: 4)
        @type  maxPerCenter: int
        """

        resIndex = self.m_sorted.resIndex()
        resModels= self.m_sorted.resModels()
        m = self.m_sorted

        self.currentAtom = 0

        groups = []
        atoms = DictList()

        for i in range( len( resIndex ) ):

            first_atom = resIndex[ i ]

            if i < len( resIndex )-1:
                last_atom  = resIndex[ i+1 ] - 1
            else:
                last_atom = len( self.a_indices ) - 1

            a = m.atoms[ first_atom ]

##             res_name  = m.atoms[ first_atom ]['residue_name']
##             segid     = m.atoms[ first_atom ]['segment_id']
##             chainId   = m.atoms[ first_atom ]['chain_id']
##             res_number= m.atoms[ first_atom ]['serial_number']

            ## position of this residue's atoms in original PDBModel (unsorted)
            a_indices = self.a_indices[ first_atom : last_atom+1 ]

            ## for each center create list of atom indices and a center atom
            if a['residue_name'] != 'GLY' and a['residue_name'] != 'ALA':

                bb_a_indices = N0.compress( resModels[i].maskBB(), a_indices)
                sc_a_indices = N0.compress(
                    N0.logical_not( resModels[i].maskBB()), a_indices )

                sc_groups = self.group( sc_a_indices, maxPerCenter )

            else:
                bb_a_indices = a_indices
                sc_groups = []

            groups += [ bb_a_indices ]
            atoms  += [ self.nextAtom(a, 'BB') ]

            i = 0
            for g in sc_groups:
                groups += [ g ]
                atoms  += [ self.nextAtom( a, 'SC%i'%i) ]
                i += 1

        self.groups = groups
        self.atoms = atoms


    def reduceXyz( self, xyz, axis=0 ):
        """
        Reduce the number of atoms in the given coordinate set. The set must
        have the same length and order as the reference model. It may have
        an additional (time) dimension as first axis.
        
        @param xyz: coordinates (N_atoms x 3) or (N_frames x N_atoms x 3)
        @type  xyz: array
        @param axis: axis with atoms (default: 0)
        @type  axis: int
        
        @return: coordinate array (N_less_atoms x 3) or
                 (N_frames x N_less_atoms x 3)
        @rtype: array
        """
        masses = self.m.atoms.get('mass')
        r_xyz = None

        for atom_indices in self.groups:

            x = N0.take( xyz, atom_indices, axis )
            m = N0.take( masses, atom_indices )

            center = N0.sum( x * N0.transpose([m,]), axis=axis) / N0.sum( m )

            if axis == 0:
                center = center[N0.NewAxis, :]

            if axis == 1:
                center = center[:, N0.NewAxis, :]

            if r_xyz is None:
                r_xyz = center

            else:
                r_xyz = N0.concatenate( (r_xyz, center), axis )

        return r_xyz


    def reduceToModel( self, xyz=None, reduce_profiles=1  ):
        """
        Create a reduced PDBModel from coordinates. Atom profiles the source
        PDBModel are reduced by averaging over the grouped atoms.
        
        @param xyz: coordinte array (N_atoms x 3) or
                    None (->use reference coordinates)
        @type  xyz: array OR None
        
        @return: PDBModel with reduced atom set and profile 'mass'
        @rtype: PDBModel
        """

        mass = self.m.atoms.get('mass')
        if xyz is None: xyz = self.m.getXyz()

        mProf = [ N0.sum( N0.take( mass, group ) ) for group in self.groups ]
        xyz = self.reduceXyz( xyz )

        result = PDBModel()

        for k in self.atoms.keys():
            result.atoms.set( k, self.atoms.valuesOf(k) )

##         result.setAtoms( self.atoms )

        result.setXyz( xyz )
        result.atoms.set( 'mass', mProf )

        if reduce_profiles:
            self.reduceAtomProfiles( self.m, result )

            result.residues = self.m.residues

        return result


    def reduceAtomProfiles( self, from_model, to_model ):
        """
        reduce all atom profiles according to the calculated map by calculating
        the average over the grouped atoms.
        
        @param from_model: model
        @type  from_model: PDBModel
        @param to_model: model
        @type  to_model: PDBModel
        """
        for profname in from_model.atoms:

            p0 =  from_model.atoms.get(profname)
            info = from_model.profileInfo( profname )

            try:
                pr = [ N0.average( N0.take( p0, group ) ) for group in self.groups ]

                to_model.atoms.set( profname, pr )
            except:
                pass
                
            to_model.atoms.setInfo( profname, **info )



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    def test_ReduceCoordinates(self):
        """ReduceCoordinates test"""

        self.m = PDBModel( T.testRoot()+'/com/1BGS.pdb' )
        self.m = self.m.compress( N0.logical_not( self.m.maskH2O() ) )

        self.m.atoms.set('test', list(range(len(self.m))))

        self.red = ReduceCoordinates( self.m, 4 )

        self.mred = self.red.reduceToModel()
        
        if self.local:
            print('\nAtoms before reduction %i'% self.m.lenAtoms())
            print('Atoms After reduction %i'% self.mred.lenAtoms())

        self.assertEqual( self.mred.lenAtoms(), 445 )

if __name__ == '__main__':

    BT.localTest()
    
