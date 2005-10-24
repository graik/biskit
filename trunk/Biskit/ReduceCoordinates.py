##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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
##
##
## $Revision$
## last $Date$
## last $Author$


from PDBModel import PDBModel
import Numeric as N
import tools as T
import molUtils as MU


class ReduceCoordinates:
    """
    Translate a PDBModel or frames from a trajectory to structure(s) with
    only one backbone and up to 2 side chain atoms per residue.
    The new atoms are the centers of mass of several atoms and carry the
    weight of the pooled atoms in an atom profile called 'mass'.

    Examples:
    reducer = ReduceCoordinates( m_ref ) ## create with reference PDBModel
    m_red = reducer.reduceToModel()   ## creates reduced PDBModel from m_ref

    OR:
    ..
    m_red_1 = reducer.reduceToModel( m1.getXyz() ) ## reduce many models
    m_red_2 = reducer.reduceToModel( m2.getXyz() ) ## with identical atoms

    OR:
    ## reduce a complete Trajectory
    reducer = ReduceCoordinates( traj.ref )
    red_ref= reducer.reduceToModel()
    frames = reducer.reduceXyz( traj.frames )
    traj_red = Trajectory( ref=red_ref )
    traj_red.frames = frames
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
        model - PDBModel, reference model defining atom content and order
        maxPerCenter - max number of atoms per side chain center atom
        """
        self.m = model
        self.__addMassProfile( self.m )

        ## sort atoms within residues into standard order
        def cmpAtoms( a1, a2 ):
            res = a1['residue_name']
            target = self.aaAtoms[ res ]
            try:
                return cmp(target.index( a1['name'] ),
                           target.index( a2['name'] ))
            except ValueError, why:
                return cmp( a1['name'], a2['name'] )
##                 s = "Unknown atom for %s %i: %s or %s" % \
##                     (res, a1['residue_number'], a1['name'], a2['name'] )
##                 raise PDBError( s )

        self.a_indices = self.m.argsort( cmpAtoms )
        self.m_sorted = self.m.sort( self.a_indices )

        ## remove H from internal model and from list of atom positions
        maskH = self.m_sorted.remove( self.m_sorted.maskH() )
        self.a_indices = N.compress( maskH, self.a_indices )

        self.makeMap( maxPerCenter )


    def __addMassProfile( self, model ):
        h = MU.atomMasses['H']
        masses = []
        for a in model.atoms:
            if a['element'] == 'H':
                masses += [ 0 ]
            else:
                masses += [ MU.atomMasses[ a['element'] ] +
                            h * MU.aaAtomsH[a['residue_name']][a['name']] ]
        
        model.setAtomProfile( 'mass', masses )


    def group( self, a_indices, maxPerCenter ):
        """
        Group a bunch of integers (atom indices in PDBModel) so that each
        group has at most maxPerCenter items.
        a_indices - list of int
        maxPerCenter - int, max entries per group
        -> list of lists of int
        """
        ## how many groups are necessary?
        n_centers = len( a_indices ) / maxPerCenter
        if len( a_indices ) % maxPerCenter:
            n_centers += 1
            
        ## how many items/atoms go into each group?
        nAtoms = N.ones(n_centers, 'i') * int(len( a_indices ) / n_centers)
        i=0
        while N.sum(nAtoms) != len( a_indices ):
            nAtoms[i] += 1
            i += 1

        ## distribute atom indices into groups
        result = []
        pos = 0
        for n in nAtoms:
            result += [ N.take( a_indices, N.arange(n) + pos) ]
            pos += n
            
        return result

    
    def nextAtom( self, resName, resNumber, name, chainId, segid):
        self.currentAtom += 1
            
        return {'residue_name':resName, 'name':name, 'type':'ATOM',
                'residue_number':resNumber,
                'serial_number': self.currentAtom,
                'segment_id': segid, 'chain_id':chainId,
                'name_original':name, 'element':'X'}
    

    def makeMap( self, maxPerCenter=4 ):
        """
        Calculate mapping between complete and reduced atom list.
        -> (list of lists of int, list of atom dictionaries)
           groups of atom indices into original model, new center atoms
        """

        resIndex = self.m_sorted.resIndex()
        resModels= self.m_sorted.resModels()
        m = self.m_sorted

        self.currentAtom = 0

        groups = []
        atoms = []

        for i in range( len( resIndex ) ):

            first_atom = resIndex[ i ]
            
            if i < len( resIndex )-1:
                last_atom  = resIndex[ i+1 ] - 1
            else:
                last_atom = len( self.a_indices ) - 1
            
            res_name = m.atoms[ first_atom ]['residue_name']
            segid    = m.atoms[ first_atom ]['segment_id']
            chainId  = m.atoms[ first_atom ]['chain_id']
            res_number = m.atoms[first_atom]['serial_number']
            
            ## position of this residue's atoms in original PDBModel (unsorted)
            a_indices = self.a_indices[ first_atom : last_atom+1 ]

            ## for each center create list of atom indices and a center atom
            if res_name != 'GLY' and res_name != 'ALA':

                bb_a_indices = N.compress( resModels[i].maskBB(), a_indices)
                sc_a_indices = N.compress(
                    N.logical_not( resModels[i].maskBB()), a_indices )

                sc_groups = self.group( sc_a_indices, maxPerCenter )

            else:
                bb_a_indices = a_indices
                sc_groups = []

            groups += [ bb_a_indices ]
            atoms  += [ self.nextAtom( res_name, res_number, 'BB',
                                       chainId, segid ) ]

            i = 0
            for g in sc_groups:
                groups += [ g ]
                atoms  += [ self.nextAtom( res_name, res_number, 'SC%i'%i,
                                           chainId, segid) ]
                i += 1

        self.groups = groups
        self.atoms = atoms


    def reduceXyz( self, xyz, axis=0 ):
        """
        Reduce the number of atoms in the given coordinate set. The set must
        have the same length and order as the reference model. It may have
        an additional (time) dimension as first axis.
        xyz - array, (N_atoms x 3) or (N_frames x N_atoms x 3)
        axis- axis with atoms
        -> array (N_less_atoms x 3) or (N_frames x N_less_atoms x 3)
        """
        masses = self.m.atomProfile('mass')
        r_xyz = None

        for atom_indices in self.groups:

            x = N.take( xyz, atom_indices, axis )
            m = N.take( masses, atom_indices )

            center = N.sum( x * N.transpose([m,]), axis=axis) / N.sum( m )

            if axis == 0:
                center = center[N.NewAxis, :]

            if axis == 1:
                center = center[:, N.NewAxis, :]

            if r_xyz == None:
                r_xyz = center

            else:
                r_xyz = N.concatenate( (r_xyz, center), axis )

        return r_xyz


    def reduceToModel( self, xyz=None, reduce_profiles=1  ):
        """
        Create a reduced PDBModel from coordinates. Atom profiles the source
        PDBModel are reduced by averaging over the grouped atoms.
        xyz - array (N_atoms x 3) or None (->use reference coordinates)
        -> PDBModel with reduced atom set and profile 'mass'
        """

        mass = self.m.atomProfile('mass')
        xyz = xyz or self.m.getXyz()

        mProf = [ N.sum( N.take( mass, group ) ) for group in self.groups ]
        xyz = self.reduceXyz( xyz )

        result = PDBModel()
        result.setAtoms( self.atoms )
        result.setXyz( xyz )
        result.setAtomProfile( 'mass', mProf )
        
        if reduce_profiles:
            self.reduceAtomProfiles( self.m, result )
            
            result.rProfiles = self.m.rProfiles

        return result


    def reduceAtomProfiles( self, from_model, to_model ):
        """
        reduce all atom profiles according to the calculated map by calculating
        the average over the grouped atoms.
        model - PDBModel
        profname - str
        do_sum - 1|0, use sum of atom values instead of average (default 0)
        """
        for profname in from_model.aProfiles:

            p0 =  from_model.atomProfile(profname)
            info = from_model.profileInfo( profname )

            pr = [ N.average( N.take( p0, group ) ) for group in self.groups ]
            
            to_model.setAtomProfile( profname, pr )
            to_model.setProfileInfo( profname, **info )
        


if __name__ == '__main__':

##     from Numeric import array_constructor

    m = PDBModel( T.testRoot()+'/com_wet/1BGS.pdb' )
    m = m.compress( N.logical_not( m.maskH2O() ) )

    m.setAtomProfile('test', range(len(m)))
    
    red = ReduceCoordinates( m, 4 )

    mred = red.reduceToModel()

    print 'Atoms before reduction %i'%m.lenAtoms()
    print 'Atoms After reduction %i'%mred.lenAtoms()
##     frames = red.reduceXyz( t.frames, axis=1 )

##     ref = red.reduceToModel( t.ref.getXyz() )

##     t.frames = frames
##     t.ref = ref

