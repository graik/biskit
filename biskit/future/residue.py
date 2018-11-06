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
Attempt at a residue view to connect to PDBModel
"""
import biskit.tools as T

import numpy as N
import weakref 
import string
import glob

def nextKey( lst, key ):
    """
    Fetch next free key in a list or dict. If 'key' already exists, return
    'key_1', if 'key_1' already exists, return 'key_2' etc.
    """
    if not key in lst:
        return key

    l = key.split( '_' )
    base = l[0]
    if len( l ) < 2:
        return nextKey( lst, '%s_%i' % (base, 1) )  ## start recursion

    ext = int( l[1] ) + 1
    return nextKey( lst, '%s_%i' % (base, ext) )  ## recursion until first free


class ResidueFactory( object ):
    """
    Produce residue instance of a certain class given a residue name
    """

    RESFOLDER = 'residues'

    @staticmethod
    def get( name ):
        """
        Emit residue instance of correct class from name
        """
        pass
    
    def newAmino( name, code=None, letter=None, atoms=[], polar=False, 
                  acceptors=[] ):
        pass

    def newNucleic( name, code=None, letter=None, atoms=[], polar=False, 
                  acceptors=[] ):
        pass


class Atom(object):
    
    name = 'CA'

    standard_charge = 0.
    element = 'C'

    fixed_hydrogens = ['H']
    variable_hydrogens = []
    
    synonymes = []
    
    def __init__( self, name, xyz=None, **kw ):
        self.name = name
        self.xyz  = xyz
        self.__dict__.update( kw )
    
    def setXyz( self, xyz ):
        self.xyz = xyz
        
    def getXyz( self ):
        return self.xyz
    
    xyz = property( getXyz, setXyz )
    
def returnNone():
    """Helper function simulating an empty weakref"""
    return None
    

class Residue(object):
    """
    Instance of an actual residue in a structure. Class variables define 
    global properties of this type of residue. Object properties / variables 
    keep track of data that are specific for the individual residue instance.
    """
    ## pre-defined residue types
    OTHER = 0
    AMINO = 1
    DNA = 2
    RNA = 3
    SUGAR = 4
    ION = 5
    SOLVENT = 6
  
    ## global residue-specific properties 
    
    name  = 'alanine'
    code  = 'ala'
    letter= 'A'

    type = AMINO
    
    code_synonymes   = []  #: [ str ], list of alternative 3-letter names
    letter_synonymes = []  #: [ str ], list of alternative 1-letter codes
    
    #: reference PDBModel for this residue; 
    #: defines standard atom order and content, and reference coordinates
    ref = None

    #: element codes for standard atoms
    elements = {'N':'N', 'CA':'C', 'C':'C', 'O':'O', 'CB':'C', 
                'H':'H', 'HA':'HA' }
    
    #: residue should be marked as hetatom
    hetatom = False

    #: residue is polar
    polar = False
    
    acceptors = ['O', 'C']
    donors    = ['H']    
 
    def __init__( self, model=None, index=None, terminal=False ):
        """
        :param model: PDBModel, parent model of this residue
        :param index: int, position of this residue in model's residue index
        """
        self.model = weakref.ref( model ) if model else returnNone
        self.i = index or 0
        self.terminal = terminal

        self.reset()

    def reset(self):
        """Reset all cached variables"""
        self._from_atom = None
        self._to_atom = None
        self._atom_index = None
        
    def _getAtomIndex( self ):
        """
        Map atom name to position *within* residue, starting with 0. This
        creates a dictionary indexed by atom names and pointing to the position
        of each atom within the residue.
        If the same atom name occurs twice or more (due to alternate positions),
        the second and third occurrence are indexed as 'name_1' and 'name_2'
        respectively.
        :return: dict connecting atom names to their position within residue
        :rtype:  {str:int}
        """
        if self._atom_index is None:

            self._atom_index = {}
            m = self.model() or self.ref
        
            _from = self.from_atom
            _to = self.to_atom
            
            nr_names = []
            for n in m.atoms['name'][_from:_to] :
                nr_names.append( nextKey( nr_names, n ) )

            self._atom_index = dict( zip( nr_names, range( _to - _from ) ) )
        
        return self._atom_index

    atom_index = property( _getAtomIndex, 
                      doc="dict, atoms indexed by name (readonly, cached)" )
    
    def _getFromAtom( self ):
        """Fetch and cache starting position of this residue in parent model"""
        if self._from_atom is None:
            if self.model() is None:
                self._from_atom = 0
            else:
                self._from_atom = self.model().resIndex()[self.i]
        return self._from_atom

    from_atom = property( _getFromAtom, 
            doc="Index of first residue atom within model (read-only, cached)")
    
    def _getToAtom( self ):
        """
        Fetch and cache ending position of this residue in parent model.
        The ending position corresponds to the first atom of the next residue
        """
        if self._to_atom is None:
            m = self.model() or self.ref
         
            rI = m.resIndex()
    
            if self.i + 1 < len( rI ):
                self._to_atom = rI[ self.i + 1 ]  # first atom of next residue
            else:
                self._to_atom = len( m )

        return self._to_atom
    
    to_atom = property( _getToAtom, 
                        doc="First atom of next residue (read-only, cached)" )

    
    def __len__( self ):
        return self.to_atom - self.from_atom
    
    def __letter2number( self, l ):
        try:
            return int( l )
        except ValueError:
            return string.ascii_letters.find( l.lower() )
    
    def __getitem__( self, key ):
        """
        res[ 'CA' ] ==> first atom with name CA
        res[ 'CA', 1 ] ==> second atom with name CA
        res[ 'CA', A ] ==> first atom with name CA
        res[ 'CA', B ] ==> second atom with name CA
        
        res['temperature_factor'] ==> profile values for residue atoms
        """
        m = self.model() or self.ref

        if type( key ) is str:
            ## return profile values for residue or residue atoms
            if key in m.atoms:
                return m.atoms[key][self.from_atom : self.to_atom]
            if key in m.residues:
                return m.residues[key][self.i]
            
            ## return cross view of single atom of name 'key'
            if key in self.atom_index:
                return m[ self.atom_index[ key ] ]
        
        if type( key ) is tuple and len( key ) == 2:
            ## interpret key as pair of atom name and alternate code
            name, alt = key
            alt = self.__letter2number( alt )
            
            try:
                if alt == 0:
                    return m[ self.atom_index[ name ] ]
                
                return m[ self.atom_index[ '%s_%i' % (name, alt) ] ]

            except IndexError as why:
                raise IndexError('atom not found: %r (%r)' % (name, alt))
        
        raise IndexError('No profile or atom name matches key ' + str( key ))
        
    def mass( self ):
        pass
    

    def _getIndex( self ):
        """Get position of this residue within the model's residue list"""
        return self._index
    
    index = property( _getIndex )  #: read-only
    

    def _getSerial( self ):
        """Get residue number within PDB file"""
        if self.model() is None:
            return 1
        return self.model().atoms['residue_number'][ self.from_atom ]

    def _setSerial( self, value ):
        """Modify residue number within PDB file (doesn't move residue!)"""
        assert self.model() is not None, 'cannot modify non-attached residue'
        self.model().atoms['residue_number'][self.from_atom : self.to_atom] \
            = [ value ] * self.to_atom - self.from_atom
    
    serial = property( _getSerial, _setSerial )
    
 
    def _getChain( self ):
        """Fetch chain ID for this residue from underlying model"""
        if self.model() is None:
            return ''
        return self.model().atoms['chain_id'][ self.from_atom ]
    
    def _setChain( self, letter ):
        """Modify chain ID of underlying model at this residue position"""
        assert len( letter ) <= 1
        assert self.model() is not None, 'cannot modify non-attached residue'
        self.model().atoms['chain_id'][self.from_atom : self.to_atom] \
            = [ letter ] * self.to_atom - self.from_atom
    
    chain = property( _getChain, _setChain )


    def _getSegid( self ):
        """Fetch segment ID for this residue from underlying model"""
        if self.model() is None:
            return ''
        return self.model().atoms['segment_id'][ self.from_atom ]
    
    def _setSegid( self, segid ):
        """Modify 4-letter segment ID of underlying model at residue's position"""
        assert len( segid ) <= 4, 'segment id must be 4 letters or less'
        assert self.model() is None, 'cannot modify non-attached residue'
        self.model().atoms['segment_id'][self.from_atom : self.to_atom] \
            = [ segid ] * self.to_atom - self.from_atom
    
    segid = property( _getSegid, _setSegid )
    

    def _getXyz( self ):
        """Fetch xyz array for this residue's atoms from underlying model"""
        m = self.model() or self.ref
        return m.getXyz()[ self.from_atom : self.to_atom ]
    
    def _setXyz( self, xyz ):
        """Modify xyz array of underlying model at this residue position"""
        assert N.shape(xyz) == (len(self.atom_index), 3)
        assert self.model() is None, 'cannot modify non-attached residue'
        self.model().getXyz()[self.from_atom : self.to_atom] = xyz

    xyz = property( _getXyz, _setXyz )
    
    
    def standardAtoms(self):
        """
        :return: [str], list of standard atom names in standard order
        """
        assert self.ref is not None, 'no reference residue defined'
        return self.ref.atomNames()


    def labelDuplicateAtoms( self ):
        """
        Ensure duplicate atoms are labeled with an alternate code. Atoms
        that have no duplicates receive alternate code ''. Atoms with one or
        more duplicates within the same residue receive alternate codes
        'A' (first one in list), 'B', and so on.
        The original alternate code field from the PDB (if any) is deleted.
        """
        import string
        
        m = self.model()
        assert m is not None, 'residue is not attached to any PDBModel'
        
        names = m.atoms['name'][self.from_atom : self.to_atom]

        ## number of occurrences of each atom name
        counts = [ names.count( n ) for n in names ]
        if max( counts ) <= 1:
            return  ## nothing to do
        
        ## dict with list of positions in which each atom is found
        positions = {}
        for i,n in enumerate(names):
            T.dictAdd( positions, n, i, forceList=True )
        
        ## create new list of alternate codes
        alt = []
        for n, c in zip( names, counts ):
            pos = positions[ n ]
            
            if c == 1:
                alt.append( '' )
            else:
                alt.append( string.ascii_letters[ c - len( pos ) ].upper() )

            pos.pop()  # remove one, discard

        m.atoms['alternate'][self.from_atom : self.to_atom] = alt
    
        
#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    M = None

    def prepare( self ):
        import biskit as B
        import biskit.tools as T

        self.M = self.M or T.load( T.testRoot() + '/lig/1A19_dry.model' )
    
    def test_atomIndex( self ):
        """ProfileCollection test"""
        m = self.M.clone()
        r = Residue( m, 0 )

        self.assertEqual( len( r.atom_index ), 15 )
        
        m.mergeResidues( 0 )
        r.reset()
        
        self.assertEqual( len( r.atom_index ), 28 )
        self.assertTrue( 'CA_1' in r.atom_index )

    def test_invalidIndex( self ):
        m = self.M.clone()
        
        r = Residue( m, m.lenResidues )
        
if __name__ == '__main__':

    BT.localTest()
    
    from biskit import *
    import biskit.tools as T
    
    m = PDBModel( T.testRoot() + '/lig/1A19_dry.model'  )
    r = Residue( m, 0 )
    
    print(r.from_atom, r.to_atom)
    print(r.chain)
    print(r.segid)
    print(r.xyz)
    print(r.serial)
    
    m.mergeResidues( 0 )
    r.reset()

    
