##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2009-05-09 14:17:28 +0200 (Sat, 09 May 2009) $
## $Revision: 742 $
"""
Attempt at a residue view to connect to PDBModel
"""
import numpy as N
import weakref 

class ResidueFactory( object ):
    """
    Produce residue instance of a certain class given a residue name
    """

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
    
    #: residue is marked as hetatom
    hetatom = False

    #: residue is polar
    polar = False
    
    acceptors = ['O', 'C']
    donors    = ['H']    
 
    def __init__( self, model=None, index=None, terminal=False ):
        """
        @param model: PDBModel, parent model of this residue
        @param index: int, position of this residue in model residue index
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
        """Map atom name to position *within* residue, starting with 0"""
        if self._atom_index is None:

            self._atom_index = {}
            m = self.model() or self.ref
        
            _from = self.from_atom
            _to = self.to_atom

            self._atom_index = dict( zip( m.atoms['name'][_from : _to],
                                     range( _to - _from ) ) )
        
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
        assert self.model() is None, 'cannot modify non-attached residue'
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
        assert self.model() is None, 'cannot modify non-attached residue'
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
        @return: [str], list of standard atom names in standard order
        """
        assert self.ref is not None, 'no reference residue defined'
        return self.ref.atomNames()

########### TEST

from Biskit import *

m = PDBModel( '3TGI' )
r = Residue( m, 3 )

print r.from_atom, r.to_atom
print r.chain
print r.segid
print r.xyz
print r.serial


    