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
"""Work in progress"""

import copy
import numpy as N

import biskit as B
import biskit.tools as T
import biskit.molUtils as MU

def index2map( index, len_i ):
    """
    For example
    index2map([3,5,10], 12) ==> [0,0,0, 1,1, 2,2,2,2,2, 3,3,3]
    
    @param index: list of starting positions, e.g. [0, 3, 8]
    @type  index: [ int ] or N.array of int
    @param len_i: length of target map, e.g. 10
    @type  len_i: int

    @return: list mapping atom positions to residue(/chain) number,
             e.g. [0,0,0, 1,1,1,1,1, 2,2] from above example
    @rtype: N.array of int (and of len_i length)
    """
    index = N.concatenate( (index, [len_i]) )
    delta = index[1:] - index[:-1] 
    return N.repeat( list(range(len(delta))), delta)


def map2index( imap ):
    """
    Identify the starting positions of each residue(/chain) from a map
    giving the residue(/chain) number of each atom.

    @param imap: something like [0,0,0,1,1,1,1,1,2,2,2,...]
    @type  imap: [ int ]

    @return: list of starting positions, e.g. [0, 3, 8, ...] in above ex.
    @rtype: N.array of int
    """
    try:
        imap = N.concatenate( (imap, [imap[-1]] ) )
        delta = imap[1:] - imap[:-1] 
        r = N.flatnonzero( delta ) + 1
        return N.concatenate( ( [0], r ) )

    except IndexError:
        ## handle empty imap parameter
        return N.zeros(0)

class Model( object ):
    """
    Model is intended to become the common base class for PDBModel and
    Polymer. 
    """
    #: [str], default profiles for atoms
    ATOM_KEYS    = []

    #: [str], default profiles for residues
    RESIDUE_KEYS = []

    #: [str], default profiles for chains
    CHAIN_KEYS = []


    def __init__(self):
        """
        Create a new empty Model instance. Don't use this constructor directly.
        """
        #: starting (atom) position of each residue
        self._resIndex = None
        #: starting position of each chain
        self._chainIndex = None

        #: values associated with atoms
        self.atoms    = B.ProfileCollection()
        #: values associated with residues
        self.residues = B.ProfileCollection()
        #: values associated with chains or molecules
        self.chains   = B.ProfileCollection()
        
        for key in self.ATOM_KEYS:
            self.atoms.set( key, [], asarray=False )

        for key in self.RESIDUE_KEYS:
            self.residues.set( key, [], asarray=False )
        
        for key in self.CHAIN_KEYS:
            self.chains.set( key, [], asarray=False )
        
        #: Meta info
        self.info = { 'date':T.dateSortString() }
        self.__version__ = B.__version__
    

    def __len__( self ):
        return self.lenAtoms()
    

    def __getitem__( self, k ):
        """
        Get atom profile or profile item or CrossView for one atom::
          m['prof1']         <==>  m.atoms.get( 'prof1' )         
          m['prof1','info1'] <==>  m.atoms.get( 'prof1','info1' )
          m[10]              <==>  CrossView( m.atoms, 10 )

        @return: profile OR meta infos thereof OR CrossView dict
        @rtype: list OR array OR any OR CrossView
        """
        if type( k ) is str:
            if k in self.atoms:
                return self.atoms.get( k )
            if k in self.residues:
                return self.residues.get( k )
            if k in self.chains:
                return self.chains.get( k )
            if k in self.info:
                return self.info[ k ]

        if type( k ) is tuple:
            return self.profileInfo( k[0] )[ k[1] ]
        
        return self.atoms[k]
    
    
    def __setitem__( self, k, v ):
        """
        Set atom profile or profile item (or meta info)::
          m['prof1'] = range(10)    <==> m.atoms.set( 'prof1', range(10) )
            OR                      <==> m.residues.set( 'prof1', range(10) )
          
          m['prof1','info1]='comment'
                             <==> m.atoms.setInfo('prof1',info1='comment')
            OR               <==> m.residues.setInfo('prof1',info1='comment')
        
          m['version'] = '1.0.0'    <==> m.info['version'] = '1.0.0'
            but only if 'version' already exists in m.info 

        @return: item
        @rtype: any        
        """
        if type( k ) is str:
            if v is not None and len( v ) == self.lenAtoms():
                return self.atoms.set( k, v )
            if v is not None and len( v ) == self.lenResidues():
                return self.residues.set( k, v )
            if v is not None and len( v ) == self.lenChains():
                return self.chains.set( k, v )

            if k in self.atoms:
                return self.atoms.set( k, v )
            if k in self.residues:
                return self.residues.set( k, v )
            if k in self.chains:
                return self.chains.set( k, v )
            if k in self.info:
                self.info[ k ] = v
            raise ProfileError('Value cannot clearly be assigned to either atom or '\
                  + 'residue or chain profiles')
                
        if type( k ) is tuple:
            key, infokey = k
            if key in self.atoms:
                self.atoms[key, infokey] = v
                return
            if key in self.residues:
                self.residues[key, infokey] = v
                return
            self.chains[key, infokey] = v
            return
        
        raise ProfileError('Cannot interpret %r as profile name or profile info record' % k)

    
    def __getslice__( self, *arg ):
        """
        Get list of CrossViews::
          m[0:100:5] <==> [ CrossView(m.atoms, i) for i in range(0,100,5) ]
        """
        return self.atoms.__getslice__( *arg )
    
    def __iter__( self ):
        return self.atoms.iterCrossViews()


    def _concatHook( self, resultModel, nextModel ):
        pass
    
    def concat( self, *models ):
        """
        Concatenate the given models (in given order) to this model.
        
        Note for developers: concat is called recursively on the growing model.
          The typical pattern for overriding this method is hence a bit
          different. See _concatHook()!
          
        @param models: models to concatenate
        @type  models: Model, Model, ...
        @return: resulting model
        @rtype: Model or subclass thereof
        """

        if len( models ) == 0:
            return self

        m = models[0]

        r = self.__class__()

        r.residues = self.residues.concat( m.residues, )
        r.atoms    = self.atoms.concat( m.atoms )

        r._resIndex   = N.concatenate(
            (self._resIndex, m._resIndex + self.lenAtoms())) 
        r._chainIndex =N.concatenate(
            (self._chainIndex, m._chainIndex +self.lenAtoms()))
        
        r.info = copy.deepcopy( self.info )
        
        self._concatHook( r, m )

        return r.concat( *models[1:] )


    def lenAtoms(self):
        """
        @return: number of atoms in this model
        @rtype: int
        """
        return self.atoms.profLength()
    
    def lenResidues( self ):
        """
        @return: number of residues in this model
        @rtype: int
        """
        if self._resIndex is None:
            return 0
        return len( self._resIndex )
    
    def lenChains( self ):
        """
        @return: number of chains in this model
        @rtype: int
        """
        if self._chainIndex is None:
            return 0
        return len( self._chainIndex )


    def resMap(  self ):
        """
        Get list to map from any atom to a continuous residue numbering
        (starting with 0).

        @return: array [00011111122223333..], residue index for each atom
        @rtype:  list of int
        """
        return index2map( self._resIndex, self.lenAtoms() )


    def chainMap( self ):
        """
        Get chain index of each atom.

        @return: array 1 x N_atoms of int, e.g. [000000011111111111122222...]
        @rtype: list of int
        """
        return index2map( self._chainIndex, self.lenAtoms() )
    
    
    def take( self, i ):
        """
        Extract a Model with a subset of atoms::
          take( atomIndices ) -> Polymer / sub-class.

        @param i: atomIndices, positions to take in the order to take
        @type  i: list/array of int

        @return: Model / sub-class
        @rtype: Model
        """
        r = self.__class__()

        r.atoms    = self.atoms.take( i )
        
        ## more tricky: rescue residue borders and extract residue profiles
        new_resmap  = N.take( self.resMap(), i )
        ## Note: this erases ordering information and fails for repeated residues
        ## -- see PDBModel version for fix
        r._resIndex = map2index( new_resmap )

        i_res      = N.take( new_resmap, r._resIndex )
        r.residues = self.residues.take( i_res )

        ## now the same with chains
        new_chainmap  = N.take( self.chainMap(), i )
        ## Note: this erases ordering information and fails for repeated residues
        ## -- see PDBModel version for fix
        r._chainIndex = map2index( new_chainmap )
        
        i_chains = N.take( new_chainmap, r._chainIndex )
        r.chains = self.chains.take( i_chains )

        ## copy non-sequential infos
        r.info = copy.deepcopy( self.info )

        return r



