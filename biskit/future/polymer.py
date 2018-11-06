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

from biskit.future.model import Model
from biskit.future.profilemirror import ProfileMirror

class Feature:
    """
    A Feature is adding information like, e.g., a structure to part of a 
    Polymer and provides a specialized view on this part of the Polymer.
    About like this::
    
        xxxxxxxxxxxxxxPolymerxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
             * ****  **********Feature1*******
                                           *******Feature2*******
    
    The idea is to map specialized content about parts of a molecule to the 
    complete thing.
    
    The most important fields of a feature are:
    
    * model -- the "parent" Polymer to which the feature belongs

    * map   -- the atom positions in the Polymer that are occupied by the
               feature
               
    * atoms -- a virtual ProfileCollection that mirrors all atom profiles of 
               the underlying parent model into  the atom space covered by
               the feature
               
    That means a feature only pretends to have its own profiles. For example,
    asking for the name of the first atom::
    
        feature['name'][0]
    
    will access the atom profile of the parent model and give the name of the 
    first atom which is covered by the feature. This is equivalent to::
    
        feature.model['name'][ feature.map[0] ]
        
    Ah well, maybe that's a sligthly complicated way of putting it... 
    """
    
    def __init__(self, model, map2model, add2model=True ):
        """
        Create a new feature for a given Polymer.
        
        @param model: existing Polymer
        @type  model: Polymer
        @param map:   map feature positions to positions in the model
        @type  map:   [ int ]
        @param add2model: register the new feature with the parent model [1]
        @type  add2model: bool
        """
        assert isinstance( model, Polymer )

        self.model = model
        self.map   = N.array( map2model, 'i' )
        
        self.atoms = ProfileMirror( self.model.atoms, self.map )

        if add2model:
            self.model.addFeature( self )


    def take( self, i, add2model=True ):
        """
        Extract part of the feature.
        
        @param i: feature indices
        @type  i: [ int ] or N.array of int
        
        @return: a new Feature connected to the same Polymer
        """
        rmap = N.take( self.map, i )
        return self.__class__( self.model, rmap, add2model=add2model )
    

    def concat( self, *others ):
        """
        Concatenate this with one or more other features (of the same 
        Polymer). Example::
        
            feature123 = feature1.concat( feature2, feature3 )
        
        @return: a new Feature
        @rtype : Feature OR sub-class (same class as this Feature)
        """
        if len( others ) == 0:
            return self

        next = others[0]

        assert isinstance( next, self.__class__ ), \
               'can only concat with other Feature instances'
        
        assert next.model is self.model, \
               'can only concat features that belong to the same Polymer.'

        rmap = N.concatenate( (self.map, next.map ) )
        
        r = self.__class__( rmap, self.model )

        return r.concat( *others[1:] )

    
    def __len__( self ):
        return len( self.map )
    
    def modelMask( self ):
        """
        Create a mask for parent Polymer, marking every position covered by
        this feature.
        @return:
        @rtype: N.array of int or bool
        """
        r = N.zeros( len( self.model ) )
        N.put( r, self.map, 1 )
        return r
        
    
    def feature2modelIndices( self, i ):
        """
        Convert positions in feature into positions in the overall model.
        """
        return N.take( self.map, i )
    

    def model2featureIndices( self, i ):
        """
        Convert model positions into the corresponding feature positions.

        Note:
          * indices that fall outside the feature are filtered out
          
        @param i: indices pointing to atoms in the parent model
        @type  i: [ int ]
        @return: feature positions corresponding to i
        @rtype : N.array of int
        """
        ## filter out positions that do not point to any feature atom
        i = N.compress( N.take( self.modelMask(), i ), i )
        
        ## shortcut
        if len(i) == 0:
            return i
        
        r = N.zeros( len( self.model ), int)

        ## map from each model position to the corresponding feature position
        ## I first tried with N.put but it doesn't work the supposed way...
        for j in range(len(self.map)):
            r[ self.map[j] ] = j

        ## extract the 
        return N.take( r, i )

    
    def _projectMap( self, model_i ):
        """
        Create a new feature -> model map for the given model positions.
        see: transfer()
        """
        ## positions in new model that still have a corresponding feature atom
        mask_f2new_model = N.take( self.modelMask(), model_i )

        ## the new feature2model map
        return N.flatnonzero( mask_f2new_model )


    def transfer( self, new_model, model_i, add2model=True ):
        """
        Transfer a copy of this feature into a model that is about to be
        extracted from the current parent model. (This method is needed for
        Polymer.take.) The method calculates the new feature2model map
        from the current map and the given model indices.

        @param new_model: new parent model derrived from the current one
        @type  new_model: Polymer
        @param model_i  : positions of the new model's atoms in the old model
        @type  model_i  : [ int ]
        @param add2model: register the Feature with the new parent model [1]
        @type  add2model: bool

        @return: Extract of this feature connected to the new model
        @rtype:  Feature or subclass thereof
        """
        f_new = self.take( self.model2featureIndices( model_i ), add2model=0 )

        f_new.map = self._projectMap( model_i )

        f_new.model = new_model
        
        if add2model:
            new_model.addFeature( f_new )

        return f_new
        
    def cloneTo( self, model=None, map2model=None, add2model=False ):
        """
        Clone this Feature to another (or the same) parent model. This method
        differs from transfer in that it needs a pre-computed model2feature
        map.
        
        @param model: new parent model, if different from current [None]
        @type  model: Polymer
        @param map2model: new feature->model map, if different from old [None]
        @param map2model: [int]
        @param add2model: register new feature with parent model [0]
        @type  add2model: bool
        
        @return: an independent copy of this feature pointing to the same model
        @rtype : Feature or subclass thereof
        """
        model = model or self.model
        if map2model is None: map2model = self.map
        
        return self.__class__( model, map2model, add2model=add2model )


class Polymer( Model ):
    
    #: default profiles for atoms
    ATOM_KEYS    = ['residue_name', 'name', 'element']

    #: default profiles for residues
    RESIDUE_KEYS = ['name']

    #: default profiles for chains
    CHAIN_KEYS = []
    
    def __init__(self, *sequences ):
        """
        Create a new Polymer.

        @param sequence: sequence, 1-letter coded or list of 3-letter names; 
                         several sequences are interpreted as distinct chains
        @type  sequences: str or [ str ]
        """
        Model.__init__( self )
        
        #: list of features mapping into this model
        self.features = []
        
        if sequences:
            self.addSequences( *sequences )
            
        self.__version__ = B.__version__

    
    def __newChain( self ):
        """
        Prepare adding a new chain/molecule to this model.
        """

        if self._chainIndex is None:
            self._chainIndex = N.array( [0] )
        else:
            self._chainIndex = N.concatenate( (self._chainIndex, 
                                               [ self.lenAtoms() ] ) )


    def __extendProfile( self, collection, key, n_items ):
        """
        Add n_items default values to profile key in collection.
        """
        collection[key].extend( [collection[key, 'default']] * n_items )


    def addSequence( self, sequence, newchain=True, default=None ):
        """
        Add sequence to this model with or without creating a new chain.
        For each residue, the necessary atoms are created from a standard 
        dictionary with standard names.
        
        @param sequence: sequence, 1-letter coded or list of 3-letter names
        @type  sequences: str or [ str ]
        @param newchain: put added residues into a new chain [True]
        @type  newchain: bool
        @param default: default value to add to existing profiles [None]
        @type  default: any
        """
        if newchain and len( sequence ) > 0:
            self.__newChain()

        rrnames = self.residues['name']    # residue names for residue profile 
        anames  = self.atoms['name']         # atom names for atom profile
        rnames  = self.atoms['residue_name'] # residue names for atom profile
        elements= self.atoms['element']      # element codes for atom profile

        resIndex= self._resIndex
        i_res = 0                      # pointer to current residue position
        if resIndex is None and len(sequence) > 0:
            resIndex = N.array( [], N.int )
        else:
            i_res = self.lenAtoms()
   
        if type( sequence ) is str:
                sequence = MU.single2longAA( sequence )
        
        atoms = []

        for res in sequence:
            i_res = i_res + len(atoms) ## number of atoms found in previous it.
            resIndex = N.concatenate( (resIndex, [ i_res ]) )

            rrnames.append( res )
        
            atoms = MU.atomDic[ res ][:-1]
            anames.extend( atoms )  ## include all atoms but last OXT
            rnames.extend( [res] * len( atoms ) )  ## append residue_names
        
            ## too simple minded but there is no good dictionary in molUtlis
            elements.extend( [ a[0] for a in atoms ] )

        ## add terminal OXT
        anames.append( MU.atomDic[res][-1] )
        rnames.append( res )
        elements.append( MU.atomDic[res][-1][0] )
        
        self._resIndex = resIndex


    def addSequences( self, *sequences ):
        """
        m.addSequence( seq_1, [ seq_2, ..., seq_n] )
        Append sequence(s) for one or more chains of residues to this model.
        
        For each residue, the necessary atoms are created from a standard 
        dictionary with standard names.
        @param sequence: sequence, 1-letter coded or list of 3-letter names; 
                         several sequences are interpreted as distinct chains
        @type  sequences: str or [ str ]
        """
        for seq in sequences:
            self.addSequence( seq )



    def addFeature( self, feature ):
        self.features.append( feature )
    

    def take( self, i ):
        """
        Extract a Polymer with a subset of atoms::
          take( atomIndices ) -> Polymer / sub-class.

        @param i: atomIndices, positions to take in the order to take
        @type  i: list/array of int

        @return: Polymer / sub-class
        @rtype: Polymer
        """
        r = Model.take( self, i )
        
        ## extract features...
        r.features = [ f.transfer( r, i ) for f in self.features ]
        
        return r

    def _concatHook( self, result, nextModel ):
                        
        for f in self.features:
            result.addFeature( f.cloneTo( result, add2model=False ) )

        for f in nextModel.features:
            fmap = f.map + self.lenAtoms()
            result.addFeature( f.cloneTo( result, fmap, add2model=False) )



    def addAtoms( self, atoms ):
        pass
    

    def sequence(self, xtable=MU.xxDic ):
        """
        Amino acid sequence in one letter code.

        @param xtable: dict {str:str}, additional residue:single_letter mapping
                       for non-standard residues (default molUtils.xxDic)
        @type  xtable: dict

        @return: 1-letter-code AA sequence (based on first atom of each res).
        @rtype: str
        """

        l = self.residues['name']
        return ''.join( MU.singleAA( l, xtable ))


if __name__ == '__main__':
    
    m = Polymer('ACAGPL','SS')

    f_ = Feature(m, list(range( 1, len(m), 2)) )

    m_ = m.take(range( 1, len(m), 2)) 
    
    m2 = m.concat( m )
    
    assert m2.features[0].atoms['name'] == m2.features[1].atoms['name']
    
    m3 = m2.take(range( 0, 52, 2))
    
    views = m.atoms[:]
    views2 = m2.atoms[:]
    views3 = m3.atoms[:]
    
    views = m2.atoms.take(range( 0, 52, 2)).toDicts()
    assert views == m3.atoms[:]

    
