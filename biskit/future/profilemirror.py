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

from biskit.profileCollection import _ViewSignal, CrossView, ProfileError

class ProfileMirror( B.ProfileCollection ):
    """
    Access only part of an underlying ProfileCollection from a Polymer.
    """
    
    def __init__( self, parentProfiles, map2model ):
        """
        @param parentProfiles: profiles of the Polymer
        @type  parentProfiles: ProfileCollection
        @param map: map positions in this ProfileCollection to positions in the
                    parent profile
        @type  map: N.array of int
        """
        #: feature -> model map from parent feature
        self.map = map2model

        #: a ProfileCollection of the Polymer that is parent of the Feature
        assert isinstance( parentProfiles, B.ProfileCollection )
        self.pc = parentProfiles

        #: default value to use if profile is removed / changed by set()
        self.default = None

        #: support CrossView de-activation
        self._viewSignal = _ViewSignal()
    
        self.initVersion = B.__version__

 
    def __len__( self ):
        """
        @return: number of profiles in the collection
        @rtype: int        
        """
        return len( self.pc )

    def keys( self ):
        return self.pc.keys()

    def has_key( self, k ):
        return k in self.pc  ## verify !!

    def __iter__(self):
        """
        Iterate over profile keys::
          for k in self  <==>  for k in p.keys()
        
        @return: list of items
        @rtype: list
        """
        return iter(self.pc)


    def getInfo( self, name ):
        """
        Use::
           getInfo( name ) -> dict with meta infos about profile::
           
        Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1|0

        @param name: profile name
        @type  name: str

        @return: dict with infos about profile
        @rtype: dict
        @raise ProfileError: if no profile is found with |name|
        """
        return self.pc.getInfo( name )


    def get( self, name, default=None ):
        
        r = self.pc.get( name, default=None )
        
        if r is None:
            r = default

        if self.pc[ name, 'isarray']:
            return N.take( r, self.map )
        
        return [ r[i] for i in self.map ]


    def setInfo( self, name, **args ):
        """
        Add/Override infos about a given profile::
          e.g. setInfo('relASA', comment='new', params={'bin':'whatif'})
        
        @raise  ProfileError: if no profile is found with |name|
        """
        self.pc.setInfo( name, **args )

    def __iszero( self, r):
        """
        checks whether given r is set to the integer 0 without being confused if r is an array
        """
        return type(r) is int and r == 0


    def set( self, name, prof, mask=None, default=None, asarray=1,
             comment=None, **moreInfo ):
        """
        Add/override a profile. The two info records 'version', 'changed' and 
        'isarray' are always modified but can be overridden by key=value pairs 
        to this function. If the profile does not yet exist in the parent 
        ProfileCollection, a new profile is created for the whole Polymer
        and the values not covered by this ProfileMirror are set to <default>.
        
        @param name: profile name (i.e. key)
        @type  name: str
        @param prof: list of values OR None
        @type  prof: [any] OR None
        @param mask: list 1 x N_items of 0|1, if there are less values than
                      items, provide mask with 0 for missing values,
                      N.sum(mask)==N_items
        @type  mask: [int]
        @param default: value for items masked and for new parent profile
                        (default: None for lists, 0 for arrays]
        @type  default: any
        @param asarray: store as list (0), as array (2) or store numbers as
                        array but everything else as list (1) (default: 1)
        @type  asarray: 0|1|2
        @param comment: goes into info[name]['comment']
        @type  comment: str
        @param moreInfo: additional key-value pairs for info[name]
        @type  moreInfo: key=value
        
        @raise ProfileError: if length of prof != length of other profiles
        @raise ProfileError: if mask is given but N.sum(mask) != len(prof)
        """
       
        ## consistency check
        if mask and N.sum(mask) != len(prof):
            raise ProfileError(
                "Mask doesn't match profile ( N.sum(mask)!=len(prof) ). " +
                "%i != %i" % (N.sum(mask), len( prof ) ) )

        r = self.pc.get( name, default=0 )
        
        ## take array status from existing profile
        if not self.__iszero(r):
            if self.pc.get( (name, 'isarray') ):
                asarray = 2
            else:
                asarray = 0

        prof = self.array_or_list( prof, asarray )

        ## use default == 0 for arrays
        if not default and isinstance( prof, N.ndarray ):
            default = 0

        self.default = default
        
        ## expand profile to have a value also for masked positions
        prof = self.expand( prof, mask, default )

        l = len( self.map )
        if l and len( prof ) != l:
            raise ProfileError( "Profile %s has wrong length." % name )
        
        
        ## create a new profile
        if self.__iszero(r):
            r = [ default ] * self.pc.profLength()
            if isinstance( prof, N.ndarray ):
                r = N.array( r )
        
        if isinstance( r, N.ndarray ):
            N.put( r, self.map, prof )
        else:
            r = N.array( r, 'O' )
            N.put( r, self.map, prof )
            r = r.tolist()
            
        self.pc.set( name, r, asarray=asarray, comment=comment, **moreInfo )
        
        
    def take( self,indices, *initArgs, **initKw ):
        """
        Take from profiles using provided indices::
          take( indices ) -> FeatureProfiles mapping to subset of parent

        @param indices: list of indices into this profile (not the Polymer)
        @type  indices: [int]

        @return: new instance that maps to a new subset of positions
        @rtype: ProfileMirror
        """
        rmap = N.take( self.map, indices )
        return self.__class__( rmap, self.pc )

    def concat( self, *profiles ):
        """
        Concatenate this with the given ProfileMirror(s).::
          p0.concat( p1 [, p2, ..]) -> single ProfileMirror with the
          same number of profiles as p0 but with the length of p0+p1+p2..

        Note: it's not yet clear what meaning this will have for ProfileMirrors

        @param profiles: profile(s) to concatenate
        @type  profiles: ProfileCollection(s)
        
        @return: concatenated profile(s)  
        @rtype: ProfileCollection / subclass
        """
        if len( profiles ) == 0:
            return self

        next = profiles[0]

        assert isinstance( next, self.__class__ ), \
               'can only concat with other ProfileMirrors'

        rmap = N.concatenate( (self.map, next.map ) )
        
        r = self.__class__( rmap, self.pc )

        return r.concat( *profiles[1:] )


    def remove( self, *key ):
        """
        Profiles can only be removed by the parent ProfileCollection
        @raises ProfileError: always
        """
        raise B.ProfileError('Cannot remove profile %r from a ProfileMirror'\
              % key)

    def clone( self ):
        """
        Clone profile mirror::
          clone() -> ProfileMirror (or sub-class)

        @return: profile mirror
        @rtype: ProfileCollection          
        """
        return self.__class__( copy.copy( self.map ), self.pc )
    
    def clear( self ):
        """
        Profiles can only be cleared by the parent ProfileCollection
        """
        raise B.ProfileError('Cannot clear profiles from a ProfileMirror')

    def profLength( self ):
        """
        Length of profile::
          profLength() -> int; length of first non-None profile or 0
        
        @return: length of first non-None profile or 0
        @rtype: int     
        """
        return len( self.map )
    


if __name__ == '__main__':

    import string

    p = B.ProfileCollection()
    p['name'] = string.ascii_letters
    p['id']   = list(range( len(string.ascii_letters)))

    ## mirror looks at every second position
    mirror = ProfileMirror( p, list(range(0, len(string.ascii_letters), 2)) )

    assert mirror['name'] == list( string.ascii_letters[::2] )
    assert N.all( mirror['id'] == list(range( 0, p.profLength(), 2)) )

    ## create a new profile
    mirror['m_id'] = list(range( mirror.profLength()))
    assert N.all( p['m_id'][::2] == mirror['m_id'] )

    mirror['name'][2] = '#'  ## does not have any effect
    mirror['name'] = ['#'] * mirror.profLength()
    assert p['name'][::2] == ['#'] * mirror.profLength()
    assert p['name'][1::2]== list( string.ascii_letters[1::2] )
    
