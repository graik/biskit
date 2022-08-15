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
Manage profiles.
"""

import numpy as N

import biskit
import biskit.tools as T
import biskit.mathUtils as M
from biskit import EHandler
from biskit.hist import density


import copy
## manual replace UserDict: http://python3porting.com/problems.html#replacing-userdict
##from UserDict import DictMixin
from collections.abc import MutableMapping as DictMixin
import weakref

try:
    import biggles
except:
    biggles = 0

class ProfileError(Exception):
    pass

class _ViewSignal:
    pass


class CrossView( DictMixin ):
    """
    CrossView instances provide a dictionary-like read- and writeable
    view on one position in a ProfileCollection.

    Example:
    ========
      Let's assume we have a ProfileCollection p containing 2 profiles
      'name' and 'mass':

        >>> p
        ProfileCollection: 2 profiles of length 5
        name
        {...}
               ['a','b','c','d','e']
        mass
        {...}
               [10, 20, 30, 40, 50]

      A CrossView on the second position of p will then look like this:

        >>> view = p[1]
        >>> view
        CrossView{'name': 'b', 'mass': 20}

      And can be used to change the values of this position in the collection:

        >>> view['name'] = 'X'

      ... is equivalent to (but slower than):

        >>> p['name'][1] = 'X'

      The view can also be used to change several profiles simultaneously:

        >>> view.update( {'name':'X', 'mass':100} )

      It can I{not} be used to delete a profile position or remove a key.

    Performance issues:
    ===================
      Creating and using a CrossView to access a single position in a
      profile is not much more costly than the direct
      p['profilename'][index] access -- Both take about 10 times longer
      than accessing a normal dictionary.
      
      However, excessive use of CrossViews can deliver a blow to your
      program's performance. Iterating over a whole ProfileCollection
      in form of CrossViews (:class:`ProfileCollection.iterCrossViews`) is
      about a factor 100 slower than direct iteration over the
      elements of a single profile (array or list). Alternatively,
      :class:`ProfileCollection.iterDicts` yields normal disconnected
      dictionaries. It is in itself not much faster (30%) than
      iterCrossViews but subsequent work with the dictionaries is, at
      least, 10 times more efficient.

    Developer Note:
    ===============
      As a safety measure, CrossView instances are / should be
      invalidated as soon as the underlying ProfileCollection is
      re-ordered or removed. To this end, the CrossView only keeps a
      weak reference to its parent collection (so that the collection
      will not be kept in memory for referencing views alone). Another weak
      reference to the collection's _viewSignal field is used to sense
      other changes -- overriding this field signals the killing
      (alive=False) of all CrossViews pointing to this collection, see
      :class:`ProfileCollection.killViews`.
    """

    def _cease( self, ref ):
        try:
            self.alive = False
        except:
            EHandler.warning('error in CrossView._cease')
            pass

    def __init__( self, parent, index ):
        """
        :param parent: ProfileCollection
        :type  parent: ProfileCollection
        :param index: position in profile collection
        :type  index: int
        """
        self.parent = weakref.ref( parent, self._cease )
        self.index  = index

        self._signal = weakref.ref( parent._viewSignal, self._cease )

        #: parent collection is alive and position is still correct
        self.alive = True

    def __getitem__( self, key ):
        assert self.alive, 'view on dead or changed profile collection'

        return self.parent().profiles[ key ][self.index]

    def __setitem__( self, key, v):
        assert self.alive, 'view on dead or changed profile collection'

        self.parent().profiles[ key ][self.index] = v
        
    def __delitem__(self, key):
        raise ProfileError('Cannot delete profile key %s (or any other) from '\
                           +' CrossView instance.' % key)
    
    def __len__(self):
        return len(self.parent().profiles.keys())

    def keys( self ):
        assert self.alive, 'view on dead or changed profile collection'

        return self.parent().profiles.keys()

##   def has_key(self, key ):
##       assert self.alive, 'view on dead or changed profile collection'
##
##        return key in self.parent().profiles

    def __contains__( self, key ):
        assert self.alive, 'view on dead or changed profile collection'

        return key in self.parent().profiles

    def __iter__( self ):
        for k in self.parent().profiles.keys():
            assert self.alive, 'view on dead or changed profile collection'

            yield k
            
    def iteritems( self ):
        """
        Iterate over key : value pairs. Keys are all keys from the parent
        ProfileCollection, values are the corresponding values at the
        position for wich this CrossView was created.

        :return: key:value pairs
        :rtype: iterator over [ ( str, any ) ]
        """
        for k in self.parent().profiles.keys():
            assert self.alive, 'view on dead or changed profile collection'

            yield (k, self.parent()[ k ][self.index])

    def toDict( self ):
        """
        Convert this CrossView into a standard dictionary that is detached
        from the underlying ProfileCollection.

        :return: dict with values indexed by parent ProfileCollection's keys
        :rtype: dict
        """
        keys = self.parent().profiles.keys()

        values = [ self.parent().profiles[k][self.index] for k in keys ]

        assert self.alive, 'view on dead or changed profile collection'

        return dict( zip( keys, values ) )

    def __repr__(self):
        if not self.alive:
            return 'Dead CrossView{ }'
        return 'CrossView' + DictMixin.__repr__(self)

    def iterkeys(self):
        return self.__iter__()


class ProfileCollection:
    """    
    Manage profiles (arrays or lists of values) for Trajectory frames or
    atoms/residues in PDBModel. ProfileCollection resembles a 2-dimensional
    array where the first axis (let's say row) is accessed by a string key
    rather than an index. Each row has an additional (meta)info dictionary
    assigned to it. The take() and concat() methods operate on the columns,
    i.e. they are applied to all profiles simultaneously.

    By default, profiles of numbers are stored and returned as Numpy.array and
    all others are stored and returned as ordinary list. This behaviour can be
    modified with the option asarray of ProfileCollection.set(). Using both
    lists and arrays is a compromise between the efficiency of Numeric arrays
    and two problems of of the old Numeric module -- (1) arrays of objects
    could not be unpickled (Numeric bug) and (2) arrays of strings would end
    up as 2-D arrays of char. The 'isarray' entry of a profile's info
    dictionary tells whether the profile is stored as array or as list. (Since
    we have now replaced Numeric by numpy, we can probably switch to the
    exclusive use of numpy arrays.)

    Acessing profiles
    =================

      ProfileCollection p can be used like a dictionary of lists::
        len( p )          -> number of profiles (== len( p.profiles ) )
        p['prof1']        -> list/array with values of profile 'prof1'
        del p['prof1']    -> remove a profile
        p['prof1'] = [..] -> add/override a profile without additional infos

        'prof1' in p      -> True, if collection contains key 'prof1'
        for k in p:       -> iterate over profile keys
        for k in p.iteritems():  -> iterate over key:profile pairs

    Accessing Metainfo
    ==================

      Each profile key also has a dictionary of meta infos assigned to it (see
      getInfo(), setInfo(), p.infos). These can be accessed like::
        p['prof1','date']   -> date of creation of profile named 'prof1'
        p.getInfo('prof1')  -> returns all info records
        p['prof1','comment'] = 'first prof'  -> add/change single info value

    CrossViews
    ==========

      ProfileCollections can also be viewed from the side (along columns) --
      :class:`CrossView`s provide a virtual dictionary with the values of all
      profiles at a certain position. Changes to the dictionary will change
      the value in the underlying profile and vice versa, for example::

        atom = p[10]        -> CrossView{'prof1' : 42.0, 'name' : 'CA', ... }
        atom['prof1'] = 33  -> same as p['prof1'][10] = 33
        p[10]['prof1']= 33  -> still the same but much slower, use instead...
        p['prof1'][10]= 33  -> doesn't invoke CrossView and is thus faster

        p[0] == p[-1]  -> True if the values of all profiles are identical at
                          first and last position

        for a in p.iterCrossViews():  -> iterates over CrossView dictionaries
        p.toCrossViews()              -> list of CrossViews for all positions

      For read-only access, normal dictionaries are faster than CrossViews::
      
        for d in p.iterDicts():       -> iterate over normal dictionaries
        p.toDicts()         -> list of normal (but disconnected) dictionaries 

      Adding a profile to a ProfileCollection will also 'magically'
      add an additional key to all existing CrossViews of
      it. CrossViews become invalid if their parent collection is
      garbage collected or re-ordered.

      **Note**: The creation of many CrossViews hampers performance.
      We provide a :class:`ProfileCollection.iterCrossView`
      iterator to loop over these pseudo dictionaries for convenience --
      direct iteration over the elements of a profile (array or list)
      is about 100 times faster. If you need to repeatedly read from
      many dictionaries, consider using :class:`ProfileCollection.toDicts`
      and cache the resulting normal (disconnected) dictionaries::
      
        cache = p.toDicts() # 'list( p.iterDicts() )' is equivalent but slower
        

    Note
    ====

      Profile arrays of float or int are automatically converted to
      arrays of type Float32 or Int32. This is a safety measure
      because we have stumbled over problems when transferring
      pickled objects between 32 and 64 bit machines. With the transition
      to numpy, this may not be needed any longer.

      See also: :class:`ProfileCollection.__picklesave_array`
    """

    __CrossView = CrossView

    def __init__( self, profiles=None, infos=None ):

        self.profiles = profiles or {}
        self.infos = infos or {}

        self.initVersion = biskit.__version__

        #: re-create this field to invalidate CrossViews! (see :class:`killViews()`)
        self._viewSignal = _ViewSignal()


    def __setstate__(self, state ):
        """
        called for unpickling the object.
        Compability fix: Convert Numeric arrays to numpy arrays 
                         -- requires old Numeric
        """
        self.__dict__ = state

        try:
            import Numeric
        except:
            return

        for k, v in self.profiles.items():

            if getattr( v, 'astype', 0) and not isinstance( v, N.ndarray):
                self.profiles[k] = N.array( v )
    

    def __getitem__( self, k ):
        """
        Get profile item::
          p['prof1']         <==>  p.get( 'prof1' )         
          p['prof1','info1'] <==>  p.get( 'prof1','info1' )
          p[10]              <==>  CrossView( p, 10 )

        :return: profile OR meta infos thereof OR CrossView dict
        :rtype: list OR array OR any OR CrossView
        """
        if isinstance(k, (int, N.integer)):
            return CrossView( self, k )

        if type(k) is slice:
            return self.__getslice__( k )

        return self.get( k )
            

    def __setitem__( self, k, v ):
        """
        Set profile item::
          p['prof1'] = range(10)      <==> p.set( 'prof1', range(10) )      
          p['prof1','info1']='comment' <==> p.setInfo('prof1',info1='comment')

        :return: item
        :rtype: any        
        """
        if type(k) == tuple:
            return self.setInfo( k[0], **{k[1]:v} )

        return self.set( k, v )


    def __delitem__( self, k ):
        """
        Delete profile item::
          del p['prof1']         <==>  p.remove( 'prof1' )          
          del p['prof1','info1'] <==>  p.remove( 'prof1', 'info1' ) 
        """
        result = self.remove( k )


    def __getslice__( self, *arg ):
        """
        Get list of CrossViews::
          p[0:100:5] <==> [ CrossView(p,i) for i in range(0,100,5) ]
        """
        s = arg[0]
        if len( arg ) > 1:
            s = slice( *arg )

        indices = range( *s.indices( self.profLength() ) )

        return [ CrossView(self,i) for i in indices ]

    
    def __len__( self ):
        """
        :return: number of profiles in the collection
        :rtype: int        
        """
        return len( self.profiles )


    def __contains__( self, k ):
        """
        Check if profile contains key::
          k in self  <==>  p.has_key( k ) 

        :return: True or False
        :rtype: 1|0  
        """
        return k in self.profiles


    def __iter__(self):
        """
        Iterate over profile::
          for k in self  <==>  for k in p.keys()
        
        :return: list of items
        :rtype: list
        """
        return iter(self.profiles)


    def keys( self ):
        return self.profiles.keys()


    def has_key( self, k ):
        return k in self.profiles


    def values( self ):
        """
        Get list of all profiles (arrays or lists of values)::
          p.values() -> [ [any], [any], ... ]

        :return: list of lists or arrays
        :rtype: [ list/array ]       
        """
        return [ self.get( key ) for key in self.keys() ]

    def hasNoneProfile(self):
        """
        Check wether any profile is None, which means it is 
        waiting to be updated from a source ProfileCollection.
        This method is written such that it is *not* triggering the 
        updating mechanism.
        :return bool
        """
        for v in self.profiles.values():
            if v is None:
                return True
        return False
            

    def items( self ):
        """
        Get list of tuples of profile names and profiles::
          p.items() -> [ (key1, [any]), (key2, [any]), ..) ]

        :return: list of tuples of profile names and profiles
        :rtype: [ ( str, list/array ) ]       
        """
        return [ (key, self.get(key)) for key in self.keys() ]

    def iteritems(self):
        """
        Iterate over (key : profile) pairs:
          >>> for key, profile in p.iteritems():
          ...
        """
        for key in self.keys():
            yield (key, self.get( key ) )


    def iterCrossViews(self):
        """
        Iterate over values of all profiles as :class:`CrossView` 'dictionaries'
        indexed by profile name, for example:
          >>> for atom in p.iterCrossViews():
          ...     print atom['name'], atom['residue_name']

        The CrossViews remain connected to the profiles and can be
        used to change values in many profiles simultaneously.
        Consider using the somewhat faster :class:`ProfileCollection.iterDicts`
        if this is not needed and speed is critical.

        :return: CrossView instances behaving like dictionaries
        :rtype: iterator over [ CrossView ]
        """
        for i in range( self.profLength() ):
            yield self.__CrossView(self, i)


    def iterDicts(self):
        """
        Iterate over (copies of) values of all profiles as normal dictionaries
        indexed by profile name, for example:
          >>> for atom in p.iterCrossViews():
          ...     print atom['name'], atom['residue_name']

        :return: dictionaries
        :rtype: iterator over { 'key1':v1, 'key2':v2 }
        """

        keys = self.keys()
        profs= zip( *[ self.get( k ) for k in keys ] )

        for values in profs:

            yield dict( zip(keys, values) )


    def toCrossViews(self):
        """
        :return: list of CrossView pseudo dictionaries
        :rtype: [ CrossView ]
        """
        return [self.__CrossView(self, i) for i in range(self.profLength()) ]

    def toDicts(self):
        """
        :return: (copies of) values of all profiles as normal dictionaries
        :rtype: [ dict ]
        """
        keys = self.keys()
        profs= zip( *[ self.get( k ) for k in keys ] )

        ## a bit faster than list( self.iterDicts() )
        return [ dict( zip(keys, values) ) for values in profs ]



    def __picklesave_array( self, prof ):
        """
        Convert integer arrays to Int32 and float arrays to
        Float32. This function is needed because of Numeric issues
        when pickles are transferred between 64 and 32 bit
        platforms. Rather than letting numpy.array select a
        platform-dependent default type, it is safer to assign an
        explicit Float32, or Int32.

        :param prof: the profile array to recast
        :type  prof: numpy.array

        :return: recast array or unchanged original
        :rtype:  Numpy.array
        """
        if prof.dtype.char in ['i','l']:
            return prof.astype( 'i' )

        if prof.dtype.char in ['f','d']:
            return prof.astype( 'f' )

        return prof


    def array_or_list( self, prof, asarray ):
        """
        Convert to array or list depending on asarray option
        
        Beware: empty lists will be upgraded to empty Float arrays.

        :param prof: profile
        :type  prof: list OR array
        :param asarray: 1.. autodetect type, 0.. force list, 2.. force array
        :type  asarray: 2|1|0
        
        :return: profile
        :rtype: list OR array
        
        :raise ProfileError:
        """
        try:

            ## autodetect type
            if asarray == 1:

                if isinstance( prof, N.ndarray ):
                    return self.__picklesave_array( prof )

                if type( prof ) is str:  # tolerate strings as profiles
                    return list( prof )
    
                if len(prof) == 0: # don't create arrays from empty lists
                    return list( prof )
                
                p = self.__picklesave_array( N.array( prof ) )
                if p.dtype.char not in ['O','c','S', 'U']: ## no char or object arrays!
                    return p

                return list( prof )

            ## force list
            if asarray == 0:

                if isinstance( prof, N.ndarray ):
                    return prof.tolist()

                return list( prof )

            ## force array
            if asarray == 2:
                if isinstance( prof, N.ndarray ):
                    return self.__picklesave_array( prof )
                
                return self.__picklesave_array( N.array( prof ) )

        except TypeError as why:
            ## Numeric bug: N.array(['','','']) raises TypeError
            if asarray == 1 or asarray == 0:
                return list( prof )

            raise ProfileError("Cannot create array from given list. %r"\
                  % T.lastError())

        raise ProfileError("%r not allowed as value for asarray" % asarray)


    def expand( self, prof, mask, default ):
        """
        Expand profile to have a value also for masked positions.

        :param prof: input profile
        :type  prof: list OR array
        :param mask: atom mask
        :type  mask: [int]
        :param default: default value
        :type  default: any
        
        :return: profile
        :rtype: list OR array
        """
        if mask is not None:

            ## optimized variant for arrays
            if isinstance( prof, N.ndarray ):
                p = N.resize( prof, (len(mask), ) )
                p[:] = default
                N.put( p, N.nonzero( mask )[0], prof )
                return p

            p = [ default ] * len( mask )
            prof.reverse()
            for i in N.nonzero( mask )[0]:
                p[i] = prof.pop()
            return p

        return prof


    def set( self, name, prof, mask=None, default=None, asarray=1,
             comment=None, **moreInfo ):
        """
        Add/override a profile. None is allowed as special purpose value - in
        which case all other parameters are ignored. Otherwise, the two info
        records 'version', 'changed' and 'isarray' are always modified but can
        be overridden by key=value pairs to this function.
        
        :param name: profile name (i.e. key)
        :type  name: str
        :param prof: list of values OR None
        :type  prof: [any] OR None
        :param mask: list 1 x N_items of 0|1, if there are less values than
                      items, provide mask with 0 for missing values,
                      N.sum(mask)==N_items
        :type  mask: [int]
        :param default: value for items masked.
                        (default: None for lists, 0 for arrays]
        :type  default: any
        :param asarray: store as list (0), as array (2) or store numbers as
                        array but everything else as list (1) (default: 1)
        :type  asarray: 0|1|2
        :param comment: goes into info[name]['comment']
        :type  comment: str
        :param moreInfo: additional key-value pairs for info[name]
        :type  moreInfo: key=value
        
        :raise ProfileError: if length of prof != length of other profiles
        :raise ProfileError: if mask is given but N.sum(mask) != len(prof)
        """
        if prof is None:
            self.profiles[ name ] = None
            self.infos[name] = {}

            if not 'changed' in moreInfo: moreInfo['changed'] = 0
            self.setInfo( name, asarray=asarray, comment=comment,
                         **moreInfo )           
            return

        ## consistency check
        if mask is not None and N.sum(mask) != len(prof):
            raise ProfileError(
                "Mask doesn't match profile ( N.sum(mask)!=len(prof) ). " +
                "%i != %i" % (N.sum(mask), len( prof ) ) )

        prof = self.array_or_list( prof, asarray )

        ## use default == 0 for arrays
        if not default and isinstance( prof, N.ndarray ):
            default = 0

        ## expand profile to have a value also for masked positions
        prof = self.expand( prof, mask, default )

        l = self.profLength()
        if l and len( prof ) != l:
            raise ProfileError( "Profile %s has wrong length." % name )

        ## collect additional infos about this profile
        info = self.infos.get( name, {} )

        info['version'] = '%s %s' % (T.dateString(), biskit.__version__ )
        if comment: info['comment'] = comment
        info['isarray'] = isinstance( prof, N.ndarray )
        info['default'] = default

        ## optional infos
        info.update( moreInfo )

        ## new profiles are always changed=1, updated profiles are checked
        if not 'changed' in moreInfo:
            if name in list(self.keys()):
                info['changed'] = self.infos[name]['changed'] or \
                                  not M.arrayEqual( self.profiles[name], prof )
            else:
                info['changed'] = 1

        ## put profile into dict
        self.profiles[ name ] = prof
        self.infos[ name ] = info


    def setInfo( self, name, **args ):
        """
        Add/Override infos about a given profile::
          e.g. setInfo('relASA', comment='new', params={'bin':'whatif'})
        
        :raise  ProfileError: if no profile is found with |name|
        """
        self.getInfo( name ).update( args )


    def setMany( self, profileDict, infos={} ):
        """
        setMany( dict, [infoDict] ) Add/Override many profiles

        :param profileDict: dict with name:profile pairs
        :type  profileDict: dict
        :param infos: info dicts for each profile, indexed by name
        :type  infos: dict of dict
        """
        for key, value in profileDict.items():
            self.set( key, value, **infos.get( key,{} ) )


    def get( self, name, default=None ):
        """
        get( profKey, [default] ) -> list of values 
        **OR** 
        get( (profKey, infoKey), [default] ) -> single value of info dict
        
        :param name: profile key or profile and info key
        :type  name: str OR (str, str)
        :param default: default result if no profile is found,
                        if None and no profile is found, raise exception
        :type  default: any
        
        :raise ProfileError: if no profile is found with |name|
        """
        ## get an info value
        if type( name ) == tuple:
            result = self.getInfo( name[0] ).get( name[1], default )

            if result is None and name[1] not in self.getInfo( name[0] ):
                raise ProfileError( 'No info value found for '+str(name[1]) )

            return result

        ## get a profile
        result = self.profiles.get( name, default )

        ## but tolerate profiles that are set to None -> return None
        if result is None and name not in self.profiles:
            raise ProfileError( 'No profile found with name '+str(name) )

        return result


    def getInfo( self, name ):
        """
        Use::
           getInfo( name ) -> dict with meta infos about profile::
           
        Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1|0

        :param name: profile name
        :type  name: str

        :return: dict with infos about profile
        :rtype: dict
        
        :raise ProfileError: if no profile is found with |name|
        """
        result = self.infos.get( name, None )

        if result is None:
            raise ProfileError( 'No profile info found with name '+str(name))

        return result


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        Convert profile into a mask based on the max and min cutoff values.
        
        :param profName: profile name
        :type  profName: str
        :param cutoff_min: lower limit
        :type  cutoff_min: float
        :param cutoff_max: upper limit
        :type  cutoff_max: float
        
        :return: mask len( get(profName) ) x 1|0
        :rtype: [1|0]
        """
        p = self.get( profName )

        cutoff_min = cutoff_min or min( p ) - 1
        cutoff_max = cutoff_max or max( p ) + 1

        return N.greater_equal( p, cutoff_min ) * N.less( p, cutoff_max )


    def take( self, indices, *initArgs, **initKw ):
        """
        Take from profiles using provided indices::
          take( indices ) -> ProfileCollection with extract of all profiles

        Any additional parameters are passed to the constructor of the
        new instance.

        :param indices: list of indices
        :type  indices: [int]

        :return: new profile from indices
        :rtype: ProfileCollection (or sub-class)
        
        :raise ProfileError: if take error
        """
        result = self.__class__( *initArgs, **initKw )

        try:
            for key in self.profiles:

                prof = self.get( key )

                if isinstance( prof, N.ndarray ):
                    result.set( key, N.take( prof, indices ) )
                else:
                    result.set( key, [ prof[i] for i in indices ], asarray=0 )

                result.setInfo( key, **copy.deepcopy(self.getInfo(key)) )
                result.setInfo( key, changed=1 )

        except Exception as why:
            raise ProfileError( "Can't take sub-profile %r: %r" % (key,why) )

        return result


    def compress( self, cond ):
        """
        Extract using a mask::
          p.compress( mask ) <==> p.take( N.nonzero( mask ) )

        :param cond: mask with 1 for the positions to keep
        :type  cond: array or list of int 
        """
        return self.take( N.flatnonzero( cond ) )


    def remove( self, *key ):
        """
        Remove profile **OR** info values of profile::
          remove( profKey ) -> 1|0, 1 if complete entry has been removed
          remove( profKey, infoKey ) -> 1|0, 1 if single info value was removed

        :param key: profile name OR name, infoKey
        :type  key: str OR str, str
        
        :return: sucess status
        :rtype: 1|0
        """
        try:
            if len( key ) == 2:
                del self.infos[ key[0] ][ key[1] ]

            else:
                del self.profiles[ key[0] ]
                del self.infos[ key[0] ]

        except KeyError:
            return 0

        return 1

    def __clonedefault(self, a, length, default):
        """
        array mirroring shape and type of <a> but with length <length>
        and filled with default value
        """
        x = a if isinstance(a, N.ndarray) else N.array(a)
        s = list(x.shape)
        s[0] = length
        
        r = N.empty(s, dtype=x.dtype)
        r.fill(default)
        
        if isinstance(a, list):
            return r.tolist()
        return r

    def concat( self, *profiles ):
        """
        Concatenate all profiles in this with corresponding profiles in the
        given ProfileCollection(s). Profiles that are not found in all
        ProfileCollections are skipped::
          p0.concat( p1 [, p2, ..]) -> single ProfileCollection with the
          same number of profiles as p0 but with the length of p0+p1+p2..

        :param profiles: profile(s) to concatenate
        :type  profiles: ProfileCollection(s)
        
        :return: concatenated profile(s)  
        :rtype: ProfileCollection / subclass
        """
        ## end recursion (no more arguments)
        if len( profiles ) == 0:
            return self

        next = profiles[0]

        r = self.__class__()
        
        ##!!! BIG FAT WARNING: empty profilecollection does not imply empty model
        ## an empty PC w/o any profiles currently doesn't know which length
        ## is is supposed to have. If profLength == 0 for real, then
        ## the next PC's profiles don't need to be skipped
        ## Otherwise,
        ## this creates too-short profiles if the PC parent model has 
        ## non-zero length and simply doesn't have any profiles registered.
        
##        ## special case 1: concat something to empty profile collection
##        if not self.keys():
##            return next.clone().concat( *profiles[1:] )
##
##        ## special case 2: concat empty profile collection to this one
##        if not next.keys():
##            return self.clone().concat( *profiles[1:] )
##                
        allkeys = M.union( list(self.profiles.keys()), list(next.keys()) )

##        for k, p in self.profiles.items():
        for k in allkeys:
            p = self.profiles.get(k, None)
            pnext = next.profiles.get(k, None)
            infos = {}
            
            if p is None:
                default = next[k,'default']
                p = self.__clonedefault(pnext, self.profLength(), default)
                infos = next.infos[k]
                
            if pnext is None:
                default = self[k,'default']
                pnext = self.__clonedefault(p, next.profLength(), default)
                infos = self.infos[k]

            try:
                if isinstance( p, N.ndarray ):
                    
                    if len(pnext) == 0:
                        pnext = pnext.astype(p.dtype)
                        
                    r.set( k, N.concatenate( (p, pnext) ), **infos )
                    
                else:
                    r.set( k, p + pnext, **infos )
            except:
                EHandler.warning("Profile %s skipped during concat." % k, 
                                 error=1)
                r.remove( k )

        return r.concat( *profiles[1:] )


    def update( self, other, stickyChanged=1, mask=None ):
        """
        Merge other ProfileCollection into this one, replacing existing
        profiles and info values. This is the obvious translation of
        dict.update(). The changed flag of each profile is set to 1 if:
           1. an existing profile is overridden with different values
           2. the profile is marked 'changed' in the other collection

        The two ProfileCollections should have the same dimension in terms of
        atoms, that is p1.profLength() == p2.profLength(). If this is not the
        case, it is possible to 'mask' atoms in p1 that are 
        are missing in p2. That means the target ProfileCollection can
        have more atoms then the other collection but not vice-versa.
        
        Example::
            p1.profLength() == 10
            p2.profLength() == 5
            p1.update( p2, mask=[0,1,0,1,0,1,0,1,0,1] )
        ...would assign the atom values of the shorter collection p2 to every
        second atom of the longer collection p1.
        If p2 has more items (atoms) per profile than p1, this would not work.
        In this case p2 first needs to be compressed to the same shape
        as p1::
            p1.profLength() == 5
            p2.profLength() == 10
            p2 = p2.compress( [0,1,0,1,0,1,0,1,0,1] )
            p1.update( p2 )
          
        :param other: profile
        :type  other: ProfileCollection
        :param stickyChanged: mark all profiles 'changed' that are marked
                              'changed' in the other collection (default: 1)
        :type  stickyChanged: 0|1
        :param mask: 1 x N_atoms array of 0|1, if the other collection has less
                     atoms than this one, mark those positions with 0 that are
                     only existing in this collection 
                     (N.sum(mask)==self.profLength())
        :type  mask: [int]
        :param mask: N.array of 0 or 1, 
        """
        for key, prof in other.items():

            info = copy.copy( other.getInfo( key ) )
            changed = info.get('changed',0)

            if stickyChanged:
                if not changed:
                    del info['changed']
            else:
                del info['changed']

            self.set( key, prof, mask=mask, **info )


    def updateMissing( self, source, copyMissing=1, allowEmpty=0,
                       setChanged=0 ):
        """
        Merge other ProfileCollection into this one but do not override
        existing profiles and info records. There is one exception:
        Empty profiles (None or []) are replaced but their info records stay
        untouched. If copyMissing=0, profiles that are existing in source but
        not in this collection, are NOT copied (i.e. only empty profiles are
        replaced).

        For each profile copied from the source the 'changed' flag is reset
        to |setChanged| (default 0), regardless whether or not the profile is
        marked 'changed' in the source collection.
        
        :param source: profile
        :type  source: ProfileCollection
        :param copyMissing: copy missing profiles that exist in source
                            (default: 1)
        :type  copyMissing: 0|1
        :param allowEmpty: still tolerate zero-length profiles after update
                           (default: 0)
        :type  allowEmpty: 0|1
        :param setChanged: label profiles copied from source as 'changed' [0]
        :type  setChanged: 0|1
        
        :raise ProfileError: if allowEmpty is 0 and some empty profiles
                             cannot be found in source
        """
        for key, prof in source.items():

            ## replace "None" profiles
            if key in self and self.profiles[ key ] in (None, []):
                self.set( key, prof, changed=setChanged )

            ## add profiles that exist in source but not yet in this collection
            if copyMissing and not key in self:
                info = copy.copy( source.getInfo( key ) )
                info['changed'] = setChanged

                self.set( key, prof, **info )

        if not allowEmpty and self.hasNoneProfile():
            for key, prof in self.profiles.items():
                if not prof:
                    raise ProfileError(('Trying to update %s profile but cannot find'\
                           + ' it in source.') % key)


    def isChanged( self, keys=None ):
        """
        :param keys: only check these profiles (default: None -> means all)
        :type  keys: [ str ] OR str
        :return: True, if any of the profiles is tagged as 'changed'
        :rtype: bool
        """
        keys = keys or list(self.keys())
        keys = T.toList( keys )  ## in case single str is given as argument 

        for k in keys:
            if self.getInfo( k )['changed']:
                return True

        return False
        

    def clone( self ):
        """
        Clone (deepcopy) profiles::
          clone() -> ProfileCollection (or sub-class)

        :return: profile
        :rtype: ProfileCollection          
        """
        r = self.__class__(profiles=copy.deepcopy(self.profiles),
                           infos=copy.deepcopy(self.infos))
        return r


    def killViews(self):
        """
        Deactivate any CrossView instances referring to this ProfileCollection.
        """
        self._viewSignal = _ViewSignal()

    def clear( self ):
        """
        Delete all::
          clear() -> None; delete all profiles and infos.
        """
        self.killViews()

        self.profiles = {}
        self.infos = {}
        

    def profLength( self, default=0 ):
        """
        Length of profile::
          profLength() -> int; length of first non-None profile or default (0)
        
        :param default: value to return if all profiles are set to None
        :type  default: any
          
        :return: length of first non-None profile or 0
        :rtype: int     
        """
        for k, p in self.profiles.items():

            if p is not None:
                return len( p )

        return default


    def plot( self, *name, **arg ):
        """
        Plot one or more profiles using Biggles::
          plot( name1, [name2, ..],[arg1=x, arg2=y]) -> biggles.FramedPlot
        
        :param name: one or more profile names
        :type  name: str
        :param arg: key=value pairs for Biggles.Curve() function
        :type  arg: 
        :raise TypeError: if profile contains non-number items

        :return: plot, view using plot.show() 
        :rtype: biggles.FramedPlot

        :raise ImportError: If biggles module could not be imported
        """
        if not biggles:
            raise ImportError('module biggles could not be imported')

        plot = biggles.FramedPlot()

        colors = T.colorSpectrum( len( name ) , '00FF00', 'FF00FF') 

        for i in range( len(name)):

            p = N.array( self.get( name[i] ) )

            if p.dtype in ['O','c']:
                raise TypeError('Cannot plot values of profile %s.' % name[i])

            # Biggles with its old Numeric cannot handle numpy arrays
            plot.add( biggles.Curve( list(range( len(p))), list(p), color=colors[i],
                                     **arg ) )

            plot.add( biggles.PlotLabel( 0.8, 0.8-i/8.0, name[i],
                                         color=colors[i]) )

        return plot

    
    def plotArray( self, *name, **arg ):
        """
        Plot several profiles as a panel of separate plots.
        :param *name: one or more profile names or tuples of profile names
        :type  *name: str or (str, str,...)
        :param xkey: profile to be used as x-axis (default: None)
        :type  xkey: str
        :param arg: key=value pairs for Biggles.Curve() function
        :type  arg: 

        :return: plot, view using plot.show() 
        :rtype: biggles.FramedPlot

        :raise TypeError: if profile contains non-number items
        :raise ImportError: If biggles module could not be imported
        """
        if not biggles:
            raise ImportError('module biggles could not be imported')

        plot = biggles.FramedArray( len(name),1 )

        if not 'size' in arg:
            arg['size'] = 1
        
        xkey = arg.get('xkey', None)
        if xkey:
            plot.xlabel = xkey
        
        for i, keys in enumerate(name):
            
            if not type(keys) is tuple:
                keys = ( keys, )
                
            for j, key in enumerate(keys):

                x = self.get( xkey, list(range(self.profLength())) )
                y = self[key]
                
                colors = [0] + T.colorSpectrum( len(keys) , '00FF00', 'FF00FF')

                plot[i,0].add( biggles.Curve( x, y, color=colors[j], **arg ) )
    
                plot[i,0].add( biggles.PlotLabel( 0.7, 0.95-j*0.1, key,
                                                  color=colors[j] ) )

        return plot
    
    def plotHistogram( self, *name, **arg ):
        """
        :param bins: number of bins (10)
        :type  bins: int
        :param ynormalize: normalize histograms to area 1.0 (False)
        :type  ynormalize: bool
        :param xnormalize: adapt bin range to min and max of all profiles (True)
        :type  xnormalize: bool
        :param xrange: min and max of bin range (None)
        :type  xrange: (float, float)
        :param steps: draw histogram steps (True)
        :type  steps: bool
        """
        if not biggles:
            raise ImportError('module biggles could not be imported')

        plot = biggles.FramedArray( len(name),1 )

        if not 'size' in arg:
            arg['size'] = 1
        
        bins = arg.get('bins', 10)
        hist = arg.get('ynormalize', False)
        steps= arg.get('steps', 1)
        histrange= arg.get('xrange', None)
        autorange= arg.get('xnormalize', histrange is None)
        
        xkey = arg.get('xkey', None)
        if xkey:
            plot.xlabel = xkey

        if autorange:
            v = []
            for keys in name:
                if not type(keys) is tuple:
                    keys = ( keys, )
                for key in keys:
                    v += list(self[key])
            histrange = ( min( v ), max( v ) )
            
        for i, keys in enumerate(name):
            
            if not type(keys) is tuple:
                keys = ( keys, )
                
            for j, key in enumerate(keys):

                h = density( self[key], nBins=bins, steps=steps, hist=hist,
                             range=histrange )
                x = h[:,0]
                y = h[:,1]
                
                colors = [0] + T.colorSpectrum( len(keys) , '00FF00', 'FF00FF')

                plot[i,0].add( biggles.Curve( x, y, color=colors[j], **arg ) )
    
                plot[i,0].add( biggles.PlotLabel( 0.7, 0.95-j*0.1, key,
                                                  color=colors[j] ) )

        return plot

        
    def __shortString( self, s, maxLen ):
        if len( s ) <= maxLen:
            return s

        return s[:maxLen-4] + '...'


    def __repr__( self ):
        """
        :return: string representation within interactive python interpreter.
        :rtype: str
        """
        s = "ProfileCollection: %i profiles of length %i\n" % \
            (len( self ), self.profLength() )
        for k in self.keys():
            s += k + '\n'
            s += str(self.getInfo(k)) + '\n'
            s += '\t' + self.__shortString( str(self.profiles[k]), 50 ) + '\n'
        return s


#############
##  TESTING        
#############
import biskit.test as  BT
        
class Test(BT.BiskitTest):
    """Test"""

    def test_ProfileCollection( self ):
        """ProfileCollection test"""
        import string

        self.p = ProfileCollection()

        self.p.set( 't1', list(range(10)), comment='test 1', option='x' )
        self.p.set( 't2', list(range(12,22)), comment='second test', option='y' )

        mask = N.zeros( 10 )
        mask[0:10:2] = 1
        l = [ s for s in string.ascii_letters[:5] ] ## list of letters

        self.p.set( 't3', l, comment='masked test', option='z',
                    mask=mask, default=99, asarray=0 )

        if self.local:
            print('\nmasked profile: ', repr( self.p['t3'] ))

        self.p = self.p.take( list(range(0,10,2)) )

        if self.local:
            print('unmasked profile: ', repr( self.p['t3'] ))

        self.p2 = ProfileCollection()
        self.p2.set( 't1', self.p['t1'], comment='overridden', changed=1 )
        self.p2.set( 't4', list(range(30, 35)), comment='added' )

        self.r = self.p.concat( self.p, self.p )  ## concatenate 3 copies of p

        self.p.update( self.p2, stickyChanged=1 )

        self.assertTrue( N.all( self.r['t1'] ==\
                      [0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8]))
    

    def test_concat(self):
        """ProfileCollection.concat test"""
        import string
        import random

        self.p3 = ProfileCollection()
        self.p3.set( 'letters', string.ascii_letters )
        self.p3.set( 'numbers', list(range(len(string.ascii_letters))) )
        self.p3.set( 'random', [ random.randint(0,10000)
                                 for i in range(len(string.ascii_letters)) ] )

        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        self.p3 = self.p3.concat( self.p3, self.p3, self.p3 )
        
        empty = ProfileCollection()
        double = self.p3.concat( self.p3 )
        self.p5 = empty.concat( self.p3 )
        self.assertTrue( N.all(self.p3['letters'] == self.p5['letters']) )
        
        self.p5 = self.p3.concat( empty, empty, self.p3 )
        self.assertTrue( N.all(double['letters'] == self.p5['letters']) )
    
    def test_concatempty(self):
        p0 = ProfileCollection()
        p1 = ProfileCollection()
        p1.set('test', [0,0], asarray=1,)


    def test_crossvies(self):
        """ProfileCollection.crossviews test"""
        import string
        
        self.p4 = ProfileCollection()
        self.p4['letters'] = string.ascii_letters
        
        self.assertTrue( self.p4[0]['letters'] == 'a' )
        views = self.p4.toCrossViews()
        
        letters = ''.join( [ v['letters'] for v in views ] )
        self.assertEqual( letters, string.ascii_letters )
        
        del self.p4
        
        self.assertTrue( views[0].alive is False ) 
        
    def test_plots(self):
        """ProfileCollection.plot* test (only in interactive mode)"""
        import random
        self.p5 = ProfileCollection()
        self.p5['range'] = list(range( 5000))
        self.p5['gamma'] = [ N.random.gamma(3.0) for i in range(5000) ]
        self.p5['normal'] = [ N.random.normal(3.0) for i in range(5000) ]
        self.p5['poisson'] = [ N.random.poisson(5.) for i in range(5000) ]

        if biggles:
            plot = self.p5.plotHistogram( 'gamma', ('normal','poisson'), 
                                        bins=20, hist=1, xnormalize=True )
            if self.local:            
                plot.show()


def clock( s ):

    import profile
    profile.run( s, 'report.out' )

    ## Analyzing
    import pstats
    p = pstats.Stats('report.out')

    ## long steps and methods calling them
    p.sort_stats('cumulative').print_stats(20)
    p.print_callers(0.0)

if __name__ == '__main__':

    BT.localTest()

