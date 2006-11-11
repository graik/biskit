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
## last $Date$
## last $Author$
## $Revision$

"""
Manage profiles.
"""

import Numeric as N
import tools as T
import mathUtils as M
from Biskit import EHandler

import copy

try:
    import biggles
except:
    biggles = 0
    EHandler.warning('Could not import module biggles.')


class ProfileError(Exception):
    pass

class ProfileCollection:
    """
    Manage profiles (arrays or lists of values) for trajectory frames
    or atoms/residues in PDBModel. ProfileCollection resembles a
    2-dimensional array where the first axis (let's say row) is
    accessed by a string key and each row has an additional info
    dictionary assigned to it. The take() and concat() methods operate
    on the columns, i.e. they are applied to all profiles at the same
    time.

    By default, profiles of numbers are stored and returned as
    Numeric.array and all others are stored and returned as ordinary
    list. This behaviour can be modified with the option asarray of
    ProfileCollection.set(). Using both lists and arrays is a
    compromise between the efficiency of Numeric arrays and the
    problem that arrays of objects cannot be unpickled (Numeric bug)
    and that arrays of strings would end up as 2-D arrays of char.
    The 'isarray' entry of a profile's info dictionary tells whether
    the profile is stored as array or as list.

    ProfileCollection p can be accessed like a dictionary of lists::
      len( p )          -> number of profiles (== len( p.profiles ) )
      p['prof1']        -> list with values of profile 'prof1'
      del p['prof1']    -> remove a profile
      p['prof1'] = [..] -> add a profile without additional infos
      for k in p        -> iterate over profile keys
      'prof1' in p      -> 1, if collection contains key 'prof1'

    But it is more than that - each key also has a dictionary of info values
    assigned to it (see getInfo(), setInfo(), p.infos). These can be accessed
    like::
      p['prof1','date']   -> date of creation of profile named 'prof1'
      p.getInfo('prof1')  -> returns all info records
      p['prof1','comment'] = 'first prof'  -> add/change single info value
    """

    def __init__( self, version=None, profiles=None, infos=None ):

        self.profiles = profiles or {}
        self.infos = infos or {}

        self.initVersion = version or self.version()


    def version( self ):
        """
        Class version.
        
        @return: class version number
        @rtype: str
        """
        return 'ProfileCollection $Revision$'


    def __getitem__( self, k ):
        """
        Get profile item::
          p['prof1']         <==>  p.get( 'prof1' )         
          p['prof1','info1]  <==>  p.get( 'prof1','info1' ) 

        @return: item
        @rtype: any
        """
        return self.get( k )


    def __setitem__( self, k, v ):
        """
        Set profile item::
          p['prof1'] = range(10)      <==> p.set( 'prof1', range(10) )      
          p['prof1','info1]='comment' <==> p.setInfo('prof1',info1='comment')

        @return: item
        @rtype: any        
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


    def __len__( self ):
        """
        Length of profile

        @return: profile length
        @rtype: int        
        """
        return len( self.profiles )


    def __contains__( self, k ):
        """
        Check if profile contains key::
          k in self  <==>  p.has_key( k ) 

        @return: True or False
        @rtype: 1|0  
        """
        return self.has_key( k )


    def __iter__(self):
        """
        Iterate over profile::
          for k in self  <==>  for k in p.keys()
        
        @return: list of items
        @rtype: list
        """
        return iter(self.profiles)


    def keys( self ):
        return self.profiles.keys()


    def has_key( self, k ):
        return self.profiles.has_key(k)


    def values( self ):
        return self.profiles.values()


    def items( self ):
        """
        Get list of tuples of profile names and profiles::
          p.items() -> [ (key1, [any]), (key2, [any]), ..) ]

        @return: list of tuples of profile names and profiles
        @rtype: list       
        """
        return self.profiles.items()


    def __array_or_list( self, prof, asarray ):
        """
        Convert to array or list depending on asarray option

        @param prof: profile
        @type  prof: list OR array
        @param asarray: option, convert to array
        @type  asarray: 1|0
        
        @return: profile
        @rtype: list OR array
        
        @raise ProfileError:
        """
        ## autodetect type
        if asarray == 1:
            if isinstance( prof, N.arraytype ):
                return prof

            p = N.array( prof )
            if p.typecode() not in ['O','c']:  ## no char or object arrays!
                return p
            return prof

        ## force list
        if asarray == 0:
            if isinstance( prof, N.arraytype ):
                return prof.tolist()
            return prof

        ## force array
        if asarray == 2:
            if isinstance( prof, N.arraytype ):
                return prof
            return N.array( prof )

        raise ProfileError, "%r not allowed as value for asarray" % asarray


    def __expand( self, prof, mask, default ):
        """
        Expand profile to have a value also for masked positions.

        @param prof: profile
        @type  prof: list OR array
        @param mask: atom mask
        @type  mask: [int]
        @param default: default value
        @type  default: any
        
        @return: profile
        @rtype: list OR array
        """
        if mask:

            ## optimized variant for arrays
            if isinstance( prof, N.arraytype ):
                p = N.resize( prof, (len(mask), ) )
                p[:] = default
                N.put( p, N.nonzero( mask ), prof )
                return p

            p = [ default ] * len( mask )
            prof.reverse()
            for i in N.nonzero( mask ):
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
        
        @param name: profile name (i.e. key)
        @type  name: str
        @param prof: list of values OR None
        @type  prof: [any] OR None
        @param mask: list 1 x N_items of 0|1, if there are less values than
                      items, provide mask with 0 for missing values,
                      N.sum(mask)==N_items
        @type  mask: [int]
        @param default: value for items masked.
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
        if prof is None:
           self.profiles[ name ] = None
           return

        ## consistency check
        if mask and N.sum(mask) != len(prof):
            raise ProfileError(
                "Mask doesn't match profile ( N.sum(mask)!=len(prof) ). " +
                "%i != %i" % (N.sum(mask), len( prof ) ) )

        prof = self.__array_or_list( prof, asarray )

        ## use default == 0 for arrays
        if not default and isinstance( prof, N.arraytype ):
            default = 0

        ## expand profile to have a value also for masked positions
        prof = self.__expand( prof, mask, default )

        l = self.profLength()
        if l and len( prof ) != l:
            raise ProfileError( "Profile %s has wrong length." % name )

        ## collect additional infos about this profile
        info = self.infos.get( name, {} )

        info['version'] = '%s %s' % (T.dateString(), self.version() )
        if comment: info['comment'] = comment
        info['isarray'] = isinstance( prof, N.arraytype )

        ## optional infos
        info.update( moreInfo )

        ## new profiles are always changed=1, updated profiles are checked
        if not 'changed' in moreInfo:
            if name in self.keys():
                info['changed'] = self.infos[name]['changed'] or \
                                  not M.arrayEqual( self[name], prof )
            else:
                info['changed'] = 1

        ## put profile into dict
        self.profiles[ name ] = prof
        self.infos[ name ] = info


    def setInfo( self, name, **args ):
        """
        Add/Override infos about a given profile::
          e.g. setInfo('relASA', comment='new', params={'bin':'whatif'})
        
        @raise  ProfileError: if no profile is found with |name|
        """
        self.getInfo( name ).update( args )


    def setMany( self, profileDict, infos={} ):
        """
        setMany( dict, [infoDict] ) Add/Override many profiles

        @param profileDict: dict with name:profile pairs
        @type  profileDict: dict
        @param infos: info dicts for each profile, indexed by name
        @type  infos: dict of dict
        """
        for key, value in profileDict.items():
            self.set( key, value, **infos.get( key,{} ) )


    def get( self, name, default=None ):
        """
        get( profKey, [default] ) -> list of values 
        B{OR} 
        get( (profKey, infoKey), [default] ) -> single value of info dict
        
        @param name: profile key or profile and info key
        @type  name: str OR (str, str)
        @param default: default result if no profile is found,
                        if None and no profile is found, raise exception
        @type  default: any
        
        @raise ProfileError: if no profile is found with |name|
        """
        ## get an info value
        if type( name ) == tuple:
            result = self.getInfo( name[0] ).get( name[1], default )

            if result is None and not self.infos[ name[0] ].has_key(name[1]):
                raise ProfileError( 'No info value found for '+str(name[1]) )

            return result

        ## get a profile
        result = self.profiles.get( name, default )

        ## but tolerate profiles that are set to None -> return None
        if result is None and not self.profiles.has_key(name):
            raise ProfileError( 'No profile found with name '+str(name) )

        return result


    def getInfo( self, name ):
        """
        getInfo( name ) -> dict with infos about profile::
          Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1|0

        @param name: profile name
        @type  name: str

        @return: dict with infos about profile
        @rtype: dict
        
        @raise ProfileError: if no profile is found with |name|
        """
        result = self.infos.get( name, None )

        if result is None:
            raise ProfileError( 'No profile info found with name '+str(name))

        return result


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        Convert profile into a mask based on the max and min cutoff values.
        
        @param profName: profile name
        @type  profName: str
        @param cutoff_min: lower limit
        @type  cutoff_min: float
        @param cutoff_max: upper limit
        @type  cutoff_max: float
        
        @return: mask len( get(profName) ) x 1|0
        @rtype: [1|0]
        """
        p = self.get( profName )

        cutoff_min = cutoff_min or min( p ) - 1
        cutoff_max = cutoff_max or max( p ) + 1

        return N.greater( p, cutoff_min ) * N.less( p, cutoff_max )


    def take( self, indices ):
        """
        Take on profile using provided indices::
          take( indices ) -> ProfileCollection with extract of all profiles

        @param indices: list of indices
        @type  indices [int]

        @return: new profile from indices
        @rtype: profile
        
        @raise ProfileError: if take error
        """
        result = self.__class__( self.version() )

        try:
            for key, prof in self.profiles.items():

                if isinstance( prof, N.arraytype ):
                    result.set( key, N.take( prof, indices ) )
                else:
                    result.set( key, [ prof[i] for i in indices ], asarray=0 )

                result.setInfo( key, **copy.deepcopy(self.getInfo(key)) )

        except Exception, why:
            raise ProfileError( "Can't take sub-profile: "+str(why) )

        return result


    def remove( self, *key ):
        """
        Remove profile B{OR} info values of profile::
          remove( profKey ) -> 1|0, 1 if complete entry has been removed
          remove( profKey, infoKey ) -> 1|0, 1 if single info value was removed

        @param key: profile name OR name, infoKey
        @type  key: str OR str, str
        
        @return: sucess status
        @rtype: 1|0
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


    def concat( self, *profiles ):
        """
        Concatenate all profiles in this with corresponding profiles in the
        given ProfileCollection(s). Profiles that are not found in all
        ProfileCollections are skipped::
          p0.concat( p1 [, p2, ..]) -> single ProfileCollection with the
          same number of profiles as p0 but with the length of p0+p1+p2..

        @param profiles: profile(s) to concatenate
        @type  profiles: profileCollection(s)
        
        @return: concatenated profile(s)  
        @rtype: profileCollection
        """

        if len( profiles ) == 0:
            return self

        next = profiles[0]

        r = self.__class__()

        for k, p in self.profiles.items():

            try:
                if isinstance( p, N.arraytype ):
                    r.set( k, N.concatenate( (p, next.get(k)) ),
                           **self.infos[k] )
                else:
                    r.set( k, p + next.get(k), **self.infos[k] )
            except:
                EHandler.warning("Can't concat profile "+k)
                r.remove( k )

        return r.concat( *profiles[1:] )


    def update( self, other, stickyChanged=1 ):
        """
        Merge other ProfileCollection into this one, replacing existing
        profiles and info values. This is the obvious translation of
        dict.update(). The changed flag of each profile is set to 1 if:
           1. an existing profile is overridden with different values
           2. the profile is marked 'changed' in the other collection
          
        @param other: profile
        @type  other: ProfileCollection
        @param stickyChanged: mark all profiles 'changed' that are marked
                              'changed' in the other collection (default: 1)
        @type  stickyChanged: 0|1
        """
        for key, prof in other.items():

            info = copy.copy( other.getInfo( key ) )
            changed = info.get('changed',0)

            if stickyChanged:
                if not changed:
                    del info['changed']
            else:
                del info['changed']

            self.set( key, prof, **info )


    def updateMissing( self, source, copyMissing=1, allowEmpty=0 ):
        """
        Merge other ProfileCollection into this one but do not replace / update
        existing profiles and info records. There is one exception:
        Empty profiles (None or []) are replaced but their info records stay
        untouched. If copyMissing=0, profiles that are existing in source but
        not in this collection, are NOT copied (i.e. only empty profiles are
        replaced).
        
        @param source: profile
        @type  source: ProfileCollection
        @param copyMissing: copy missing profiles that exist in source
                            (default: 1)
        @type  copyMissing: 0|1
        @param allowEmpty: tolerate zero-length profiles after update
                           (default: 0)
        @type  allowEmpty: 0|1
        
        @raise ProfileError: if allowEmpty is 0 and some empty profiles
                             cannot be found in source
        """
        for key, prof in source.items():

            ## replace "None" profiles
            if key in self and not self[ key ]:
                self.set( key, prof )

            ## add profiles that exist in source but not yet in this collection
            if copyMissing and not key in self:
                info = copy.copy( source.getInfo( key ) )
                del info['changed']

                self.set( key, prof, **info )

        if not allowEmpty and ( None in self.values() or [] in self.values() ):
            for key, prof in self.items():
                if not prof:
                    raise ProfileError, \
                          ('Trying to update %s profile but cannot find'\
                           + ' it in source.') % key


    def clone( self ):
        """
        Clone (deepcopy) profile::
          clone() -> ProfileCollection (or sub-class, actually a deepcopy)

        @return: profile
        @rtype: profileCollection          
        """
        return copy.deepcopy( self )


    def clear( self ):
        """
        Delete all::
          clear() -> None; delete all profiles and infos.
        """
        self.profiles = {}
        self.infos = {}


    def profLength( self ):
        """
        Length of profile::
          profLength() -> int; length of first non-None profile or 0
        
        @return: length of first non-None profile or 0
        @rtype: int     
        """
        for k, p in self.items():

            if p != None:
                return len( p )

        return 0


    def plot( self, *name, **arg ):
        """
        Plot one or more profiles using Biggles::
          plot( name1, [name2, ..],[arg1=x, arg2=y]) -> biggles.FramedPlot
        
        @param name: one or more profile names
        @type  name: str
        @param arg: key=value pairs for Biggles.Curve() function
        @type  arg: 
        @raise TypeError: if profile contains non-number items

        @return: plot, view using plot.show() 
        @rtype: biggles.FramedPlot

        @raise ImportError: If biggles module could not be imported
        """
        if not biggles:
            raise ImportError, 'module biggles could not be imported'

        plot = biggles.FramedPlot()

        colors = T.colorSpectrum( len( name ) , '00FF00', 'FF00FF') 

        for i in range( len(name)):

            p = N.array( self.get( name[i] ) )

            if p.typecode() in ['O','c']:
                raise TypeError, 'Cannot plot values of profile %s.' % name[i]

            plot.add( biggles.Curve( range( len(p) ), p, color=colors[i],
                                     **arg ) )

            plot.add( biggles.PlotLabel( 0.8, 0.8-i/8.0, name[i],
                                         color=colors[i]) )

        return plot


    def __shortString( self, s, maxLen ):
        """
        """
        if len( s ) <= maxLen:
            return s

        return s[:maxLen-4] + '...'


    def __repr__( self ):
        """
        @return: string representation within interactive python interpreter.
        @rtype: str
        """
        s = "ProfileCollection: %i profiles of length %i\n" % \
            (len( self ), self.profLength() )
        for k in self.keys():
            s += k + '\n'
            s += str(self.infos[k]) + '\n'
            s += '\t' + self.__shortString( str(self.profiles[k]), 50 ) + '\n'
        return s


#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype: int
        """
        import string

        p = ProfileCollection()

        p.set( 't1', range(10), comment='test 1', option='x' )
        p.set( 't2', range(12,22), comment='second test', option='y' )

        mask = N.zeros( 10 )
        mask[0:10:2] = 1
        l = [ s for s in string.letters[:5] ] ## list of letters
        p.set( 't3', l, comment='masked test', option='z',
               mask=mask, default=99, asarray=0 )

        if local: print repr( p['t3'] )

        p = p.take( range(0,10,2) )

        if local: print repr( p['t3'] )

        p2 = ProfileCollection()
        p2.set( 't1', p['t1'], comment='overridden', changed=1 )
        p2.set( 't4', range(30, 35), comment='added' )

        r = p.concat( p, p )  ## concatenate 3 copies of p

        p.update( p2, stickyChanged=1 )

        if local:
            globals().update( locals() )
            
        return r['t1']


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return N.array([0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8])
    
        

if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()


