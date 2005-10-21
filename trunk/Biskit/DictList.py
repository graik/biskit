##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$
"""organise, sort, and filter list of dictionaries or similar objects"""

import sys
import Numeric as N
import types
import biggles 
import copy
import random

import Biskit.tools as t
from Biskit import EHandler, BisList, ConditionError, AmbiguousMatch,\
     ItemNotFound, BisListError


class DictList( BisList, list ):
    """
    List of dictionaries (or similar objects). Special support is
    given to access and use entries in the dictionaries for filtering
    and sorting. Some care is taken to avoid the adding of
    non-allowed objects (but somehow it is, certainly, still
    possible).

    The list can be sorted and filtered by any dictionary key,
    values within the dictionaries are easily accessible and can be
    plotted against each other.

    Overriding / Extending:
    
    The class is designed to be easily adapted for holding more complex
    objects. In order to do so, the attribute item_type should be
    set to the allowed class. All objects with a
    get( key, default ) and a keys()
    method should then work out of the box. If the two methods are missing,
    getItemValue and getItemKeys must be overriden. 
    """

    def __init__(self, lst=[], item_type=dict ):
        """
        lst        - [ dict ]
        item_type  - type, class of allowed items [ dict ]
        !! raise BisListError, if list contains non-item_type item.
        """
        BisList.__init__( self )
        
        self.item_type = item_type

        if lst != []:
            self.extend( lst )
      

    def version( self ):
        return 'DictList $Revision$'


    def getItemValue( self, item, key, default=None ):
        """
        Get a value from a given item (dictionary). Override this
        method to use the DictList for entries without get() function.
        item    - any, possible entry of this list
        key     - any, dictionary key
        default - any, return value if key is not found [None]
        -> any
        """
        return item.get( key, default )


    def getItemKeys( self, item ):
        """
        Get the attribute keys used by a certain item. Override this
        method to use the DictList for entries without keys() function.
        item    - any, possible entry of this list
        -> [ any ], list of keys
        """
        return item.keys()


    def getValue( self, i, key, default=None ):
        """
        Get the value of a dictionary entry of a list item. 
        i       - int, position in list
        key     - any, dictionary key
        default - any, return value if key is not found [None]
        -> any
        """
        return self.getItemValue( self[i], key, default )
        

    def checkType( self, v ):
        """Make sure v is a dictionary"""
        if not isinstance(v, self.item_type):
            raise BisListError(
                str( v ) + " not allowed. DictList requires "+\
                str(self.item_type))


    def _processNewItem( self, v, i=None ):
        """
        Called before an item is added to the list. Override but call.
        v - dict, value
        i - int, anticipated index (ignored in this implementation)
        """
        self.checkType( v )
        return v

    def _processNewItems( self, lst, indices=None ):
        """
        Called before list of items is added to the list.
        For efficiency, the single new items are only processed if lst is
        not already an instance of this class (DictList by default).
        Override but call.
        lst     - [ any ], list of new items
        indices - [ int ], anticipated positions of the new items [None]
                  defaults to indices following the current end of list
        """
        if not isinstance( lst, list ):
            raise BisListError("Wrong argument type: "+str(type(lst)))

        if not isinstance( lst, self.__class__ ):

            indices = indices or range( len(self), len(self)+len(lst) )

            lst = [self._processNewItem(v, i) for v,i in zip( lst, indices )]
            
        return lst


    def __setitem__(self, i, v ):
        """ lst[i] = v """
        v = self._processNewItem( v, i )
        list.__setitem__( self, i, v)

    def __getslice__( self, i, j ):
        """-> new instance with only the given range of items"""
        r = self.__class__(super(DictList, self).__getslice__(i,j))
        return r

    def extend( self, lst ):
        """extend( list ). Add all items to (the end of) this instance"""
        lst = self._processNewItems( lst )
        list.extend( self, lst )

    def append( self, v ):
        """append( dict ). Append item to the end of this list."""
        v = self._processItem( v, len( self ) )
        list.append( self, v )

    def take( self, indices ):
        """
        indices - array/list of int, list positions
        -> DictList (or sub-class) with specified items
        """
        r = self.__class__( [ self[i] for i in indices ] )

        return r

    def keys( self ):
        """
        -> [ any ], attribute keys used by the current items.
        """
        r = []

        for x in self:
            for k in self.getItemKeys( x ):
                if not k in r:
                    r.append( k )

        return r

    def argsortRandom( self ):
        """ argsortRandom() -> [ int ], indices in random order."""
        r = range( len( self ) )
        random.shuffle( r )
        return r


##
## TESTING
##

if __name__ == '__main__':

    print "DOING something"

    l = DictList()

    for i in range( 10 ):
        d = {'random':random.random(), 'name':'A'}
        l += [ d ]

        
    p = l.plotArray( 'index', 'random', 'random' )

    p.show()
