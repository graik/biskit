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

## last $Author$
## last $Date$
## $Revision$
"""
organise, sort, and filter list of dictionaries or similar objects
"""

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

    I{Overriding / Extending:}
    The class is designed to be easily adapted for holding more complex
    objects. In order to do so, the attribute item_type should be
    set to the allowed class. All objects with a get( key, default )
    and a keys() method should then work out of the box. If the two
    methods are missing, getItemValue and getItemKeys must be overriden. 
    """

    def __init__(self, lst=[], item_type=dict ):
        """
        @param lst: list of dictionaries
        @type  lst: [ dict ]
        @param item_type: class of allowed items [ dict ]
        @type  item_type: type
        
        @raise BisListError: if list contains non-item_type item.
        """
        BisList.__init__( self )

        self.item_type = item_type

        if lst != []:
            self.extend( lst )


    def version( self ):
        """
        Version of class.
        
        @return: version of class
        @rtype: str
        """
        return 'DictList $Revision$'


    def getItemValue( self, item, key, default=None ):
        """
        Get a value from a given item (dictionary). Override this
        method to use the DictList for entries without get() function.
        
        @param item: possible entry of this list
        @type  item: any
        @param key: dictionary key
        @type  key: any
        @param default: return value if key is not found (default: None)
        @type  default: any
        
        @return: any
        @rtype: any
        """
        return item.get( key, default )


    def getItemKeys( self, item ):
        """
        Get the attribute keys used by a certain item. Override this
        method to use the DictList for entries without keys() function.
        
        @param item: possible entry of this list
        @type  item: any
        
        @return: list of keys
        @rtype: [ any ],
        """
        return item.keys()


    def getValue( self, i, key, default=None ):
        """
        Get the value of a dictionary entry of a list item.
        
        @param i: position in list
        @type  i: int
        @param key: dictionary key
        @type  key: any
        @param default: return value if key is not found (default: None)
        @type  default: any
        
        @return: any
        @rtype: any
        """
        return self.getItemValue( self[i], key, default )


    def checkType( self, v ):
        """
        Make sure v is a dictionary

        @raise BisListError: if not a dictionary
        """
        if not isinstance(v, self.item_type):
            raise BisListError(
                str( v ) + " not allowed. DictList requires "+\
                str(self.item_type))


    def _processNewItem( self, v, i=None ):
        """
        Called before an item is added to the list. Override but call.
        
        @param v: value
        @type  v: dict
        @param i: anticipated index (ignored in this implementation)
        @type  i: int

        @return: value
        @rtype: dict      
        """
        self.checkType( v )
        return v


    def _processNewItems( self, lst, indices=None ):
        """
        Called before list of items is added to the list.
        For efficiency, the single new items are only processed if lst is
        not already an instance of this class (DictList by default).
        Override but call.
        
        @param lst: list of new items
        @type  lst: [ any ]
        @param indices:  anticipated positions of the new items (default: None)
                         defaults to indices following the current end of list
        @type  indices: [ int ]

        @return: list
        @rtype: list
        """
        if not isinstance( lst, list ):
            raise BisListError("Wrong argument type: "+str(type(lst)))

        if not isinstance( lst, self.__class__ ):

            indices = indices or range( len(self), len(self)+len(lst) )

            lst = [self._processNewItem(v, i) for v,i in zip( lst, indices )]

        return lst


    def __setitem__(self, i, v ):
        """
        Set value v of position i.
          >>> lst[i] = v

        @param i: list index
        @type  i: int
        @param v: value
        @type  v: any
        """
        v = self._processNewItem( v, i )
        list.__setitem__( self, i, v)


    def __getslice__( self, i, j ):
        """
        Return new instance with only the given range of items.
        
        @param i: start list index
        @type  i: int
        @param j: end list index
        @type  j: int
        
        @return: new instance with only the given range of items
        @rtype: instance
        """
        r = self.__class__(super(DictList, self).__getslice__(i,j))
        return r


    def extend( self, lst ):
        """
        extend( list ). Add all items to (the end of) this instance

        @param lst: list of new items
        @type  lst: [ any ]
        """
        lst = self._processNewItems( lst )
        list.extend( self, lst )


    def append( self, v ):
        """
        append( dict ). Append item to the end of this list.

        @param v: value
        @type  v: any
        """
        v = self._processItem( v, len( self ) )
        list.append( self, v )


    def take( self, indices ):
        """
        @param indices: list positions
        @type  indices: array/list of int
        
        @return: DictList (or sub-class) with specified items
        @rtype: DictList
        """
        r = self.__class__( [ self[i] for i in indices ] )

        return r


    def keys( self ):
        """
        @return: attribute keys used by the current items.
        @rtype: [ any ],
        """
        r = []

        for x in self:
            for k in self.getItemKeys( x ):
                if not k in r:
                    r.append( k )

        return r


    def argsortRandom( self ):
        """
        Random sort.
         
        @return: indices in random order.
        @rtype: [ int ]
        """
        r = range( len( self ) )
        random.shuffle( r )
        return r


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
        l = DictList()

        for i in range( 10 ):
            d = {'random':random.random(), 'name':'A'}
            l += [ d ]

        p = l.plotArray( 'index', 'random', 'random' )

        if local:
            p.show()
            globals().update( locals() )
            
        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
       

if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()


