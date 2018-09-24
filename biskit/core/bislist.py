## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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
organise, sort, and filter collection of dictionaries or similar objects
S{->} abstract base class (aka interface)

Currently only used as base for DictList.
"""

import biskit
import biskit.core.oldnumeric as N0
import types

import biskit.tools as t
import biskit.colorspectrum as C
from biskit import EHandler
from biskit.errors import BiskitError

try:
    import biggles
except:
    biggles = 0


class BisListError( BiskitError ):
    pass

class ConditionError( BisListError ):
    pass

class AmbiguousMatch( ConditionError ):
    pass

class ItemNotFound( ConditionError):
    pass

class BisList:
    """
    An *abstract* base class that lays out the basic interface for
    collections of dictionary-like objects. To the outside, it behaves
    like a list, albeit one with a second dimension addressed by
    'keys'. You could think of BisList as a kind of two-dimensional
    array or a dictionary of lists. However, no assumptions are yet
    made as to the internal data structure. BisList provides uniform
    methods for sorting, filtering, extraction and combination of
    sub-lists (take, compress), and plotting.

    Classes derived from BisList have to override several
    methods to be functional (a NotImplementedError is raised otherwise):

      * getValue, extend, append, take, keys,
      * __len__, __setitem__, __getslice__

    The latter 3 are not yet defined in BisList (no
    NotImplementedError) but are nevertheless required. They can also
    be provided from the python built-in list type via multiple
    inheritence (see L{Biskit.DictList}, for an example).

    That means there are two ways of implementing BisList.
      1. via multiple inheritence from BisList (first!) and list
         -> only L{getValue}, L{extend}, L{append} and L{take} need to be
         overriden.
      2. inheritence from BisList only
         -> the __xxx__ methods have to be implemented, too.

    See L{DictList} for an example of strategy 1.
    """

    def __init__(self): 
        """
        Override but call.
        """
        self.initVersion = t.dateString() + ' ' + biskit.__version__


    def getValue( self, i, key, default=None ): # abstract
        """
        Get the value of a dictionary entry of a list item.
        B{Override!}

        @param i: position in collection
        @type  i: int
        @param key: attribute key
        @type  key: any
        @param default: return value if key is not found [None]
        @type  default: any

        @return: any
        @rtype: any
        """
        raise NotImplementedError

##     def __setitem__(self, i, v ): # abstract
##         """
##         c.__setitem( i, v ) <==> c[i] = v
##         Override!
##         """
##         raise NotImplemented

##     def __getslice__( self, i, j ): # abstract
##         """
##         c.__getslice__( i, j ) <==> c[ i : j ]
##         Override!
##         i - int, starting position (including)
##         j - int, ending position (excluding)
##         -> new instance with only the given range of items
##         """
##         raise NotImplemented

##     def __len__( self ): # abstract
##         """
##         c.__len__() <==> len( c )
##         Override!
##         """
##         raise NotImplemented


    def __add__( self, other ):
        """
        c.__add__( other ) <==> c + other

        @param other: other instance
        @type  other: instance

        @return: new instance with one collection appended to the other
        @rtype: any
        """
        r = self.__class__( self )
        r.extend( other )
        return r


    def __iadd__( self, other ):
        """
        c.__iadd__( other ) <==> c += other

        @param other: other instance
        @type  other: instance

        @return: this instance with other appended
        @rtype: any
        """
        self.extend( other )
        return self


    def extend( self, other ): # abstract
        """
        Add all items of other to (the end of) this instance.
        B{Override!}

        @param other: other instance
        @type  other: instance        
        """
        raise NotImplementedError


    def append( self, v ): # abstract
        """
        Append a single item to the end of this list.
        B{Override!}

        @param v: any (left to the implementing class)
        @type  v: any
        """
        raise NotImplementedError


    def keys( self ):
        """
        @return: attribute keys used by the current items.
        @rtype: [ any ]
        """
        raise NotImplementedError


    def argsort( self, sortKey, cmpfunc=None ):
        """
        Sort by values of a certain item attribute.

        @param sortKey: attribute key
        @type  sortKey: any
        @param cmpfunc: used for comparing values; cmpfunc(v1,v2)
                        -> -1,0,1 [built-in cmp]
        @type  cmpfunc: function

        @return: indices after sorting (the collection itself is not sorted)
        @rtype: [ int ]
        """
        ## cmp vanished in python 3.x (but still available in past.builtins)
        cmpfunc = cmpfunc or ( lambda x, y: (x > y) - (x < y))

        pairs = [(self.getValue(i,sortKey),i) for i in range(0, len(self))]
        pairs.sort( cmpfunc )
        return [ x[1] for x in pairs ]


    def take( self, indices, deepcopy=0 ): # abstract
        """
        Extract certain items in a certain order.
        B{Override!}

        @param indices: positions
        @type  indices: [ int ]
        @param deepcopy: deepcopy items (default: 0)
        @type  deepcopy: 0|1

        @return: new instance (or sub-class) with specified items
        @rtype: instance
        """
        raise NotImplementedError


    def compress( self, mask, deepcopy=0 ):
        """
        Extract certain items.

        @param mask: mask of positions; len( mask ) == len( self )
        @type  mask: [ 1|0 ]
        @param deepcopy: deepcopy items (default: 0)
        @type  deepcopy: 1|0

        @return: new instance (or sub-class) with specified items
        @rtype: instance
        """
        return self.take( N0.nonzero( mask ), deepcopy=deepcopy )


    def sortBy( self, sortKey, cmpfunc=None ):
        """
        Use::
          sortBy( sortKey ) -> new instance sorted by item[ sortKey ]

        @param sortKey: key for item attribute
        @type  sortKey: any
        @param cmpfunc: comparison function
        @type  cmpfunc: function

        @return: new instance (or sub-class) sorted by item
        @rtype: instance        
        """
        ## cmp vanished in python 3.x (but still available in past.builtins)
        cmpfunc = cmpfunc or (lambda x, y: (x > y) - (x < y))
        
        return self.take( self.argsort( sortKey, cmpfunc ))


    def valuesOf(self, key, default=None, indices=None, unique=0 ): 
        """
        Get all values assigned to a certain key of all or some
        items. If unique==0, the result is guaranteed to have the same
        length as the collection (or the list of given
        indices). Missing values are replaced by default (None).

        @param key: key for item attribute
        @type  key: any
        @param default: default value if key is not found (default: None)
        @type  default: any
        @param indices: indices defining a subset of this list (default: None)
        @type  indices: list of int OR None
        @param unique: report each value only once (set union), (default 0)
        @type  unique: 1|0

        @return: list of values
        @rtype: list
        """
        l = self
        if indices is not None:
            l = self.take( indices )

        if not unique:
            return [ self.getValue( i,key,default) for i in range(len(l)) ]

        r = []
        for i in range( len(l) ):
            if self.getValue( i, key, default) not in r:
                r += [ self.getValue( i, key, default ) ]
        return r


    def filterRange( self, key, vLow, vHigh ):
        """
        Get indices of items where vLow <= item[ key ] <= vHigh.

        @param key: item attribute
        @type  key: any
        @param vLow: lower bound
        @type  vLow: any
        @param vHigh: upper bound
        @type  vHigh: any

        @return: array of int
        @rtype: array
        """
        vLst = self.valuesOf( key )

        maskL = N0.greater_equal( vLst, vLow )
        maskH = N0.less_equal( vLst, vHigh )

        return N0.nonzero( maskL * maskH )


    def filterEqual( self, key, lst ):
        """
        Get indices of items for which item[ key ] in lst.

        @param key: item attribute
        @type  key: any
        @param lst: [ any ], list of allowed values
        @type  lst: list

        @return: array of int
        @rtype: array
        """
        mask = [ self.getValue( i,key) in lst for i in range( len(self)) ]
        return N0.nonzero( mask )


    def filterFunct( self, f ):
        """
        Get indices of items for which f( item ) == 1.

        @param f: f must take a single item as argument and return 1 or 0
        @type  f: function

        @return: array of int
        @rtype: array
        """
        mask = [ f( c ) for c in self ]
        return N0.nonzero( mask )


    def filter( self, key, cond ):
        """
        Extract items matching condition.

        @param key: item attribute  (not used if cond is function )
        @type  key: any
        @param cond: conditon::
                     - (vLow, vHigh) -> vLow <= item[ key ] <= vHigh
                     - list          -> item[ key ] in cond
                     - function      -> cond( c ) == 1
        @type  cond: any

        @return: new instance (or sub-class)
        @rtype: instance

        @raise ConditionError: if cond is neither list nor tuple nor function: 
        """
        indices = None

        if type( cond ) == tuple:

            indices = self.filterRange( key, cond[0], cond[1] )

        if type( cond ) == list:

            indices = self.filterEqual( key, cond )

        if type( cond ) == types.FunctionType:

            indices = self.filterFunct( cond )

        if indices is None:
            try:
                indices = self.filterEqual( key, [cond] )
            except:
                raise ConditionError( "Can't interprete filter condition.")

        return self.take(indices)


    def argmax( self, key ):
        """
        @param key: item attribute
        @type  key: any

        @return: index of item with highest item[key] value
        @rtype: int
        """
        vLst = self.valuesOf( key )
        return N0.argmax( vLst )


    def max( self, key ):
        """
        @param key: item attribute
        @type  key: any

        @return: item with highest item[key] value
        @rtype: float
        """
        return self[ self.argmax(key) ]


    def argmin( self, key ):
        """
        @param key: item attribute
        @type  key: any

        @return: index of item with lowest item[infokey] value
        @rtype: int
        """
        vLst = self.valuesOf( key )
        return N0.argmin( vLst )


    def min( self, key ):
        """
        @param key: item attribute
        @type  key: any

        @return: item with lowest item[key] value
        @rtype: float
        """
        return self[ self.argmin( key ) ]


    def getIndex( self, key, value ):
        """
        @param key: item attribute
        @type  key: any
        @param value: item value
        @type  value: any        

        @return: position of item for which item[key] == value
        @rtype: int

        @raise AmbiguousMatch: ItemNotFound,
               if there are more or less than 1 matches
        """
        l = self.filterEqual( key, [ value ] )

        if len( l ) == 1:
            return l[0]

        if len( l ) > 1:
            raise AmbiguousMatch('More than one Complexes match.')

        raise ItemNotFound("No matching item.")


    def getItem( self, key, value ):
        """
        @param key: item attribute
        @type  key: any
        @param value: item value
        @type  value: any

        @return: item for which item[key] == value
        @rtype: any

        @raise AmbiguousMatch: ItemNotFound,
               if there are more or less than 1 matches
        """
        return self[ self.getIndex( key, value ) ]


    def toDict( self, key ):
        """
        Convert collection into dict indexed by the value of a certain item
        attribute. If several items have the same value, the result will have
        a list registered for this key.

        C{ EXAMPLE: lst.toDict('soln')                 }
        C{ -> {soln1:Item, soln3:Item, solnN:Item}     }

        @param key: item attribute
        @type  key: any

        @return: { info1:dict, info2:dict, info3:[dict, dict].. }
        @rtype: dict
        """
        result = {}
        for i in range( len(self)):
            t.dictAdd( result, self.getValue( i, key), self[i] )

        return result


    def toList( self ): 
        """
        @return: simple python list of items
        @rtype: [ item ]
        """
        return list( self )


    def __maskNone( self, l1, l2 ):
        """
        Take out positions from l1 and l2 that are None in either of them.

        @param l1: first list
        @type  l1: list
        @param l2: second list
        @type  l2: list

        @return: modified lists
        @rtype: (l1, l2)
        """
        r1, r2 = [],[]

        for i in range( len(l1)):
            if l1[i] is not None and l2[i] is not None:
                r1 += [ l1[i] ]
                r2 += [ l2[i] ]

        return r1, r2


    def plot( self, xkey, *ykey, **arg ):
        """
        Plot pairs of item values. The additional arg arguments are handed
        over to biggles.Points(). The special xkey value 'index' uses the
        position of each item as x-axis. If only one key is given,
        it is taken as ykey and the x-axis is the index of each item
        (xkey='index').

        C{ EXAMPLE: plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) }
        C{         -> biggles.FramedPlot }

        @param xkey: key for x-values
        @type  xkey: any
        @param ykey: key for y-values
        @type  ykey: any
        @param arg: arguments handed over to biggles.Points()
        @type  arg: any

        @return: Biggles.FramedPlot, display with show() !
        @rtype:  Biggles.FramedPlot
        """
        if not biggles:
            raise ImportError('biggles module could not be imported.')

        if len(ykey) == 0:
            xkey, ykey = 'index', [ xkey ]

        plot = biggles.FramedPlot()

        plot.xlabel = xkey

        colors = C.colorRange( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = list(range( len( self )))
            if xkey != 'index':
                x = self.valuesOf( xkey )

            y = self.valuesOf( ykey[i] )

            x, y = self.__maskNone( x, y )

            plot.add( biggles.Points( x, y, color=colors[i], **arg ) )

            plot.add( biggles.PlotLabel( 0.2, 0.95-i/8.0, ykey[i],
                                         color=colors[i] ) )

        return plot


    def plotArray( self, xkey, *ykey, **arg ):
        """
        Plot pairs of item values.

        C{ EXAMPLE: plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) }
        C{         -> biggles.FramedPlot                           }

        @param xkey: key for x-values
        @type  xkey: any
        @param ykey: key for y-values
        @type  ykey: any
        @param arg: arguments handed over to biggles.Points()
        @type  arg: any

        @return: Biggles.FramedArray, display with show()
        @rtype: Biggles.FramedArray
        """
        if not biggles:
            raise ImportError('biggles module could not be imported.')

        if len(ykey) == 0:
            xkey, ykey = 'index', [ xkey ]

        plot = biggles.FramedArray( len(ykey),1 )

        plot.xlabel = xkey

        colors = C.colorRange( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = list(range( len( self )))
            if xkey != 'index':
                x = self.valuesOf( xkey )

            y = self.valuesOf( ykey[i] )

            x, y = self.__maskNone( x, y )

            plot[i,0].add( biggles.Points( x, y, color=colors[i], **arg ) )

            plot[i,0].add( biggles.PlotLabel( 0.2, 0.95, ykey[i],
                                              color=colors[i] ) )

        return plot

##############
## empty test
##############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, the BisList is tested in L{Biskit.DictList}"""
    pass
