##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$
"""
organise, sort, and filter collection of dictionaries or similar objects
-> abstract base class (aka interface)
"""

import Numeric as N
import types
import biggles
import copy
import random

import Biskit.tools as t
from Biskit import EHandler
from Errors import BiskitError

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
    collections of dictionary-like objects.  To the outside, it
    behaves like a list. However, no assumptions are yet made as to
    the internal data structure.

    Classes derived from AbstractDictList have to override several
    methods to be functional (a NotImplementedError is raised otherwise):

    getValue, extend, append, take, keys,
    __len__, __setitem__, __getslice__
    
    The latter 3 are not yet defined in BisList (no
    NotImplementedError) but are nevertheless required. They can also
    be provided from the built-in list type via multiple inheritence (see
    DictList, for an example).

    That means there are two ways of implementing BisList:
    1) via multiple inheritence from BisList (first!) and list
    -> only getValue, extend, append, take need to be overriden
    2) inheritence from BisList only
    -> the __xxx__ methods have to be implemented too.

    See DictList for an example of strategy 1).
    """

    def __init__(self): 
        """
        Override but call.
        """
        self.initVersion = t.dateString() + ' ' + self.version()

    def version( self ):
        """
        -> str, CVS version of this class (see also attribute initVersion)
        """
        return 'BisList $Revision$'


    def getValue( self, i, key, default=None ): # abstract
        """
        Get the value of a dictionary entry of a list item.
        Override!
        i       - int, position in collection
        key     - any, attribute key
        default - any, return value if key is not found [None]
        -> any
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
        -> new instance with one collection appended to the other
        """
        r = self.__class__( self )
        r.extend( other )
        return r

    def __iadd__( self, other ):
        """
        c.__iadd__( other ) <==> c += other 
        -> this instance with other appended
        """
        self.extend( other )
        return self

    def extend( self, other ): # abstract
        """
        Add all items of other to (the end of) this instance.
        other - AbstractDictList (left to the implementing class)
        Override!
        """
        raise NotImplementedError

    def append( self, v ): # abstract
        """
        Append a single item to the end of this list.
        Override!
        v - any (left to the implementing class)
        """
        raise NotImplementedError

    def keys( self ):
        """
        -> [ any ], attribute keys used by the current items.
        """
        raise NotImplementedError

    def argsort( self, sortKey, cmpfunc=cmp ):
        """
        Sort by values of a certain item attribute.
        sortKey - any, attribute key
        cmpfunc - function, used for comparing values; cmpfunc(v1,v2) -> -1,0,1
                  [built-in cmp]
        -> [ int ], indices after sorting (the collection itself is not sorted)
        """
        pairs = [(self.getValue(i,sortKey),i) for i in range(0, len(self))]
        pairs.sort( cmpfunc )
        return [ x[1] for x in pairs ]


    def take( self, indices, deepcopy=0 ): # abstract
        """
        Extract certain items in a certain order.
        Override!
        indices  - [ int ], positions
        deepcopy - 0|1, deepcopy items [0]
        -> new instance (or sub-class) with specified items
        """
        raise NotImplementedError

    def compress( self, mask, deepcopy=0 ):
        """
        Extract certain items.
        mask     - [ 1|0 ], mask of positions; len( mask ) == len( self )
        deepcopy - 1|0, deepcopy items [0]
        """
        return self.take( N.nonzero( mask ), deepcopy=deepcopy )

    def sortBy( self, sortKey, cmpfunc=cmp ):
        """
        sortBy( sortKey ) -> new instance sorted by item[ sortKey ]
        """
        return self.take( self.argsort( sortKey, cmpfunc ))

    def valuesOf(self, key, default=None, indices=None, unique=0 ): #
        """
        Get all values assigned to a certain key of all or some
        items. The result is guaranteed to have the same length as the
        collection (or the list of given indices). Missing values are
        replaced by default (None).
        key     - any, key for item attribute
        default - any, default value if key is not found [None]
        indices - list of int OR None(=all), indices [None]
        unique  - 1|0, report each value only once (set union), (default 0)
        -> list of values
        """
        l = self
        if indices != None:
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
        key - any, item attribute
        vLow, vHigh - any, lower and upper bound
        -> array of int
        """
        vLst = self.valuesOf( key )

        maskL = N.greater_equal( vLst, vLow )
        maskH = N.less_equal( vLst, vHigh )

        return N.nonzero( maskL * maskH )


    def filterEqual( self, key, lst ):
        """
        Get indices of items for which item[ key ] in lst.
        key - any, item attribute
        lst - [ any ], list of allowed values
        -> array of int
        """
        mask = [ self.getValue( i,key) in lst for i in range( len(self)) ]
        return N.nonzero( mask )

         
    def filterFunct( self, f ):
        """
        Get indices of items for which f( item ) == 1.
        f - function, f must take a single item as argument and return 1 or 0
        -> array of int
        """
        mask = [ f( c ) for c in self ]
        return N.nonzero( mask )


    def filter( self, key, cond ):
        """
        Extract items matching condition.
        key     - item attribute  (not used if cond is function )
        cond    - (vLow, vHigh) -> vLow <= item[ key ] <= vHigh
                - list          -> item[ key ] in cond
                - function      -> cond( c ) == 1
        -> new instance (or sub-class)
        !! ConditionError if cond is neither list nor tuple nor function
        """
        indices = None
        
        if type( cond ) == tuple:
            
            indices = self.filterRange( key, cond[0], cond[1] )

        if type( cond ) == list:

            indices = self.filterEqual( key, cond )

        if type( cond ) == types.FunctionType:

            indices = self.filterFunct( cond )
            
        if indices == None:
            try:
                indices = self.filterEqual( key, [cond] )
            except:
                raise ConditionError( "Can't interprete filter condition.")

        return self.take(indices)

    def argmax( self, key ):
        """
        key - any, item attribute
        -> int, index of item with highest item[key] value
        """
        vLst = self.valuesOf( key )
        return N.argmax( vLst )

    def max( self, key ):
        """
        key - any, item attribute
        -> any, item with highest item[key] value
        """
        return self[ self.argmax(key) ]


    def argmin( self, key ):
        """
        key - any, item attribute
        -> int, index of item with lowest item[infokey] value
        """
        vLst = self.valuesOf( key )
        return N.argmin( vLst )

    def min( self, key ):
        """
        key - any, item attribute
        -> any, item with lowest item[key] value
        """
        return self[ self.argmin( key ) ]

    def getIndex( self, key, value ):
        """
        -> int, position of item for which item[key] == value
        !! AmbiguousMatch, ItemNotFound,
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
        -> any, item for which item[key] == value
        !! AmbiguousMatch, ItemNotFound,
           if there are more or less than 1 matches
        """
        return self[ self.getIndex( key, value ) ]
    
    def toDict( self, key ):
        """
        Convert collection into dict indexed by the value of a certain item
        attribute. If several items have the same value, the result will have
        a list registered for this key.
        key - any, item attribute
        
        EXAMPLE: lst.toDict('soln') -> {soln1:Item, soln3:Item, solnN:Item}

        -> dict, { info1:dict, info2:dict, info3:[dict, dict].. }
        """
        result = {}
        for i in range( len(self)):
            t.dictAdd( result, self.getValue( i, key), self[i] )

        return result


    def toList( self ): 
        """ -> [ item ], simple python list of items """
        return list( self )


    def __maskNone( self, l1, l2 ):
        """
        Take out positions from l1 and l2 that are None in either of them.
        -> (l1, l2) modified lists
        """
        r1, r2 = [],[]

        for i in range( len(l1)):
            if l1[i] != None and l2[i] != None:
                r1 += [ l1[i] ]
                r2 += [ l2[i] ]

        return r1, r2

    def plot( self, xkey, *ykey, **arg ):
        """
        plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) -> biggles.FramedPlot
        Plot pairs of item values. The additional arg arguments are handed
        over to biggles.Points(). The special xkey value 'index' uses the
        position of each item as x-axis. If only one key is given,
        it is taken as ykey and the x-axis is the index of each item
        (xkey='index').
        -> Biggles.FramedPlot, display with show() !
        """
        if len(ykey) == 0:
            xkey, ykey = 'index', [ xkey ]

        plot = biggles.FramedPlot()

        plot.xlabel = xkey

        colors = t.colorSpectrum( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = range( len( self ) )
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
        plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) -> biggles.FramedPlot
        Plot pairs of item values.
        -> Biggles.FramedArray, display with show()
        """
        if len(ykey) == 0:
            xkey, ykey = 'index', [ xkey ]

        plot = biggles.FramedArray( len(ykey),1 )

        plot.xlabel = xkey

        colors = t.colorSpectrum( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = range( len( self ) )
            if xkey != 'index':
                x = self.valuesOf( xkey )

            y = self.valuesOf( ykey[i] )

            x, y = self.__maskNone( x, y )

            plot[i,0].add( biggles.Points( x, y, color=colors[i], **arg ) )

            plot[i,0].add( biggles.PlotLabel( 0.2, 0.95, ykey[i],
                                         color=colors[i] ) )

        return plot

