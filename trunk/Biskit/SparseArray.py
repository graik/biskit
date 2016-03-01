## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## last $Author$
## last $Date$
## $Revision$

"""
Memory saving representation of a sparse array.
"""

import numpy as N
import types
import copy

class SparseArrayError( Exception ):
    pass


## to be transferred into Biskit.tools
def isType( o, t ):
    """
    Test for correct type or correct class::
      isType( o, type_or_class ) -> 1|0

    @param o: object to test
    @type  o: any
    @param t: type OR class
    @type  t: any
    @return: result of test
    @rtype: 1|0
    """
    tt = type(o)
    if tt == types.TypeType:
        return type( o ) == t
    if tt == types.ClassType:
        return isinstance( o, t )
    raise Exception, 'unsupported argument type: %s.' % str(tt)


## to be transferred into Biskit.tools
def getType( o ):
    """
    Get the type or (if o is instance) class of an object::
      getType( o ) -> type or Class

    @param o: get the type OR class of this object
    @type  o: any    
    """
    r = type( o )
    if r == types.InstanceType:
        return o.__class__
    return r


class SparseArray:
    """
    Memory - saving representation of a sparse array.

    Most positions of a sparse array should have a certain default
    value (e.g. 0) so that it is more efficient to store the
    non-default values and their positions rather than storing
    all the repeating zeros.

    A typical application example are atom contact matrices that are
    large (N_atoms x N_atoms) but consist mostly of zeros.
    """

    def __init__( self, array_or_shape=0, typecode='f', default=0. ):
        """
        Create a sparse array from a normal Numeric array or create
        an empty sparse array with a certain shape.

        @param array_or_shape: craeate sparse array from::
                                 ( int, ),   shape of array
                                  OR int,     length of array
                                  OR array,   numeric (dense) array
        @type  array_or_shape: ( int, ) OR int OR array
        @param typecode: single char type of values ['f' OR input type]
        @type  typecode: str
        @param default: value of majority of array content (default: 0.)
        @type  default: any
        """
        self.values  = []
        self.indices = []
        self.__default = default
        self.__typecode= typecode

        self.__last_pos = 0  ## cache last position manipulated
        self.__last_i = 0    ## cache last index manipulated

        a = array_or_shape
        atype = type( a )

        if atype is tuple: self.shape = a
        else:
            if atype is int: self.shape = ( a, )
            else:
                if atype is N.ndarray or atype is list:
                    self.shape = N.shape( a )
                else:
                    raise SparseArrayError, '%s argument not allowed.' %\
                          str(atype)

        self.is1D = len( self.shape ) == 1

        if not self.is1D:
            ## multidimensional
            self.__default = SparseArray( self.shape[1:], typecode, default )
            self.__typecode = 'SA'

        if atype is N.ndarray or atype is list :
            self.__setAll( a )


    def __setAll_1D( self, a ):
        """
        Replace content of this sparseArray with values from Numeric array
        or list of numbers -- only for 1-dimensional arrays.

        @param a: array OR list
        @type  a: array OR [ number ]
        """
        if type( a ) is list:
            a = N.array( a, self.__typecode )

        if self.shape != a.shape:
            raise SparseArrayError, 'dimensions not aligned'

        self.indices = N.nonzero( N.logical_not( N.equal(a, self.__default) ) )[0]
        self.indices = self.indices.tolist()

        self.values = N.take( a, self.indices, 0 )
        self.values = self.values.tolist()


    def __setAll( self, a ):
        """
        Replace content of this array with values from (multidimensional)
        Numeric array of numbers or list of lists (of lists..) of numbers.

        @param a: array OR  list of lists
        @type  a: array OR [ [ number ] ]
        """
        if self.is1D:
            return self.__setAll_1D( a )

        if type(a) in [ N.ndarray, list ]:
            if len(a) != self.shape[0]:
                raise SparseArrayError, 'dimensions not aligned'

            ## typecode and default value for lowest dimension
            typecode = self.typecode()
            default  = self.default()

            for i in range( len(a) ):
                x = SparseArray( a[i], typecode, default )
                if not x.__eq__( self.__default ):
                    self.indices.append( i )
                    self.values.append( x )

    def typecode( self ):
        """
        Get type code of object.

        @return: typecode of lowest dimension
        @rtype: str
        """
        if self.is1D:
            return self.__typecode
        return self.__default.typecode()


    def default( self ):
        """
        Get default value, defined in L{__init__}.

        @return: default value for array elements (of lowest dimension)
        @rtype: number        """
        if self.is1D:
            return self.__default
        ## multidimensional
        return self.__default.default()


    def nondefault( self ):
        """
        Get a 1D list of indices that have a non-default value in a raveled
        version of this array. If L.default()==0 this would be equivalent to
        nonzero( ravel( L.toarray() ) ) (except that the Numeric array is
        never constructed).

        @return: list of indices with none default values
        @rtype: [ int ]
        """
        if self.is1D:
            return self.indices

        ## multidimensional
        r = []
        len_axis_B = self.shape[1]
        for (i,a) in zip( self.indices, self.values ):
            r += (N.array( a.nondefault() ) + len_axis_B * i ).tolist()

        return r


    def nondefault_values( self ):
        """
        Get a 1D-list of all values != L.default() in the order that they
        would have in a raveled array. If L.default()==0 this would be
        equivalent to take( L.toarray(), nonzero( ravel( L.toarray() ) ) )
        (except that the Numeric array is never constructed).

        @return: list of none default values
        @rtype: [ number ]
        """
        if self.is1D:
            return self.values

        ## multidimensional
        r = []
        for a in self.values:
            r.extend( a.nondefault_values() )
        return r


    def toarray( self ):
        """
        Reconstruct dense array::
          L.toarray() -> Numeric.array, normal dense array

        @return: normal dense array
        @rtype: array
        """
        if self.default() is 0:
            a = N.zeros( ( self.shape ), self.typecode() )
        else:
            a = N.ones( (self.shape ), self.typecode() ) * self.default()

        N.put( a, self.nondefault(), self.nondefault_values() )

        return a


    def tolist( self ):
        """
        Reconstruct list from dense array::
          L.tolist() -> [ ([) number (]) ], same as L.toarray().tolist()

        @return: list from normal dense array
        @rtype: list        
        """
        return self.toarray().tolist()


    def put( self, i, v ):
        """
        Replace one or several values, L.put( i, v )

        @param i: indices
        @type  i: int OR [ int ]
        @param v: values
        @type  v: any OR [ any ]
        """
        if type( i ) in [ list, N.ndarray ]:
            self.__setMany( i, v )
        else:
            self.__setitem__( i, v )


    def __setMany( self, indices, values ):
        """
        Add / replace values of the array.

        @param indices: indices, [ int ] OR Numeric.array('i')
        @type  indices: [int]
        @param values: values, [ any ] OR Numeric.array
        @type  values: [any] OR array
        """
        if type( values ) not in [ list, N.ndarray ]:
            values = [ values ] * len( indices )

        map( self.__setitem__, indices, values )


    def __pos( self, i ):
        """
        Target position of given index in index list.

        @param i: index
        @type  i: int

        @return: target position of given index in index list
        @rtype: int

        @raise IndexError: if out of bounds
        """
        if i > self.shape[0] or i < 0:
            raise IndexError, "index %i out of bounds" % i

        if i >= self.__last_i:
            pos = self.__last_pos
        else:
            pos = 0

        l = len( self.indices )
        while pos != l and self.indices[ pos ] < i:
            pos += 1

        if pos != l:
            self.__last_i   = self.indices[ pos ]
            self.__last_pos = pos
        else:
            self.__last_i = 0
            self.__last_pos = 0

        return pos


    def __setitem__( self, i, v ):
        """
        Set position specifyed by the index i to value v.

        @param i: index 
        @type  i: int OR (int,int)
        @param v: value
        @type  v: any

        @raise SparseArrayError: if dimensions not aligned
        @raise SparseArrayError: if no sequence value
        """
        itype = type( i )
        vtype = getType( v )

        if itype is tuple and len( i ) == 1:
            i = i[0]
            itype = int

        if itype is int:
            try:
                if not self.is1D:
                    if len( v ) != self.shape[1]:
                        raise SparseArrayError, 'dimensions not aligned.'
                    if vtype is not SparseArray:
                        v = SparseArray( v, self.typecode(), self.default() )

                return self.__setitem_1D( i, v )

            except TypeError:
                raise SparseArrayError, 'sequence value required.'

        if itype is tuple:
            x, y = i[0], i[1:]

            a = self.__getitem__( x )
            a.__setitem__( y, v )
            self.__setitem__( x, a )


    def __setitem_1D( self, i, v ):
        """
        Set position specifyed by the index i to value v in a 1D array.

        @param i: array index
        @type  i: int
        @param v: value
        @type  v: any

        @raise IndexError: if i < 0 or i > len( this )
        """
        if i in self.indices:
            pos = self.indices.index( i )

            if v != self.__default:
                self.values[ pos ] = v
            else:
                del self.values[ pos ]
                del self.indices[ pos ]

        else:
            if v != self.__default:

                pos = self.__pos( i )

                self.indices.insert( pos, i )
                self.values.insert( pos, v )


    def __getitem_1D( self, i ):
        """
        Value for specified position::
          this[ i ] -> number OR SparseArray

        @param i: array index
        @type  i: int         

        @raise IndexError: if i < 0 or i >= len( this )
        """
        if i >= self.shape[0] or i < 0:
            raise IndexError, "index %i out of bounds" % i

        if i in self.indices:
            pos = self.indices.index( i )
            return self.values[ pos ]

        if self.is1D:
            return self.__default

        return copy.deepcopy( self.__default )


    def __getitem__( self, i ):
        """
        Value for specified position.

        @param i: array index
        @type  i: int OR (int,int)
        """
        itype = type( i )

        if itype is tuple and len( i ) == 1:
            i = i[0]
            itype = int

        if itype is int:
            return self.__getitem_1D( i )

        if itype is tuple:
            x, y = i[0], i[1:]
            return self.__getitem__( x ).__getitem__( y )


    def __len__( self ):
        """
        Get length: supports len( this )

        @return: object length
        @rtype: int
        """
        return self.shape[0]


    def __eq__( self, o ):
        """
        Comparison equal: supports this == other

        @return: result of comparison
        @rtype: 0|1
        """
        if not isinstance( o, self.__class__ ):
            return 0
        if self.shape != o.shape:
            return 0
        return self.values == o.values and self.indices == o.indices


    def __ne__( self, o ):
        """
        Comparison not equal: supports this != other

        @return: result of comparison
        @rtype: 0|1        
        """
        if not isinstance( o, self.__class__ ):
            return 1
        if self.shape != o.shape:
            return 1
        return self.values != o.values or self.indices != o.indices


    def __getslice__( self, a, b ):
        """
        Slice a sparce array::
          this[ a : b ] -> SparseArray

        @return: sliced sparse array
        @rtype: SparseArray
        """
        shape = ( abs( b - a ), ) + self.shape[1:]
        result = self.__class__( shape, self.__typecode, self.__default )

        pos_low  = self.__pos( a )
        pos_high = self.__pos( b )

        result.put( N.array( self.indices[pos_low : pos_high] ) - a,
                    self.values[pos_low : pos_high] )
        return result


    def __contains__( self, v ):
        """
        Sparse array contains value: supports v in this -> 0|1

        @param v: value
        @type  v: any      

        @return: result of comparison
        @rtype: 0|1
        """
        return ( cmp( v, self.__default ) == 0 ) or (v in self.values)


    def count(self, v ):
        """
        Count the occuravces of value in sparse array::
          count( value ) -> int, number of occurences of value.

        @param v: value
        @type  v: any      

        @return: number of occurances
        @rtype: int
        """
        if v == self.__default:
            return len( self ) - len( self.values )
        return self.values.count( v )


    def index( self, v ):
        """
        position of first occurrence of value::
          index( value ) -> int

        @param v: value to look for
        @type  v: any

        @return: index of first occurance
        @rtype: int

        @raise ValueError: if value is not contained in this list
        """
        if v in self.values:
            return self.indices[ self.values.index( v ) ]

        if v != self.__default:
            raise ValueError, "SparseArray.index(): value not in list"

        i = 0
        while i in self.indices:
            i += 1

        if i == len( self ):
            raise ValueError, "SparseArray.index(): value not in list"

        return i


    def empty( self ):
        """
        Check if object is empty.

        @return: true if lenght is 0
        @rtype: 1|0
        """
        return len( self.indices ) is 0


    def __delitem__( self, i ):
        """
        Delete value at index i: supports del this[i]

        @note: del this[int:int] is not supported

        @param i: index
        @type  i: int     
        """
        if i >= self.shape[0] or i < 0:
            raise IndexError, "index %i out of bounds" % i

        if i in self.indices:
            pos = self.indices.index( i )

            del self.indices[ pos ]
            del self.values[ pos ]

        else:
            pos = self.__pos( i )

        if pos < len( self.indices ) - 1:
            self.indices = self.indices[:pos] +\
                [ p-1 for p in self.indices[pos:] ]

        self.shape = (self.shape[0] - 1,) + self.shape[1:]


##     def reverse( self ):
##         raise Exception, "not implemented."


    def insert( self, i, v ):
        """
        Insert value before index::
          this.insert(index, value) 

        @param v: value to insert
        @type  v: any
        @param i: index
        @type  i: int

        @raise IndexError: if i < 0 or i > len( this )
        """
        pos = self.__pos( i )

        inc = 0

        if v != self.__default:
            self.indices.insert( pos, i )
            self.values.insert( pos, v )
            inc = 1

        if pos < len( self.indices ):
            self.indices = self.indices[:pos+inc] +\
                [ p+1 for p in self.indices[pos+inc:] ]

        self.shape = (self.shape[0]+1, ) + self.shape[1:]


##     def pop( self ):
##         raise Exception, "not implemented."


##     def remove( self, i ):
##         raise Exception, "not implemented."


##     def __add__( self, lst ):
##         raise Exception, "not implemented."


##     def __iadd__( self, lst ):
##         raise Exception, "not implemented."


    def append( self, v, axis=0 ):
        """
        Enlarge this array along an axis by one value::
          L.append( value, [ axis=0 ])

        @param v: v MUST be a sparse array if the given
                  axis is not the last axis of this array
        @type  v: number OR SparseArray
        @param axis: dimension along wich v should be added
        @type  axis: int

        @raise SparseArrayError: if dimension errors
        """
        if axis == 0:

            t_default, t_v = getType( self.__default ), getType( v )
            if t_default == SparseArray and t_v == N.ndarray:
                v = SparseArray( v, self.typecode(), self.default()  )

            if getType(v) != t_default:
                raise SparseArrayError, 'Cannot append %s to array of %s.' \
                      % ( str(type(v)), str(t_default) )

            if v != self.__default:
                self.indices.append( self.shape[0] )
                self.values.append( v )

            self.shape = ( self.shape[0]+1, ) + self.shape[1:]

        else:
            if len( v ) != self.shape[0]:
                raise SparseArrayError, 'dimensions not aligned'

            for i in range( self.shape[0] ):
                v_i = v[i]
                a_i = self.__getitem__( i )
                is0 = a_i.empty()   ## -> a_i is a copy of self.__default

                if v_i != self.__default.__default:
                    a_i.append( v_i, axis=axis-1 )
                    if is0:         ## otherwise change would be lost
                        self.__setitem__( i, a_i )

            d = self.__default
            x = axis - 1
            d.shape = d.shape[:x] + (d.shape[x]+1,) + d.shape[x+1:]

            self.shape = self.shape[:axis] + ( self.shape[axis]+1, ) +\
                self.shape[axis+1:]


    def extend( self, lst ):
        """
        Extend list by appending elements from the iterable::
          L.extend(iterable)

        @param lst: list OR SparseArray with extend values
        @type  lst: [any] OR SparseArray
        """
        if not isinstance( lst, SparseArray ):
            lst = SparseArray( lst, default=self.__default )

        l = self.shape[0]
        self.indices.extend( [ i + l for i in lst.indices ] )
        self.values.extend( lst.values )
        self.shape = ( l + lst.shape[0], ) + self.shape[1:]


##     def __setitem__( self, i, v ):

##         if getattr( v, 'shape', None ) is None:
##             return SparseList.__setitem__( i, v )

##         pass



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_SparseArray(self):
        """SparseArray test"""

        a = N.zeros( (6,), N.float32 )

        self.sa = SparseArray( a.shape )
        self.sa[3] = 1.
        self.sa[5] = 2.

        b = N.zeros( (5, 6), N.float32 )
        b[0,1] = 3.
        b[0,2] = 4
        b[4,2] = 5
        b[3,0] = 6

        self.sb = SparseArray( b )

        self.sb.append( self.sa )

        if self.local:
            print self.sa.toarray()

        self.assert_( N.all( self.sb.toarray() == self.EXPECTED) )


    EXPECTED = N.array([[ 0.,  3.,  4.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  0.,  0.,  0.],
                        [ 6.,  0.,  0.,  0.,  0.,  0.],
                        [ 0.,  0.,  5.,  0.,  0.,  0.],
                        [ 0.,  0.,  0.,  1.,  0.,  2.]])


if __name__ == '__main__':

    BT.localTest()






