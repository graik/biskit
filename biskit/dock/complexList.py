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
## organise, sort, and filter Complexes
##

"""
List of Complex objects.
"""

import biskit.core.oldnumeric as N0
import types
import random

import biskit
import biskit.tools as t
from biskit import PDBError, EHandler
from biskit.errors import BiskitError

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.dock; __package__ = "biskit.dock"

from .complex import Complex
from .complexModelRegistry import ComplexModelRegistry

try:
    import biggles
except:
    biggles = 0


class ComplexListError( BiskitError ):
    pass

class ConditionError( ComplexListError ):
    pass

class AmbiguousMatch( ConditionError ):
    pass

class ComplexNotFound( ConditionError):
    pass


class ComplexList( list ):
    """
    List of Complex objects. Special support is given to access and use
    entries in the Complexes' info dictionaries for filtering and sorting.
    Some care is taken to avoid the adding of non-Complex objects (but
    somehow it is, certainly, still possible).

    The different rec- and lig_models of all Complexes are centrally kept in a
    ComplexModelRegistry. Before adding a Complex, we check
    whether its rec- or lig_model are equivalent (same fileName and unchanged)
    to one in the registy. If so, they are replaced to
    avoid that each Complex points to its own copy of essentially the same
    PDBmodel. This does only work with PDBModels that have been pickled to
    disc (see PDBModel.saveAs ) and haven't been changed since. It's a very
    good idea to do this if you want to perform any distributed calculation
    on the ComplexList. Saved PDBModels cause only little PVM traffic (not much
    more than the file name is transmitted). By contrast, unsaved ones will
    severly slow down the job distribution.

    @todo: Removing items with pop(), del, remove() etc. will not remove
           unused PDBModels from rec_models or lig_models. 
    """

    def __init__(self, lst=[] ):
        """
        @param lst: list of Complexes
        @type  lst: [Complex]
        """
        ## non-redundant rec/lig_models of all complexes indexed by file name
        self.models = ComplexModelRegistry()

        self.initVersion = t.dateString() + ' ' + biskit.__version__

        if lst != []:
            self.extend( lst )


    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        self.__dict__ = state
        ## backwards compability
        self.__defaults() 


    def __defaults( self ):
        self.models = getattr( self, 'models', ComplexModelRegistry() )
        if getattr( self, 'rec_models', 0) != 0:
            EHandler.warning(
                're-creating model registry..re-pickle this list!')
            for c in self.toList():
                self.models.addComplex( c )
            del self.rec_models
            del self.lig_models


    def checkType( self, v ):
        """
        Make sure v is a Complex.

        @raise ComplexListError: if not a Complex
        """
        if not isinstance(v, Complex):
            raise ComplexListError(
                str( v ) + " not allowed. ComplexList requires "+str(Complex))


    def checkTypes( self, lst ):
        """
        Make sure lst is a list of Complexes.

        @raise  ComplexListError: if list contains non-Complex item.
        """
        if not isinstance( lst, list ):
            raise ComplexListError("Wrong argument type: "+str(type(lst)))

        if not isinstance( lst, ComplexList ):
            for v in lst:
                self.checkType( v )


    def __setitem__(self, i, v ):
        """
        Set value v of position i.
          >>> lst[i] = v

        @param i: list index
        @type  i: int
        @param v: value
        @type  v: any
        """
        self.checkType( v )

        if i < list.__len__( self ):
            self.models.removeComplex( self[i] )
        self.models.addComplex( v )

        list.__setitem__( self, i, v)


    def __add__( self, lst ):
        """
        New ComplexList with the two lists.

        @return: new instance with simply one list appended to the other
        @rtype:  ComplexList
        """
        r = self.__class__( self )
        r.extend( lst )
        return r


    def __iadd__( self, lst ):
        """
        List appendet to this ComplexList.

        @return: this instance with lst appended
        @rtype: ComplexList
        """
        self.extend( lst )
        return self


    def ligModels( self ):
        """
        Get all shared ligand PDBModels. Stray models (changed or unpickled)
        are not returned.

        @return: list of ligand models
        @rtype: [PDBModel]
        """
        return self.models.ligModels()


    def recModels( self ):
        """
        Get all shared receptor PDBModels. Stray models (changed or unpickled)
        are not returned.

        @return:  list of receptor models
        @rtype: [PDBModel]
        """
        return self.models.recModels()


    def __add_once( self, item, lst ):
        """
        Add eithem to list of it is not already there.
        """
        if lst is None:
            lst = []
        if not item in lst:
            lst.append( item )
        return lst


    def strayModels( self ):
        """
        Look for models that are not in model registry.

        @return: { int|str:[ PDBModel ] }, { int|str:[ PDBModel ] }
        @rtype: dict

        @note: mostly for DEBUGGING
        """
        stray_ligs = {}
        stray_recs = {}
        known_recs = self.recModels()
        known_ligs = self.ligModels()
        for c in self:
            if c.rec_model not in known_recs:
                key = c.get( 'model1', None ) or c.rec_model.fileName \
                    or 1
                stray_recs[ key ] = self.__add_once( c.rec_model,
                                                     stray_recs.get( key, []) )

            if c.lig_model not in known_ligs:
                key = c.get( 'model2', None ) or c.lig_model.fileName \
                    or 1
                stray_ligs[ key ] = self.__add_once( c.lig_model,
                                                     stray_ligs.get( key, []) )

        return stray_recs, stray_ligs


    def extend( self, lst ):
        """
        extend( list ). Add all items to (the end of) this instance
        """
        self.checkTypes( lst )

        list.extend( self, lst )

        for v in lst:
            self.models.addComplex( v )


    def append( self, v ):
        """
        append( Complex ). Append Complex to the end of this list.
        """
        self.checkType( v )
        self.models.addComplex( v )
        list.append( self, v )


    def __getslice__( self, i, j ):
        """
        Slice list.

        @param i: start index
        @type  i: int
        @param j: end index
        @type  j: int

        @return: new instance with only the given range of items
        @rtype: ComplexList
        """
        r = self.__class__(super(ComplexList, self).__getslice__(i,j))
        return r


    def argsortRandom( self ):
        """
        Indices for key in random order::
          argsortRandom() -> [ int ], indices in random order.

        @return: indices in random order
        @rtype: [int]
        """
        pairs = [(random.random(), i) for i in range(0, len(self))]
        pairs.sort()
        return [ x[1] for x in pairs ]


    def argsort( self, sortKey ):
        """
        Indices sort order for values of key::
          argsort( str_sortKey ) -> [ int ], indices after sorting

        @param sortKey: key to use for sorting
        @type  sortKey: str       

        @return: indices after sorting
        @rtype: [int]
        """
        pairs = [(self[i].info.get(sortKey), i) for i in range(0, len(self))]
        pairs.sort()
        return [ x[1] for x in pairs ]


    def sortBy( self, sortKey ):
        """
        Sort ComplexList by key::
          sortBy( str_sortKey ) -> ComplexList (or sub-class)
                                   sorted by info[ str_sortKey ]

        @param sortKey: key to use for sorting
        @type  sortKey: str

        @return: sorted ComplexList
        @rtype: ComplexList
        """
        return self.take( self.argsort( sortKey ))


    def valuesOf(self, infoKey, default=None, indices=None, unique=0 ):
        """
        Get all values of a certain info record of all or some Complexes.

        @param infoKey: key for info dict
        @type  infoKey: str
        @param default: default value if infoKey is not found (None)
        @type  default: any
        @param indices: list of int OR None(=all), indices of Complexes (None)
        @type  indices: [int] OR None
        @param unique: report each value only once (set union), (default 0)
        @type  unique: 1|0

        @return: list of values
        @rtype: [any]
        """
        l = self
        if indices is not None:
            l = N0.take( N0.array(l,'O'), indices )

        if not unique:
            return [ c.info.get(infoKey, default) for c in l ]

        r = []
        for c in l:
            if c.info.get(infoKey, default) not in r:
                r += [ c.info.get( infoKey ) ]
        return r


    def filterRange( self, infoKey, vLow, vHigh ):
        """
        Get indices of Complexes where vLow <= c.info[ infoKey ] <= vHigh.

        Use::
           filterRange( str_infoKey, vLow, vHigh )

        @param infoKey: key for info dict
        @type  infoKey: str
        @param vLow: upper value limit
        @type  vLow: float
        @param vHigh: lower value limit
        @type  vHigh: float

        @return: array of int
        @rtype: [int]
        """
        vLst = self.valuesOf( infoKey )

        maskL = N0.greater_equal( vLst, vLow )
        maskH = N0.less_equal( vLst, vHigh )

        return N0.nonzero( maskL * maskH )


    def filterEqual( self, infoKey, lst ):
        """
        Get indices of Complexes where c.info[ infoKey ] in lst.

        Use::
           filterEqual( infoKey, lst )

        @param infoKey: key for info dict
        @type  infoKey: str
        @param lst: list of values to look for
        @type  lst: [any]

        @return: array of int
        @rtype: [int]
        """
        mask = [ c.info.get( infoKey ) in lst for c in self ]
        return N0.nonzero( mask )


    def filterFunct( self, f ):
        """
        Get indices of Complexes where f( c ) == 1.

        Use::
           filterFunct( f )

        @param f: filterFunct
        @type  f: function

        @return: array of int
        @rtype: [int]
        """
        mask = [ f( c ) for c in self ]
        return N0.nonzero( mask )


    def filter( self, infoKey, cond ):
        """
        Complexes matching condition.

        @param infoKey: key of Complex.info dict
                        (not used if cond is function )
        @type  infoKey: str
        @param cond: filter condition::
                - (vLow, vHigh) -> vLow <= c.info[ infoKey ] <= vHigh
                - list          -> c.info[ infoKey ] in cond
                - function      -> cond( c ) == 1
        @type  cond: any

        @return: ComplexList (or sub-class)
        @rtype: ComplexList

        @raise ConditionError: if cond is neither list nor tuple nor function 
        """
        indices = None

        if type( cond ) == tuple:

            indices = self.filterRange( infoKey, cond[0], cond[1] )

        if type( cond ) == list:

            indices = self.filterEqual( infoKey, cond )

        if type( cond ) == types.FunctionType:

            indices = self.filterFunct( cond )

        if indices is None:
            try:
                indices = self.filterEqual( infoKey, [cond] )
            except:
                raise ConditionError( "Can't interprete filter condition.")

        return self.take(indices)


    def argmax( self, infoKey ):
        """
        Get index of complex c with highest c.infos[infokey] value

        @param infoKey: key for info dict
        @type  infoKey: str

        @return: index of complex c with highest c.infos[infokey] value
        @rtype: int
        """
        vLst = self.valuesOf( infoKey )
        return N0.argmax( vLst )


    def max( self, infoKey ):
        """
        Get higest c.infos[infokey] value.

        @param infoKey: key for info dict
        @type  infoKey: str

        @return: with highest c.info['infoKey'] value
        @rtype: Complex
        """
        return self[ self.argmax(infoKey) ]


    def argmin( self, infoKey ):
        """
        Get index of complex c with lowest c.infos[infokey] value

        @param infoKey: key for info dict
        @type  infoKey: str

        @return: index of complex c with lowest c.infos[infokey] value
        @rtype: int
        """
        vLst = self.valuesOf( infoKey )
        return N0.argmin( vLst )

    def min( self, infoKey ):
        """
        Get lowest c.infos[infokey] value.

        @param infoKey: key for info dict
        @type  infoKey: str

        @return: with lowest c.info['infoKey'] value
        @rtype: Complex
        """
        return self[ self.argmin( infoKey ) ]


    def getIndex( self, infoKey, value ):
        """
        Get list position of Complex where c.info['infoKey'] == value

        @param value: vaue to look for
        @type  value: any
        @param infoKey: key for info dict
        @type  infoKey: str

        @return: position in ComplexList where c.info['infoKey'] == value
        @rtype: int

        @raise AmbiguousMatch: ComplexNotFound,
                               if there are more or less than 1 matches
        """
        l = self.filterEqual( infoKey, [ value ] )

        if len( l ) == 1:
            return l[0]

        if len( l ) > 1:
            raise AmbiguousMatch('More than one Complexes match.')

        raise ComplexNotFound("No matching Complex.")


    def getItem( self, infoKey, value ):
        """
        Get Complex from ComplexList where c.info['infoKey'] == value

        @param value: vaue to look for
        @type  value: any
        @param infoKey: key for info dict
        @type  infoKey: str        

        @return: Complex where c.info['infoKey'] == value
        @rtype: Complex

        @raise AmbiguousMatch: ComplexNotFound,
                               if there are more or less than 1 matches
        """
        return self[ self.getIndex( infoKey, value ) ]


    def take( self, indices ):
        """
        Take the complexes specified by indices.

        @param indices: array/list of int, list positions
        @type  indices: [int]

        @return: ComplexList with all items specified.
        @rtype: ComplexList
        """
        r = self.__class__( [ self[i] for i in indices ] )

        return r


    def toDict( self, infoKey ):
        """
        Convert list into dict indexed by a certain info-record value.
        If several Complexes have the same value, the result will have
        a list registered for this key.

        EXAMPLE:
          >>> clst.toDict('soln') -> {1:Complex, 3:Complex, solnN:Complex}

        @param infoKey: key of info dict in Complexes
        @type  infoKey: str

        @return: { info1:Complex, info2:Complex, info3:[Complex, Complex].. }
        @rtype: dict
        """
        result = {}
        for c in self:
            t.dictAdd( result, c.info[infoKey], c )

        return result


    def toList( self ):
        """
        Convert ComplexList to simple python list of Complexes

        @return: simple python list of Complexes
        @rtype: [Complex]
        """
        return list( self )


    def __maskNone( self, l1, l2 ):
        """
        Take out positions from l1 and l2 that are None in either of them.

        @param l1: first ComplexList
        @type  l1: ComplexList
        @param l2: second ComplexList
        @type  l2: ComplexList

        @return: (l1, l2) modified lists
        @rtype: ComplexList, ComplexList
        """
        r1, r2 = [],[]

        for i in range( len(l1)):
            if l1[i] is not None and l2[i] is not None:
                r1 += [ l1[i] ]
                r2 += [ l2[i] ]

        return r1, r2


    def plot( self, xkey, *ykey, **arg ):
        """
        Plot pairs of info values. The additional arg arguments are handed
        over to biggles.Points().::
          plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) -> biggles.FramedPlot

        @param xkey: key specifying x-values
        @type  xkey: str
        @param ykey: key specifying y-values
        @type  ykey: str OR [str]
        @param arg: additional biggles arguments
        @type  arg: key=value

        @return: biggles plot object
        @rtype: biggles.FramedPlot()
        """
        if not biggles:
            raise ImportError('biggles module could not be imported.')
        
        plot = biggles.FramedPlot()

        plot.xlabel = xkey

        colors = t.colorSpectrum( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = self.valuesOf( xkey )
            y = self.valuesOf( ykey[i] )

            x, y = self.__maskNone( x, y )

            plot.add( biggles.Points( x, y, color=colors[i], **arg ) )

            plot.add( biggles.PlotLabel( 0.2, 0.95-i/8.0, ykey[i],
                                         color=colors[i] ) )

        return plot


    def plotArray( self, xkey, *ykey, **arg ):
        """
        Plot pairs of info values.::
          plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) -> biggles.FramedArray

        @param xkey: key specifying x-values
        @type  xkey: str
        @param ykey: key specifying y-values
        @type  ykey: str OR [str]
        @param arg: additional biggles arguments
        @type  arg: key=value

        @return: biggles plot object
        @rtype: biggles.FramedArray        
        """
        if not biggles:
            raise ImportError('biggles module could not be imported.')
        
        plot = biggles.FramedArray( len(ykey),1 )

        plot.xlabel = xkey

        colors = t.colorSpectrum( len( ykey ) )

        if not 'size' in arg:
            arg['size'] = 1

        for i in range( len(ykey)):

            x = self.valuesOf( xkey )
            y = self.valuesOf( ykey[i] )

            x, y = self.__maskNone( x, y )

            plot[i,0].add( biggles.Points( x, y, color=colors[i], **arg ) )

            plot[i,0].add( biggles.PlotLabel( 0.2, 0.95, ykey[i],
                                              color=colors[i] ) )

        return plot


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_ComplexList(self):
        """Dock.ComplexList test"""
        self.cl = t.load( t.testRoot() + "/dock/hex/complexes.cl" )

        ## number of clusters among the 100 best (lowest rmsd) solutions
        self.cl_sorted = self.cl.sortBy( 'rms' )
        self.hex_clst = self.cl_sorted.valuesOf( 'hex_clst',
                                                 indices=list(range(100)),
                                                 unique=1 )

        if self.local:
            self.p = self.cl.plot( 'rms', 'hex_eshape', 'hex_etotal' )
            self.p.show()

        self.assertEqual( len( self.hex_clst ), 36)

if __name__ == '__main__':

    BT.localTest()

