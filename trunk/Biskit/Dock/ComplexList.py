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
## organise, sort, and filter Complexes

## last $Author$
## last $Date$
## $Revision$

import Numeric as N
import types
import biggles
import copy
import random

import Biskit.tools as t
from Biskit import BiskitError, PDBError, EHandler

from Biskit.Dock.Complex import Complex
from Biskit.Dock.ComplexModelRegistry import ComplexModelRegistry


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

    ToDo:
    Removing items with pop(), del, remove() etc. will not remove unused
    PDBModels from rec_models or lig_models. 
    """

    def __init__(self, lst=[] ):
        """
        lst        - list of Complexes
        rec_models - internal rec_models dict, fname : PDBModel
        lig_models - internal lig_models dict
        If both rec_models and lig_models are given, the model synchronisation
        and type checking is skipped. (for internal use, see take() )
        !! raise ComplexListError, if list contains non-Complex item.
        """
        ## non-redundant rec/lig_models of all complexes indexed by file name
        self.models = ComplexModelRegistry()

        self.initVersion = t.dateString() + ' ' + self.version()

        if lst != []:
            self.extend( lst )
      

    def version( self ):
        return 'ComplexList $Revision$'


    def __setstate__(self, state ):
        """called for unpickling the object."""
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
        """Make sure v is a Complex"""
        if not isinstance(v, Complex):
            raise ComplexListError(
                str( v ) + " not allowed. ComplexList requires "+str(Complex))


    def checkTypes( self, lst ):
        """ !! raise ComplexListError, if list contains non-Complex item. """
        if not isinstance( lst, list ):
            raise ComplexListError("Wrong argument type: "+str(type(lst)))

        if not isinstance( lst, ComplexList ):
            for v in lst:
                self.checkType( v )

    def __setitem__(self, i, v ):
        self.checkType( v )

        if i < list.__len__( self ):
            self.models.removeComplex( self[i] )
        self.models.addComplex( v )

        list.__setitem__( self, i, v)

    def __add__( self, lst ):
        """-> new instance with simply one list appended to the other"""
        r = self.__class__( self )
        r.extend( lst )
        return r

    def __iadd__( self, lst ):
        """-> this instance with lst appended"""
        self.extend( lst )
        return self

    def ligModels( self ):
        """
        Get all shared ligand PDBModels. Stray models (changed or unpickled)
        are not returned.
        -> [ PDBModel ]
        """
        return self.models.ligModels()

    def recModels( self ):
        """
        Get all shared receptor PDBModels. Stray models (changed or unpickled)
        are not returned.
        -> [ PDBModel ]
        """
        return self.models.recModels()

    def __add_once( self, item, lst ):
        if lst == None:
            lst = []
        if not item in lst:
            lst.append( item )
        return lst

    def strayModels( self ):
        """
        mostly for DEBUGGING
        Look for models that are not in model registry.
        -> { int|str:[ PDBModel ] }, { int|str:[ PDBModel ] }
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
        """extend( list ). Add all items to (the end of) this instance"""
        self.checkTypes( lst )

        list.extend( self, lst )

        for v in lst:
            self.models.addComplex( v )

    def append( self, v ):
        """append( Complex ). Append Complex to the end of this list."""
        self.checkType( v )
        self.models.addComplex( v )
        list.append( self, v )

    def __getslice__( self, i, j ):
        """-> new instance with only the given range of items"""
        r = self.__class__(super(ComplexList, self).__getslice__(i,j))
        return r


    def argsortRandom( self ):
        """ argsortRandom() -> [ int ], indices in random order."""
        pairs = [(random.random(), i) for i in range(0, len(self))]
        pairs.sort()
        return [ x[1] for x in pairs ]
        

    def argsort( self, sortKey ):
        """ argsort( str_sortKey ) -> [ int ], indices after sorting """
        pairs = [(self[i].info.get(sortKey), i) for i in range(0, len(self))]
        pairs.sort()
        return [ x[1] for x in pairs ]

    def sortBy( self, sortKey ):
        """
        sortBy( str_sortKey ) -> ComplexList (or sub-class) sorted by
        info[ str_sortKey ]
        """
        return self.take( self.argsort( sortKey ))


    def valuesOf(self, infoKey, default=None, indices=None, unique=0 ):
        """
        Get all values of a certain info record of all or some Complexes.
        infoKey - str, key for info dict
        default - any, default value if infoKey is not found (None)
        indices - list of int OR None(=all), indices of Complexes (None)
        unique  - 1|0, report each value only once (set union), (default 0)
        -> list of values
        """
        l = self
        if indices != None:
            l = N.take( l, indices )

        if not unique:
            return [ c.info.get(infoKey, default) for c in l ]

        r = []
        for c in l:
            if c.info.get(infoKey, default) not in r:
                r += [ c.info.get( infoKey ) ]
        return r


    def filterRange( self, infoKey, vLow, vHigh ):
        """
        filterRange( str_infoKey, vLow, vHigh )
        Get indices of Complexes where vLow <= c.info[ infoKey ] <= vHigh.
        -> array of int
        """
        vLst = self.valuesOf( infoKey )

        maskL = N.greater_equal( vLst, vLow )
        maskH = N.less_equal( vLst, vHigh )

        return N.nonzero( maskL * maskH )


    def filterEqual( self, infoKey, lst ):
        """
        filterEqual( infoKey, lst )
        Get indices of Complexes where c.info[ infoKey ] in lst.
        -> array of int
        """
        mask = [ c.info.get( infoKey ) in lst for c in self ]
        return N.nonzero( mask )
        
         
    def filterFunct( self, f ):
        """
        filterFunct( f )
        Get indices of Complexes where f( c ) == 1.
        -> array of int
        """
        mask = [ f( c ) for c in self ]
        return N.nonzero( mask )
        

    def filter( self, infoKey, cond ):
        """
        Complexes matching condition.
        infoKey - key of Complex.info dict (not used if cond is function )
        cond    - (vLow, vHigh) -> vLow <= c.info[ infoKey ] <= vHigh
                - list          -> c.info[ infoKey ] in cond
                - function      -> cond( c ) == 1
        -> ComplexList (or sub-class)
        !! ConditionError if cond is neither list nor tuple nor function
        """
        indices = None
        
        if type( cond ) == tuple:
            
            indices = self.filterRange( infoKey, cond[0], cond[1] )

        if type( cond ) == list:

            indices = self.filterEqual( infoKey, cond )

        if type( cond ) == types.FunctionType:

            indices = self.filterFunct( cond )
            
        if indices == None:
            try:
                indices = self.filterEqual( infoKey, [cond] )
            except:
                raise ConditionError( "Can't interprete filter condition.")

        return self.take(indices)


    def argmax( self, infoKey ):
        """
        infoKey - str
        -> int, index of complex c with highest c.infos[infokey] value
        """
        vLst = self.valuesOf( infoKey )
        return N.argmax( vLst )

    def max( self, infoKey ):
        """
        infoKey - str
        -> Complex, with highest c.info['infoKey'] value
        """
        return self[ self.argmax(infoKey) ]


    def argmin( self, infoKey ):
        """
        infoKey - str
        -> int, index of complex c with lowest c.infos[infokey] value
        """
        vLst = self.valuesOf( infoKey )
        return N.argmin( vLst )

    def min( self, infoKey ):
        """
        infoKey - str
        -> Complex, with lowest c.info['infoKey'] value
        """
        return self[ self.argmin( infoKey ) ]


    def getIndex( self, infoKey, value ):
        """
        -> int, list position of Complex where c.info['infoKey'] == value
        !! AmbiguousMatch, ComplexNotFound,
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
        -> Complex where c.info['infoKey'] == value
        !! AmbiguousMatch, ComplexNotFound,
           if there are more or less than 1 matches
        """
        return self[ self.getIndex( infoKey, value ) ]
    

    def take( self, indices ):
        """
        indices - array/list of int, list positions
        -> ComplexList with all items specified.
        """
        r = self.__class__( [ self[i] for i in indices ] )

        return r


    def toDict( self, infoKey ):
        """
        Convert list into dict indexed by a certain info-record value.
        infoKey - key of info dict in Complexes
        If several Complexes have the same value, the result will have
        a list registered for this key.
        
        EXAMPLE: clst.toDict('soln') -> {1:Complex, 3:Complex, solnN:Complex}

        -> dict, { info1:Complex, info2:Complex, info3:[Complex, Complex].. }
        """
        result = {}
        for c in self:
            t.dictAdd( result, c.info[infoKey], c )

        return result

    def toList( self ):
        """-> [ Complex ], simple python list of Complexes"""
        return list( self )

    def __maskNone( self, l1, l2 ):
        """
        take out positions from l1 and l2 that are None in either of them.
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
        Plot pairs of info values. The additional arg arguments are handed
        over to biggles.Points().
        """
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
        plot( xkey, [ykey1, ykey2..],[arg1=x, arg2=y]) -> biggles.FramedPlot
        Plot pairs of info values.
        """
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
        
    
#############################
## TEST #####################

if __name__ == '__main__':

    print "DOING something"

    l = t.Load(t.testRoot() + "/dock/hex01/complexes.cl")

    cl = ComplexList( l )

    p = cl.plot( 'rms', 'hex_eshape', 'hex_etotal' )
