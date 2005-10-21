## organise, sort, and filter Complexes in refinement
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

from ComplexList import ComplexList
from Complex import Complex as ProtComplex
from ComplexEvolving import ComplexEvolving

class ComplexEvolvingList( ComplexList ):
    """
    List of ComplexEvolving instances.
    @todo implement plotting functions for evolving Complexes
    @todo right now normal Complexes are tolerated
    @todo adapt model management
    @see ComplexEvolving
    """

    def __init__(self, lst=[] ):
        """
        lst - list of Complexes
        !! raise ComplexListError, if list contains non-Complex item.
        """
        ComplexList.__init__( self, lst )
        

    def version( self ):
        return 'ComplexEvolvingList $Revision$'


    def checkType( self, v ):
        """Make sure v is a ComplexEvolving"""
        if not isinstance(v, ComplexEvolving):
            raise ComplexListError(
                str( v ) + " not allowed. ComplexList requires "+
                str(ComplexEvolving))
        

    def allVersionList( self ):
        """
        Get all versions of each Complex as a seperate Complex instance.
        -> ComplexList of normal Complex instances
        """
        r = ComplexList()

        for c in self:
            try:
                r += c.toList()
            except:
                r += [ c ]
            
        return r


    def toComplexList( self, version=-1 ):
        """
        Get a ComplexList that contains only a single version of each Complex.
        version - int, version in history, -1 == last [-1]
        -> ComplexList
        """
        return ComplexList( [ c[ version ] for c in self ] )


    def toList( self, version=None ):
        """
        Get a simple python list of Complexes. If version==None, the list
        contains ComlexEvolving instances with all versions, otherwise
        the list contains Complex instances representing a single version.
        version - int, version in history, -1 == last, None == all [None]
        -> [ Complex ]
        """
        if version is None:
            return ComplexList.toList( self )

        return [ c[ version ] for c in self ]


    def valuesOf( self, infoKey, version=None, default=None,
                  indices=None, unique=0 ):
        """
        Get all values of a certain info record of all or some Complexes.
        infoKey - str, key for info dict
        version - int, index in history or None (=current) [None]
        default - any, default value if infoKey is not found [None]
        indices - list of int OR None(=all), indices of Complexes [None]
        unique  - 1|0, report each value only once (set union), [0]
        -> list of values
        """
        l = self
        if indices != None:
            l = N.take( l, indices )

        if not unique:
            if version is None:
                return [ c.get(infoKey, default) for c in l ]
            return [ c[version].get( infoKey, default) for c in l ]

        r = []
        for c in l:
            if version is not None:
                c = c[ version ]

            if c.info.get(infoKey, default) not in r:
                r += [ c.info.get( infoKey ) ]

        return r

    
### TEST ####

if __name__ == '__main__':

    import Biskit.tools as t
    
    l = t.Load( "~/interfaces/c15/dock_xray/hex01/complexes.cl")

    cl = ComplexList( l )

    c = ComplexEvolving( cl[0].rec(), cl[0].lig(), cl[0],
                         info={'comment':'test1'})
    c = ComplexEvolving( c.rec(), c.lig(), c,
                         info={'comment':'test2'})

    cl = ComplexEvolvingList( [c, c] )
