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
List of ComplexEvolving instances.
"""

from Biskit.Dock.ComplexList import ComplexList, ComplexListError
from Biskit.Dock.ComplexEvolving import ComplexEvolving
import numpy as N

class ComplexEvolvingList( ComplexList ):
    """
    List of ComplexEvolving instances.
    Used for organising, sorting, and filtering Complexes during refinement.

    @todo: implement plotting functions for evolving Complexes
    @todo: right now normal Complexes are tolerated
    @todo: adapt model management
    @see: L{Dock.ComplexEvolving}
    """

    def __init__(self, lst=[] ):
        """
        @param lst: list of Complexes
        @type  lst: [ComplexEvolving]

        @raise ComplexListError: if list contains non-Complex item.
        """
        ComplexList.__init__( self, lst )


    def version( self ):
        """
        Version of Dock.Complex

        @return: version of class
        @rtype: str
        """        
        return 'ComplexEvolvingList $Revision$'


    def checkType( self, v ):
        """
        Make sure v is a ComplexEvolving

        @param v: any
        @type  v: any

        @raise ComplexListError: if list contains non-Complex item.
        """
        if not isinstance(v, ComplexEvolving):
            raise ComplexListError(
                str( v ) + " not allowed. ComplexList requires "+
                str(ComplexEvolving))


    def allVersionList( self ):
        """
        Get all versions of each Complex as a seperate Complex instance.

        @return: ComplexList of normal Complex instances
        @rtype: ComplexList
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

        @param version: version in history, -1 == last [-1] (default: -1)
        @type  version: int

        @return: ComplexList
        @rtype: ComplexList
        """
        return ComplexList( [ c[ version ] for c in self ] )


    def toList( self, version=None ):
        """
        Get a simple python list of Complexes. If version==None, the list
        contains ComplexEvolving instances with all versions, otherwise
        the list contains Complex instances representing a single version.

        @param version: version in history, -1 == last, None == all
                        (default: None)
        @type  version: int

        @return: python list of Complexes
        @rtype: [ Complex ]
        """
        if version is None:
            return ComplexList.toList( self )

        return [ c[ version ] for c in self ]


    def valuesOf( self, infoKey, version=None, default=None,
                  indices=None, unique=0 ):
        """
        Get all values of a certain info record of all or some Complexes.

        @param infoKey: key for info dict
        @type  infoKey: str
        @param version: index in history or None (=current) (default: None)
        @type  version: int
        @param default: default value if infoKey is not found (default: None)
        @type  default: any
        @param indices: list of int OR None(=all), indices of Complexes
                        (default: None)
        @type  indices: [int] OR None
        @param unique: report each value only once (set union), (default: 0)
        @type  unique: 1|0

        @return: list of values
        @rtype: [any]
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


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_ComplexEvolvingList(self):
        """Dock.ComplexEvolvingList test"""

        import Biskit.tools as t
        from Biskit.Dock import ComplexEvolving
        from Biskit.Dock import ComplexEvolvingList

        ## original complex
        cl = t.load(  t.testRoot() + "/dock/hex/complexes.cl" )

        ## first evolution step
        c = ComplexEvolving( cl[0].rec(), cl[0].lig(), cl[0],
                             info={'comment':'test1'})

        ## second evolution step
        c = ComplexEvolving( c.rec(), c.lig(), c,
                             info={'comment':'test2'})

        ## create an evolving complex list
        cl = ComplexEvolvingList( [c, c] )

        if self.local:
            ## last version of all complexes in list
            print cl.valuesOf('comment')

            ## version 1
            print cl.valuesOf('comment', version=1)

            ## the first complex in the list
            print cl[0].valuesOf('comment')

            globals().update( locals() )

        self.assertEqual( (cl.valuesOf('comment'), cl[0].valuesOf('comment')),
                          (['test2', 'test2'], [None, 'test1', 'test2']) )

if __name__ == '__main__':

    BT.localTest()


