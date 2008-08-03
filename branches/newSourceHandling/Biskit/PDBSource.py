##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2008-07-28 10:10:27 +0200 (Mon, 28 Jul 2008) $
## $Revision: 625 $

"""
Keep track of the source file/instance/url/code of a PDBModel.
"""

import os.path

import Biskit as B
import PDBParserFactory as PF

class PDBSourceError( Exception ):
    """Something is wrong with the source of a PDBModel."""
    pass

class PDBSource( object ):
    """
    Each PDBModel instance tracks its own source in order to avoid the
    redundant storage of xyz or profile data.
    
    This is a helper class for PDBModel and should rarely be used directly.
    """
    
    def __init__( self, source, check=True ):
        """
        @param source: PDBModel instance, or file or url location of parent 
                       model
        @type  source: str OR PDBModel or None or LocalPath
        """
        self._source = None
        if type( source ) == str and os.path.isfile( source ):
            self._source = B.LocalPath( source )
        else:
            self._source = source
        
        self.__parser = None
        self.__exists = None  # the source exists, None meaning not checked yet
        self.__ispermanent = None # the source is permanent, None .. not checked

    def getParser( self ):
        """
        @return: PDBParser that supports the given source
        @raise PDBSourceError, if there is no supporting parser found.
        """
        if self.__parser is None:
            if self._source is None:
                raise PDBSourceError, 'no source available'
            
            self.__parser = PF.PDBParserFactory.getParser( self._source )

        return self.__parser

    parser = property( getParser ) # gives PDBParser instance or PDBSourceError

    def exists( self ):
        """@return: True, if the source is available in memory or elsewhere."""
        if self.__exists is None:
            if self._source is None:
                self.__exists = False
            else:
                try:
                    self.__exists = self.parser.sourceExists( self._source )
                except PDBSourceError, e:
                    self.__exists = False
        return self.__exists
    
    def isPermanent( self ):
        """
        Check whether the source is in a permanent location (i.e. an
        existing file or url).
        """
        if self.__ispermanent is None:
            if self._source is None or type( self._source ) is B.PDBModel:
                self.__ispermanent = False
            else:
                self.__ispermanent = self.exists()
        
        return self.__ispermanent
    
    def isEmpty( self ):
        """
        @return: True, if there has been any source specified (existing or not)
        """
        return self._source is None
    
    def isInstance( self ):
        """@return: True, if the source is pointing to a PDBModel in memory"""
        return isinstance( self._source, B.PDBModel )
    
    
    def isPDBFile( self ):
        return isinstance( self.parser, PF.PDBParseFile )
    
    def getFile( self ):
        """
        @return: str, path to existing source file
        @raise: PDBSourceError if there is no existing source file
        """
        if type( self._source ) is B.LocalPath and self.exists():
            return self._source.local()
        elif isinstance( self._source, B.PDBModel):
            return self._source.source.getFile()
        
        raise PDBSourceError, 'no source file available for %r' % self._source
    
    file = property( getFile ) # path to existing source file or PDBSourceError

    def getLocalPath( self ):
        """
        @return: Biskit.LocalPath instance to existing or non-existing source
        @raise: PDBSourceError if source is not a file path
        """
        if type( self._source ) is B.LocalPath:
            return self._source
        elif isinstance( self._source, B.PDBModel):
            return self._source.source.getLocalPath()
        
        raise PDBSourceError, 'source %r is not a file path' % self._source

    
    def __repr__( self ):
        return 'PDBSource ' + repr( self._source)
    
    def __str__( self ):
        return str( self._source )

 #############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    def test_PDBPSource_file( self ):
        """PDBSource from file test"""
        f = os.path.join( B.tools.testRoot(), 'rec/1A2P_dry.model')
        self.s = PDBSource( f )
        self.assertTrue( type(self.s.parser) is PF.PDBParsePickle )
        self.assertTrue( self.s.exists() )
        self.assertTrue( self.s.isPermanent() )
        self.assertEqual( f, self.s.file )

    def test_PDBPSource_instance( self ):
        """PDBSource from instance test"""
        f = os.path.join( B.tools.testRoot(), 'rec/1A2P_dry.model')
        self.s = PDBSource( B.PDBModel(f) )
        self.assertTrue( type( self.s.parser ) is PF.PDBParseModel )
        self.assertTrue( self.s.exists() )
        self.assertFalse( self.s.isPermanent() )
        self.assertEqual( f, self.s.file )


if __name__ == '__main__':

    BT.localTest()
