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
.. seealso:: `biskit.core.PDBParser`,`biskit.core.PDBParseFile`,
             `biskit.core.PDBParseModel`, `biskit.core.PDBParsePickle`
"""

from biskit.core.pdbparseFile   import PDBParseFile
from biskit.core.pdbparseModel  import PDBParseModel
from biskit.core.pdbparsePickle import PDBParsePickle
from biskit.core.pdbparseNCBI   import PDBParseNCBI
from biskit.core.pdbparser import PDBParserError

class PDBParserFactory:
    """
    Provide the right PDBParser for different structure sources.
    """

    @staticmethod
    def getParser( source ):
        """
        getParser( source ) -> PDBParser; Fetch a Parser for the source.

        The method is static and should be called directly with the class::

           p = ParserFactory.getParser( 'myfile.pdb' )

        :param source: structure source (PDB file, PDBModel, pickled model)
        :type source: str | LocalPath | PDBModel

        :return: a parser that should be able to handle the given source
        :rtype: PDBParser (child)

        :raise PDBError: if no compatible parser is found
        """

        if PDBParseFile.supports( source ):
            return PDBParseFile()

        if PDBParseModel.supports( source ):
            return PDBParseModel()

        if PDBParseNCBI.supports( source ):
            return PDBParseNCBI()

        if PDBParsePickle.supports( source ):
            return PDBParsePickle()
            
        raise PDBParserError('Format of %r is not recognized.' % source)

#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """nothing to test"""
    pass
