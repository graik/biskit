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
Parse a pickled PDBModel from disc into a new PDBModel instance 
"""
import biskit.tools as T
import biskit as B

from biskit.core.pdbparser import PDBParserError
from biskit.core.pdbparseModel import PDBParseModel
from biskit.core.localpath import LocalPath

class PDBParsePickle( PDBParseModel ):
    """
    Parse a pickled PDBModel from disc into a new PDBModel instance 
    """

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.supports( model ):
        >>>     ...
        
        :return: True if the given source is supported by this parser
                 implementation (equivalent to isinstance( source, PDBModel) )
        :rtype: bool
        """
        return (type(source) is str) or isinstance(source, LocalPath)

    @staticmethod
    def description():
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.description():
        >>>     ...

        :return: short free text description of the supported format
        :rtype: str
        """
        return 'pickled PDBModel (file)'
        

    def update( self, model, source, skipRes=None, updateMissing=0, force=0,
                headPatterns=[] ):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.
        Profiles that are taken from the source are labeled 'changed'=0.
        The same holds for coordinates (xyzChanged=0).
        However, existing profiles or coordinates or fields remain untouched.

        :param model: existing model
        :type  model: PDBModel
        :param source: source PDB file or pickled PDBModel or PDBModel object
        :type  source: str || file || PDBModel
        :param skipRes: list residue names that should not be parsed
        :type  skipRes: [ str ]
        :param updateMissing: check source for additional profiles [0] 
        :type  updateMissing: 1|0
        """
        try:
            if force or updateMissing or self.needsUpdate( model ):

                s = T.load( source )

                super( PDBParsePickle, self ).update(
                    model, s, skipRes=skipRes, updateMissing=updateMissing,
                    force=force )
                              
        except Exception as why:
            raise PDBParserError("Cannot unpickle source model from %s, "\
                   % str(source) + "Reason:\n" + str(why))

        model.setSource( source )


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_PDBParsePickle( self ):
        """PDBParsePickle test"""

        import biskit.core.oldnumeric as N0

        ## loading output file from X-plor
        if self.local:
            print('Loading pickled model ..')

        self.p = PDBParsePickle()
        self.m = self.p.parse2new( T.testRoot('rec/1A2P_dry.model'))

        self.assertAlmostEqual( N0.sum( self.m.centerOfMass() ),
                                114.18037, 5)

if __name__ == '__main__':
   
    BT.localTest()
