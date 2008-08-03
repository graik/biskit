## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

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
## last $Author$
## last $Date$
## $Revision$
"""
Parse a pickled PDBModel from disc into a new PDBModel instance 
"""
import os
import Biskit.tools as T
import Biskit as B

from PDBParser import PDBParserError
from PDBParseModel import PDBParseModel

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
        
        @return: True if the given source is supported by this parser
                 implementation (equivalent to isinstance( source, PDBModel) )
        @rtype: bool
        """
        return ((type(source) is str) or isinstance(source, B.LocalPath)) \
               and os.path.isfile( str( source ) )

    @staticmethod
    def description():
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.description():
        >>>     ...

        @return: short free text description of the supported format
        @rtype: str
        """
        return 'pickled PDBModel (file)'

    
    @staticmethod
    def sourceExists( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.sourceExists('~/myfile.model'):
        >>>     ...

        @return: whether or not the given source is available for reading
        @rtype: bool
        """
        if type( source ) is str:
            return os.path.exists( T.absfile( source ) )
        if type( source ) is B.LocalPath:
            return source.exists()
        
        raise PDBParserError, 'unsupported source type: %r' % source


    def update( self, model, source, skipRes=None, updateMissing=0, force=0 ):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.
        Profiles that are taken from the source are labeled 'changed'=0.
        The same holds for coordinates (xyzChanged=0).
        However, existing profiles or coordinates or fields remain untouched.

        @param model: existing model
        @type  model: PDBModel
        @param source: source PDB file or pickled PDBModel or PDBModel object
        @type  source: str || file || PDBModel
        @param skipRes: list residue names that should not be parsed
        @type  skipRes: [ str ]
        @param updateMissing: check source for additional profiles [0] 
        @type  updateMissing: 1|0
        """
        try:
            if force or updateMissing or self.needsUpdate( model ):

                if type( source ) is not str:
                    source = source.local()
                s = T.Load( source )

                super( PDBParsePickle, self ).update(
                    model, s, skipRes=skipRes, updateMissing=updateMissing,
                    force=force )
                              
        except Exception, why:
            print T.lastErrorTrace()
            raise PDBParserError, "Cannot unpickle source model from %s, "\
                   % str(source) + "Reason:\n" + str(why)

        model.setSource( source )


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_PDBParsePickle( self ):
        """PDBParsePickle test"""

        import numpy.oldnumeric as N

        ## loading output file from X-plor
        if self.local:
            print 'Loading pickled model ..'

        self.p = PDBParsePickle()
        self.m = self.p.parse2new( T.testRoot()+'/rec/1A2P_dry.model')

        self.assertAlmostEqual( N.sum( self.m.centerOfMass() ),
                                 114.18036673321811, 7)

if __name__ == '__main__':
   
    BT.localTest()
