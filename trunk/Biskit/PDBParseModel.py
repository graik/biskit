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
## last $Author$
## last $Date$
## $Revision$
"""
Parse a in-memory PDBModel instance into a new PDBModel

@see L{PDBModel}
@see L{PDBParserFactory}
"""
import Numeric as N

import Biskit.tools as T
import Biskit as B
from PDBParser import PDBParser, PDBParserError


class PDBParseModel( PDBParser ):

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
        return isinstance( source, B.PDBModel )

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
        return 'in-memory instances of PDBModel'
        

    def update( self, model, source, skipRes=None, lookHarder=0 ):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.
        Override!
        @param model: existing model
        @type  model: PDBModel
        @param source: PDBModel object
        @type  source: str | file | PDBModel
        @param skipRes: list residue names that should not be parsed
        @type  skipRes: [ str ]
        """
        try:
            ## atoms and/or coordinates need to be updated from PDB
            if self.needsUpdate( model ):

                s = source

                model.fileName = model.fileName or s.fileName

                model.pdbCode  = model.pdbCode or s.pdbCode

                model.atoms = model.atoms or s.getAtoms()

                model.xyz = model.xyz or s.getXyz()

                model.__terAtoms = getattr(model, '_PDBModel__terAtoms',[])or \
                                   getattr(s,'_PDBModel__terAtoms',[])

                model.rProfiles.updateMissing( s.rProfiles,
                                               copyMissing=lookHarder)
                model.aProfiles.updateMissing( s.aProfiles,
                                               copyMissing=lookHarder)

                if skipRes:
                    model.removeRes( skipRes )
                               
        except B.ProfileError, why:
            EHandler.warning("Cannot read/update profiles from source: %r"\
                             %why)

        model.setSource( source )

class Test:
    """
    Test class
    """

    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: coordinates of center of mass
        @rtype:  array
        """

        ## loading output file from X-plor
        if local:
            print 'Loading pdb file ..'

        p = PDBParseModel()
        m = p.parse2new( B.PDBModel(T.testRoot()+'/rec/1A2P.pdb') )

        if local:
            globals().update( locals() )

        return N.sum( m.centerOfMass() )


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: coordinates of center of mass
        @rtype:  array
        """
        return N.sum( N.array([ 29.53385022,  46.39655482,  37.75218589]))
        

if __name__ == '__main__':

    test = Test()

    assert abs( test.run( local=1 ) - test.expected_result() ) < 1e-8
