##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
"""
@see L{Biskit.PDBParseFile}
@see L{Biskit.PDBParseModel}
"""
import Biskit
import Biskit.tools as T
from PDBParseFile   import PDBParseFile

class BioUnitError( Biskit.BiskitError ):
    pass

class BioUnit:
    """
    A class for assembling multimers according to BIOMT data
    """
    
    def __init__( self, model ):
        """
        @param model: PDBModel with BIOMT record
        @type  model: Biskit.PDBModel
        """
        self.biomt = model.info['BIOMT']
        self.model = model

    def makeMultimer (self, biomoleculeNum):
        """
        @param biomoleculeNum: ID of the biomolecule (from BIOMT record)
        @type  biomoleculeNum: int
        @return PDBModel, with the bio-molecule as specified in BIOMT
        """
        try:
            (chainIds, transformMatrices) = self.biomt[biomoleculeNum]
        except:
            raise BioUnitError, \
                  'This unit does not have biomolecule #%i' % biomoleculeNum
            ## return model

        mask = self.model.maskFrom( 'chain_id', chainIds )
        chains = self.model.compress( mask )
        monomers = [ chains.transform ( rt ) for rt in transformMatrices ]
        result = monomers[0].concat( *monomers[1:] )
        return result

    def moleculeList (self):
        return self.biomt.keys()

#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_BioUnit( self ):
        """BioUnit test"""

        if self.local:
            print 'Loading pdb file ..'

        self.p = PDBParseFile()
        self.m = self.p.parse2new( T.testRoot('biounit/2V4E.pdb') )
        self.m.report()

        self.b = BioUnit(self.m)
        if self.local:
            print 'unit has', len(self.b.moleculeList()), 'molecules'

        self.mul = self.b.makeMultimer(1)
        self.mul.report()

if __name__ == '__main__':

    BT.localTest()

    ## test concat performance
    import Biskit.tools as T

    chains = [ m.takeChains( [i] ) for i in range( 16 ) ]

    T.profile( 'm = chains[0].concat( *chains[1:] )' )
    
