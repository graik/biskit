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

from PDBParseFile   import PDBParseFile

class BioUnit:
    """
    A class for assembling multimers according to BIOMT data
    """
    def __init__( self, model ):
        self.biomt = model.info['BIOMT']

    def makeMultimer (self, model, biomoleculeNum):
        try:
            (chainIds, transformMatrices) = self.biomt[biomoleculeNum]
        except:
            print 'This unit does not have biomolecule #', (biomoleculeNum)
            return model

        mask = model.maskFrom( 'chain_id', chainIds )
        chains = model.compress( mask )
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

        ## loading output file from X-plor
        if self.local:
            print 'Loading pdb file ..'

        self.p = PDBParseFile()
        self.m = self.p.parse2new( '../../2V4E.pdb')
        self.m.report()

        self.b = BioUnit(self.m)
        print 'unit has', len(self.b.moleculeList()), 'molecules'

        self.mul = self.b.makeMultimer(self.m, 1)
        self.mul.report()

if __name__ == '__main__':

    BT.localTest()
