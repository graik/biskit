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
import numpy as N

import Biskit as B
import Biskit.tools as T

class BioUnitError( B.BiskitError ):
    pass

class BioUnit:
    """
    A class for assembling multimers according to BIOMT data
    """
    
    def __init__( self, model, biomt_dict ):
        """
        @param model: parent PDBModel
        @type  model: Biskit.PDBModel
        """
        self.model = model
        self.postprocess(biomt_dict)

    def postprocess(self, biomt_dict):
        """
        postprocess BIOMT data
        """
        self.biomol = {}
        for biomol_num in biomt_dict.keys():
            (chain_ids, rt_matrices) = biomt_dict[biomol_num]
            biomol_key = 'biomol'+str(biomol_num)
            self.biomol[biomol_key] = rt_matrices
            atom_mask = self.model.maskFrom( 'chain_id', chain_ids )
            self.model.residues[biomol_key] = self.model.atom2resMask( atom_mask )

        
    def makeMultimer (self, biomol_id=None):
        """
        @param biomol_id: ID of the biomolecule
        @type  biomol_id: string
        @return PDBModel, with the bio-molecule as specified in BIOMT
        """
        try:
            biomol_id = biomol_id or self.biomol.keys()[0]
            rt_matrices = self.biomol[biomol_id]
        except:
            raise BioUnitError, \
                  'This unit does not have biomolecule ' % biomol_id

        try:
            atom_mask = self.model.res2atomMask ( self.model[biomol_id] )
            chains = self.model.compress( atom_mask )
        except:
            raise BioUnitError, \
                  'Parent model does not have biomolecule ' % biomol_id

        #mask = self.model.maskFrom( 'chain_id', chainIds )
        #chains = self.model.compress( mask )
        monomers = [ chains.transform ( rt ) for rt in rt_matrices ]
        result = monomers[0].concat( *monomers[1:] )
        return result

    def moleculeIndex (self):
        """
        @return string list, the ids for multimer atom masks
        """
        return self.biomol.keys()

    def keys():
        """instead of moleculeIndex() ??"""
        return self.biomol.keys()

    def __len__(self):
        """Return number of biomolecules as length"""
        return len(self.biomol)

    def __getitem__(self, key):
        """biounit[key] ==> biounit.biomol[key]"""
        return self.biomol[ key ]


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

        self.p = B.PDBParseFile.PDBParseFile()
        
        self.m = self.p.parse2new( T.testRoot('biounit/2V4E.pdb') )
        self.m.report()
        
        if self.local:
            print 'unit has', len(self.m.biounit.moleculeIndex()), 'multimers'

        self.mul = self.m.biounit.makeMultimer('biomol1')
        self.mul.report()
        
        ##2V4E has 2 disjoint multimers
        
        ##self.sanitycheck = self.mul.biounit.makeMultimer('biomol2')

if __name__ == '__main__':

    BT.localTest()

    ## test concat performance
    ## import Biskit.tools as T

    ## chains = [ m.takeChains( [i] ) for i in range( 16 ) ]

    ## T.profile( 'm = chains[0].concat( *chains[1:] )' )
    
