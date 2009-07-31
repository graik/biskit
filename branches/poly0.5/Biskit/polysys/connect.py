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
##
## last $Author: graik $
## last $Date: 2009-07-02 17:36:36 +0200 (Thu, 02 Jul 2009) $
## $Revision: 810 $

import os.path as osp
from Biskit import PDBModel, BiskitError, PDBError
import Biskit.tools as T
from Biskit.Residue import Residue ## experimental
import numpy as N

class FusionError( BiskitError ):
    pass


class Fusion( object ):
    
    adapters = \
        { 'helix'   : T.dataRoot('polysys/fusion/adapter_helix.model'),
          'sheet'   : T.dataRoot('polysys/fusion/adapter_sheet.model'),
          'extended': T.dataRoot('polysys/fusion/adapter_extended.model')
          }
    
    adapter_root = T.dataRoot('/polysys/fusion/')
    
    def __init__( self, adapter='extended' ):
        """
        Create new fusion object.
        @param adapter: adapter type to use for fusion peptide bond orientation.
                        (one of: 'helix', 'sheet', 'extended'(default)) or
                        custom PDBModel with two residues.
        @type  adapter: str OR PDBModel
        """
        try:
            if type( adapter ) is str and not osp.isfile( adapter ) :
                self.adapter_model = PDBModel( self.adapters[ adapter ] )
            elif type( adapter ) is str and osp.isfile( adapter ):
                self.adapter_model = PDBModel( adapter )
            elif isinstance( adapter, PDBModel ):
                self.adapter_model = adapter
            else:
                raise FusionError, 'invalid adapter argument: %r' % adapter

        except PDBError, why:
            raise FusionError, 'cannot load adapter: ' + why
    
    def select_atoms( self, m, rindex=0, names=[''] ):
        """
        Select atoms within a single residue
        """
        rmask = N.zeros( m.lenResidues() )
        rmask[rindex] = 1
        
        amask = N.zeros( m.lenAtoms() )
        for atom in names:
            amask += m['name'] == atom

        return amask * m.res2atomMask( rmask )
       
        
    def remove_atoms( self, m, rindex=-1, names=['O', 'OT1', 'OT2', 'OXT'] ):
        """
        Remove atoms with given name from given residue.
        """     
        amask = self.select_atoms( m, rindex, names )
        m.remove( N.flatnonzero( amask ) )

        
    def fuse( self, m1, m2 ):
        """
        """
        m1 = m1.clone()
        m2 = m2.clone()
        
        self.remove_atoms( m1, rindex=-1, names=['O', 'OT1', 'OT2', 'OXT'] )
        self.remove_atoms( m2, rindex=0,  names=['HT1', 'HT2', 'HT3'] )
        
        ## fit adapter to C-terminal of model 1
        mask_1  = self.select_atoms( m1, rindex=-1, names=['N', 'CA', 'C'] )
        mask_adapt = self.select_atoms( self.adapter_model, 
                                   rindex=0, names=['N', 'CA', 'C'] )
        
        ref = m1.compress( mask_1 )
        
        assert( isinstance( self.adapter_model, PDBModel ) ) 
        self.adapter_model.fit( ref, mask=mask_adapt ) 
        
        ## fit N terminal model2 to C terminal of adapter        
        mask_2 = self.select_atoms( m2, rindex=0, names=['N', 'CA', 'C'] )
        mask_adapt = self.select_atoms( self.adapter_model, 
                                        rindex=1, names=['N', 'CA', 'C'] )
        
        ref = self.adapter_model.compress( mask_adapt )
        
        m2.fit( ref, mask=mask_2 )
        
        ## concatenate models
        assert( isinstance( m1, PDBModel ) )
        r = m1.concat( m2 )
        r['chain_id'] = r['chain_id'][0] * len( r )
        
        return r
    
if __name__ == '__main__':
    
    m1 = PDBModel( T.testRoot( 'polysys/spectrinr16.pdb' ) )
    
    f = Fusion( adapter='helix' )
    r = f.fuse( m1, m1 )