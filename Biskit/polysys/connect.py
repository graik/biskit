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
from Biskit import PDBModel
import Biskit.tools as T

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

        if type( adapter ) is str and not os.path.isfile( adapter ) :
            self.adapter_model = PDBModel( self.adapters[ adapter ] )
        elif type( adapter ) is str and os.path.isfile( adapter ):
            self.adapter_model = PDBModel( adapter )
        elif isinstance( adapter, PDBModel ):
            self.adapter_model = adapter
        else:
            raise FusionError, 'invalid adapter argument: %r' % adapter
    
    def fuse( self, m1, m2 ):
        """
        """
        rmask = N.zeros( m1.lenResidues() )
        rmask[-1] = 1
        
        amask = N.zeros( m1.lenAtoms() )
        for atom in ['O', 'OT1', 'OT2', 'OXT']:
            amask += m1['name'] == atom

        amask = 
        
        return m1.concat( m2 )
