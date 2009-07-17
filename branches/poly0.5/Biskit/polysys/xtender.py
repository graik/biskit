##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner

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

## last $Revision: 680 $
## last $Date: 2009-04-25 21:06:41 +0200 (Sat, 25 Apr 2009) $
## last $Author: graik $
"""Generate peptide chain structure"""

import Biskit.molUtils as MU
import Biskit.tools as t
from Biskit import Xplorer, PDBModel

import tempfile

class Xtender( Xplorer ):
    """
    Generate structure for a loosely extended stretch of amino acids along x-axis.
    (Uses xplor in Python mode)
    """

    xplor_script = """
import psfGen, protocol, pdbTool
import protocol

simWorld.setRandomSeed(5521)

protocol.initTopology(('protein'))
protocol.initParams(('protein'))

seq = "%(seq_long)s"
psfGen.seqToPSF(seq,seqType='prot',startResid=1, segName='')

protocol.genExtendedStructure(maxFixupIters=50)

xplor.command( "minimize powell nstep=%(min_cycles)i nprint 5 end")

pdbTool.PDBTool("%(f_pdb)s").write()
"""

    def __init__( self, sequence='GSGSGS', f_pdb='', min_cycles=5000, **kw ):
        """
        @param sequence: str, single-letter code amino acid sequence of the linker
        @param f_pdb:    str, file name for result pdb (otherwise deleted)
        @param min_cycles: int, number of minimization steps
        """
        Xplorer.__init__( self, template=self.xplor_script, args='-py', **kw )

        self.seq_short = sequence
        self.seq_long  = ' '.join( MU.single2longAA( sequence ) )
        
        self.f_pdb = f_pdb or tempfile.mktemp( '_extender.pdb' )
        self.remove_pdb = not (self.f_pdb)

        self.min_cycles = min_cycles


    def cleanup( self ):

        Xplorer.cleanup( self )
        if not (self.debug or self.remove_pdb):
            t.tryRemove( self.f_pdb )

    def finish( self ):
        self.result = PDBModel( self.f_pdb )
        self.result.disconnect()


if __name__ == '__main__':

    x = Xtender( 'AGTSYTWDS', min_cycles=5000 )
    model = x.run()
    

    
    
