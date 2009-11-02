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
    Generate structure for a stretch of amino acids, loosely extended along 
    x-axis. Xtender wraps the XPlor genExtendedStructure command. The resulting
    chain is not fully extended and may still contain some vdw clashes --
    the default 5000 minimization steps resolve most but not necessarily all
    of these clashes. (Uses xplor in Python mode)
    
    Requires: Xplor-NIH
    """

    xplor_script = """
import psfGen, protocol, pdbTool
import protocol

simWorld.setRandomSeed(5521)
#~ protocol.topology['protein'] = '/home/victor/xplor-nih-2.23/toppar/topallh6x.pro' 
#~ protocol.topology['protein'] = '/home/victor/xplor-nih-2.23/toppar/charmm22/topallh22x.pro' 
#~ protocol.parameters['protein'] = '/home/victor/xplor-nih-2.23/toppar/paramallh3x.pro'
#~ protocol.parameters['protein'] = '/home/victor/xplor-nih-2.23/toppar/charmm22/parallh22x.lip'
protocol.initTopology(('protein'))
protocol.initParams(('protein'))

seq = "%(seq_long)s"
psfGen.seqToPSF(seq,seqType='prot',startResid=1, segName='')

protocol.genExtendedStructure(maxFixupIters=50)

protocol.fixupCovalentGeom('all')

xplor.command( "minimize powell nstep=%(min_cycles)i nprint 10 end")

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

    x = Xtender( 10 * 'AGTSYTWDS', min_cycles=5000 )
    model = x.run()
    
    x2 = Xtender( 3 * 'GSGSGSGSGSGSG', min_cycles=5000 )
    model2 = x2.run()

    
    
