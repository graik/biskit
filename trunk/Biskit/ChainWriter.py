## ChainWriter:
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
## Write cleaned peptide_chains as PDB for import into XPlor
##
## $Version: $
## last $Date$
## last $Author$

from Scientific.IO.PDB import *
import commands, os

import tools as T

class ChainWriter:
    """
    Take chain from chainCleaner; write single PDB with
    unique segementID and no chainID
    """

    def __init__(self, path):
        """
        Take chains from ChainCleaner and write pdb files.
        File names are created from segid of each chain + '_seg.pdb'
        path - string, output path for PDB files
        """
        self.path = T.absfile( path )

    def _startPDB(self, chain, fname):
        """
        Create pdb file and write header.
        chain - Scientific.IO.PDB.PeptideChain
        fname - string, file name
        -> PDBFile, handle of open file
        """
        f = PDBFile(fname, 'w')
        try:
            for c in chain.comments:    # write comments, if any
                f.writeComment(c)
        except:
            pass
        return f

    def removeTER(self, fname):
        """
        Remove TER record from PDB.
        fname - string, name of existing file.
        """
        try:
            command = 'egrep -v "^TER " ' + fname + '> temp.pdb'
            commands.getstatusoutput(command)
            os.rename('temp.pdb', fname)
        except (OSError):
            T.errWriteln("Error removing 'TER' statement from %s: ")
            T.errWriteln( T.lastError() )


    def writeChain(self, chain):
        """
        Write single chain as PDB. File name will be segid + '_seg.pdb'.
        chain - Scientific.IO.PDB.PeptideChain
        """
        if (chain <> None):
            try:

                fname = self.path + '/' + chain.segment_id + "_seg.PDB"
                file = self._startPDB(chain, fname)     # include comments
                chain.writeToFile(file)     # create PDB
                file.close()                # make sure file is complete!
                self.removeTER(file.file.file.name)   # remove TER record from PDB
                return 1        # wrote a chain

            except (IOError):
                T.errWriteln("Error writing chain to file %s:" % fname)
                T.errWriteln( T.lastError() )
                
        return 0            # false, no more chains to write

##############################
## TESTING
##############################

if __name__ == '__main__':

    from ChainCleaner import ChainCleaner
    from ChainSeparator import ChainSeparator

    cleaner = ChainCleaner( ChainSeparator( \
        T.testRoot()+'/com_wet/1BGS_edited.pdb',
        T.testRoot()+'/com_wet') )
    
    writer = ChainWriter( T.testRoot()+'/com_wet' )
    print 'Writing separated, cleaned chains to disk...'
    print writer.writeChain( cleaner.next() )
    print writer.writeChain( cleaner.next() )
    print writer.writeChain( cleaner.next() )

    print "DONE."

    
