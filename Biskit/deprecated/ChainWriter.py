## ChainWriter:
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
Write cleaned peptide_chains as PDB for import into XPlor

This is vintage code. See L{Biskit.PDBCleaner} for a more recent
version. The writing of single chains could easily be done with PDBModel.

@todo: re-implement pdb2xplor with PDBCleaner / PDBModel
"""
from Scientific.IO.PDB import *
import commands
import os.path

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

        @param path: output path for PDB files
        @type  path: string
        """
        self.path = T.absfile( path )


    def _startPDB(self, chain, fname):
        """
        Create pdb file and write header.

        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: chain object
        @param fname: file name
        @type  fname: string

        @return: handle of open file
        @rtype: PDBFile
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

        @param fname: name of existing file.
        @type  fname: string
        """
        try:
            #command = 'egrep -v "^TER " ' + fname + '> temp.pdb'
            path = os.path.dirname( fname )
            command = 'egrep -v "^TER " %s > %s/temp.pdb'%( fname, path )
            commands.getstatusoutput(command)
            os.rename('%s/temp.pdb'%path, fname)
        except (OSError):
            T.errWriteln("Error removing 'TER' statement from %s: ")
            T.errWriteln( T.lastError() )


    def writeChain(self, chain):
        """
        Write single chain as PDB. File name will be segid + '_seg.pdb'.

        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: chain object
        """
        if (chain <> None):
            try:

                fname = self.path + '/' + chain.segment_id + "_seg.PDB"
                file = self._startPDB(chain, fname)     # include comments
                chain.writeToFile(file)     # create PDB
                file.close()                # make sure file is complete!
                self.removeTER(file.file.file.name) #remove TER record from PDB
                return 1        # write a chain

            except (IOError):
                T.errWriteln("Error writing chain to file %s:" % fname)
                T.errWriteln( T.lastError() )

        return 0            # false, no more chains to write


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test ChainWriter"""   

    def prepare(self):
        self.fname =   T.testRoot() + '/rec/1A2P_rec_original.pdb'
        self.outPath = T.tempDir()

    def test_ChainWriter( self ):
        """ChainWriter test"""

        from ChainCleaner import ChainCleaner
        from ChainSeparator import ChainSeparator

        self.cleaner = ChainCleaner( ChainSeparator( self.fname,
                                                     self.outPath ) )

        self.writer = ChainWriter( self.outPath )

        all_msg = []
        for i in range(3):
            msg = self.writer.writeChain( self.cleaner.next() )
            if self.local:
                print 'Writing separated, cleaned chain to disk...%i'%msg
            all_msg += [ msg ]

        self.assertEquals( all_msg, [1, 0, 0] )

    def cleanUp(self):
        T.tryRemove( self.outPath + '/1A2P_waters.pdb' )
        T.tryRemove( self.outPath + '/1A2A_seg.PDB' )



if __name__ == '__main__':

    BT.localTest()

