#!/usr/bin/python
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
## $Revision$
## last $Date$
## last $Author$

#============================================
# Blast 2 sequences against each other,
# Return sequence identity
# largely obsolete - use BioPython instead
# Raik Gruenberg
#============================================

import tools as T       # errWriteln()
import commands         # getstatusoutput
import re               # reg. Expressions
import os		# os.remove()
import settings

##################################################
# methods

class Blast2Seq:
    """Determine sequence identity between 2 protein sequences"""
    
    def __init__(self):
        self._definePatterns()
        if self._fileInPath('bl2seq'):
            self.blastBin = 'bl2seq'
        else:
            self.blastBin = settings.bl2seq_bin
        self.blastout = ''
        self.inp1 = '/tmp/seq1.fasta'
        self.inp2 = '/tmp/seq2.fasta'


    def _fileInPath(self, fname):
        """Check if fname exists anywhere in system path."""

        paths = os.environ['PATH'].split( ":" )

        for p in paths:
            if os.path.exists( p + '/'+fname ):
                return 1

        return 0
    

    def _storeSequences(self, seq1, seq2):
        """create temporary files for bl2seq"""
        f1 = open(self.inp1, 'w')
        f2 = open(self.inp2, 'w')
        f1.write(">Sequence 1\n"+seq1)
        f2.write(">Sequence 2\n"+seq2)
        f1.close()
        f2.close()

    def _cleanup(self):
        """remove temporary files"""
        try:
            os.remove(self.inp1)
            os.remove(self.inp2)
        except:
            T.errWriteln("Blast2Seq._cleanup(): Error while removing temporary files.")
        
    def runBlast(self, seq1, seq2):
        """Perform Blast search of seq against other sequence."""
        self._storeSequences(seq1, seq2)    # create temp files
        # IDENTITY matrix supplied externally
        blastcmd = self.blastBin + ' -i ' + self.inp1 + ' -j ' + self.inp2 +\
                   ' -M BLOSUM62 -p blastp -F F '
        (status, blastout) = commands.getstatusoutput(blastcmd)
        self._cleanup()                     # remove temp files
        return  self.filterBlastHit(blastout)# Bsp: {'aln_id':1.0, 'aln_len':120}

            
    def _definePatterns(self):
        """Pre-compile RegEx patterns"""
        self.ex_identity = re.compile('.+ Identities = (\d+)/(\d+) ') # Blast Identities
        self.ex_expect = re.compile('.+ Expect = ([\d\-e\.]+)')

    def filterBlastHit( self, hitStr):
        """Extract sequence identity and overlap length from one single
        bl2seq hit"""
        hitStr = hitStr.replace( '\n', '#' ) # get rid of line breaks
        try:
            _id = self.ex_identity.match(hitStr)       # Blast Identity
            if (_id == None):
                return {'res_id':0, 'aln_id':0, 'aln_len':0}
            res_identical, res_number = int(_id.group(1)), int(_id.group(2))
            id = 1.0 * res_identical/res_number
            return {'res_id':res_identical, 'aln_id':id, 'aln_len':res_number}
        except:
            T.errWriteln("Error 1 in filterBlastHit:")
            T.errWriteln("Hit Parsed: " + hitStr)
            return {}


####################################
## TESTING
####################################

if __name__ == '__main__':
    blaster = Blast2Seq()
    print blaster.runBlast("AAAFDASEFFGIGHHSFKKEL","AAAFDASEFFGIGHHSAKK")
