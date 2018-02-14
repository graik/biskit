#!/usr/bin/python
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
##

"""
Blast 2 sequences against each other,
Return sequence identity

@note: argueably obsolete - use BioPython instead
"""

from Biskit import Executor, TemplateError
import tools as T          
import tempfile, re
import os.path
from Biskit.Errors import BiskitError


class Blast2SeqError( BiskitError ):
    pass


class Blast2Seq( Executor ):
    """
    Determine sequence identity between 2 protein sequences
    """

    def __init__(self, seq1, seq2, **kw ):
        """
        @param seq1: sequence string
        @type  seq1: str
        @param seq2: sequence string
        @type  seq2: str
        """
        self.seq1 = seq1
        self.seq2 = seq2

        self.inp1 = tempfile.mktemp('_seq1.fasta')
        self.inp2 = tempfile.mktemp('_seq2.fasta')

        # Blast Identities and Expext value
        self.ex_identity = re.compile('.+ Identities = (\d+)/(\d+) ') 
        self.ex_expect = re.compile('.+ Expect = ([\d\-e\.]+)')

        blastcmd = '-i %s -j %s  -M BLOSUM62 -p blastp -F F'\
                 %(self.inp1, self.inp2)

        Executor.__init__( self, 'bl2seq', blastcmd, catch_out=1, **kw )


    def prepare( self ):
        """
        create temporary fasta files for bl2seq
        """
        f1 = open(self.inp1, 'w')
        f2 = open(self.inp2, 'w')
        f1.write(">Sequence 1\n"+self.seq1)
        f2.write(">Sequence 2\n"+self.seq2)
        f1.close()
        f2.close()


    def cleanup(self):
        """
        remove temporary files
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.inp1 )
            T.tryRemove( self.inp2 )


    def filterBlastHit( self ):
        """
        Extract sequence identity and overlap length from one single
        bl2seq hit

        @return: blast results, e.g. {'aln_id':1.0, 'aln_len':120}
        @rtype: dict
        """
        ## check that the outfut file is there and seems valid
        if not os.path.exists( self.f_out ):
            raise Blast2SeqError,\
                  'Hmmersearch result file %s does not exist.'%self.f_out

        if T.fileLength( self.f_out ) < 10:
            raise Blast2SeqError,\
                  'Hmmersearch result file %s seems incomplete.'%self.f_out

        out = open( self.f_out, 'r' )
        hitStr = out.read()
        out.close()

        # get rid of line breaks
        hitStr = hitStr.replace( '\n', '#' )

        try:
            _id = self.ex_identity.match(hitStr)  # Blast Identity
            if (_id is None):
                return {'res_id':0, 'aln_id':0, 'aln_len':0}
            res_identical, res_number = int(_id.group(1)), int(_id.group(2))
            id = 1.0 * res_identical/res_number
            return {'res_id':res_identical, 'aln_id':id, 'aln_len':res_number}
        except:
            T.errWriteln("Error 1 in filterBlastHit:")
            T.errWriteln("Hit Parsed: " + hitStr)
            return {}


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.filterBlastHit( )


#############
##  TESTING        
#############

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test Blast2Seq"""

    TAGS = [BT.EXE]

    def test_all( self ):
        """Blast2Seq test """
        self.blaster = Blast2Seq("AAAFDASEFFGIGHHSFKKEL",
                                 "AAAFDASEFFGIGHHSAKK") 

        self.r = self.blaster.run()

        if self.local:
            print self.r

        self.assertEqual( self.r,{'aln_len': 19, 'aln_id': 0.94736842105263153,
                                  'res_id': 18} )


if __name__ == '__main__':

    BT.localTest()



