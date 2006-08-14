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

"""
Blast 2 sequences against each other,
Return sequence identity

@note: argely obsolete - use BioPython instead
"""

import tools as T       # errWriteln()
import commands         # getstatusoutput
import re               # reg. Expressions
import os		# os.remove()
import settings

##################################################
# methods

class Blast2Seq:
    """
    Determine sequence identity between 2 protein sequences
    """

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
        """
        Check if fname exists anywhere in system path.

        @param fname: filename
        @type  fname: str

	@return: file status
	@rtype: 1|0
        """
        paths = os.environ['PATH'].split( ":" )

        for p in paths:
            if os.path.exists( p + '/'+fname ):
                return 1

        return 0


    def _storeSequences(self, seq1, seq2):
        """
        create temporary files for bl2seq

	@param seq1: sequence string
	@type  seq1: str
	@param seq2: sequence string
	@type  seq2: str	
        """
        f1 = open(self.inp1, 'w')
        f2 = open(self.inp2, 'w')
        f1.write(">Sequence 1\n"+seq1)
        f2.write(">Sequence 2\n"+seq2)
        f1.close()
        f2.close()


    def _cleanup(self):
        """
        remove temporary files
        """
        try:
            os.remove(self.inp1)
            os.remove(self.inp2)
        except:
            T.errWriteln("Blast2Seq._cleanup(): Error while removing temporary files.")


    def runBlast(self, seq1, seq2):
        """
        Perform Blast search of seq against other sequence.

	@param seq1: sequence string
	@type  seq1: str
	@param seq2: sequence string
	@type  seq2: str

	@return: blast results, e.g. {'aln_id':1.0, 'aln_len':120}
	@rtype: dict
        """
	# create temp files
        self._storeSequences(seq1, seq2)
	
        # IDENTITY matrix supplied externally
        blastcmd = self.blastBin + ' -i ' + self.inp1 + ' -j ' + self.inp2 +\
                   ' -M BLOSUM62 -p blastp -F F '
        (status, blastout) = commands.getstatusoutput(blastcmd)
	
	# remove temp files
        self._cleanup()
	
        return  self.filterBlastHit(blastout)


    def _definePatterns(self):
        """
        Pre-compile RegEx patterns 
        """
	# Blast Identities and Expext value
        self.ex_identity = re.compile('.+ Identities = (\d+)/(\d+) ') 
        self.ex_expect = re.compile('.+ Expect = ([\d\-e\.]+)')


    def filterBlastHit( self, hitStr):
        """
        Extract sequence identity and overlap length from one single
        bl2seq hit

	@param hitStr: blast result file
	@type  hitStr: str

	@return: blast results, e.g. {'aln_id':1.0, 'aln_len':120}
	@rtype: dict
	"""
	# get rid of line breaks
        hitStr = hitStr.replace( '\n', '#' )
	
        try:
            _id = self.ex_identity.match(hitStr)  # Blast Identity
            if (_id == None):
                return {'res_id':0, 'aln_id':0, 'aln_len':0}
            res_identical, res_number = int(_id.group(1)), int(_id.group(2))
            id = 1.0 * res_identical/res_number
            return {'res_id':res_identical, 'aln_id':id, 'aln_len':res_number}
        except:
            T.errWriteln("Error 1 in filterBlastHit:")
            T.errWriteln("Hit Parsed: " + hitStr)
            return {}


#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """

    def run( self ):
        """
        run function test

        @return: balast alignment result
        @rtype: dict
        """
        blaster = Blast2Seq()
        
        return blaster.runBlast("AAAFDASEFFGIGHHSFKKEL",
                                "AAAFDASEFFGIGHHSAKK")


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: alignment data
        @rtype:  str
        """
        return {'aln_len': 19, 'aln_id': 0.94736842105263153, 'res_id': 18}

            

if __name__ == '__main__':

    test = Test()

    result = test.run()
    print result
    
    assert result == test.expected_result()





