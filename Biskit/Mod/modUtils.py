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
utility funtions for Mod package
"""

import os.path as osp
import Biskit.molUtils as MU
import Biskit.tools as T

def parse_tabbed_file( fname ):
    """
    Parse the chaim index file written by TemplateSearcher
    L{TemplateSearcher.F_CHAIN_INDEX}.

    @param fname: name of file to parse 
    @type  fname: str

    @return: key : value mapping
    @rtype: {key:value}
    """
    f = open( fname )

    result = {}
    for l in f:
        if not l[0] == '#':

            try:
                fpdb, chain_id = l.split()

                if not osp.exists(fpdb):
                    fpdb = osp.join( osp.split( fname )[0], fpdb )
##                     fpdb = '%s/%s'%(T.testRoot(), fpdb)

                if not len(fpdb) == 0:
                    result[ fpdb ] = chain_id

            except:
                ## no chain given
                fpdb = l.strip()

                if not osp.exists(fpdb):
                    fpdb = osp.join( osp.split( fname )[0], fpdb )

                result[ fpdb ] = ''

    f.close()

    return result


def format_fasta(seq, width=60):
    """
    Transform a given sequence to fasta format

    @param seq: sequence
    @type  seq: str
    @param width: length of a line in characters (default: 60)
    @type  width: int

    @return: string sequence in fasta format
    @rtype: str
    """
    fasta_sequence = ""

    for i in xrange(0,len(seq),width):            
        fasta_sequence += seq[i:i+width]

        if(i+width>=len(seq)):
            pass
        else:
            fasta_sequence += "\n"

    return fasta_sequence




#############
##  TESTING        
#############
import Biskit.test as BT    

class Test(BT.BiskitTest):
    """
    Test class
    """

    def test_modUtiles(self):
        """Mod.modUtils test"""
        vf = format_fasta( 'ABC'*100 )

        self.assertEqual( vf, self.EXPECT )

    EXPECT = """ABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABC
ABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABC
ABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABC
ABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABC
ABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABCABC"""


if __name__ == '__main__':

    BT.localTest()
