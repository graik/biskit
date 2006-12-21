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
## last $Author$
## last $Date$
## $Revision$
"""
utility funtions for Mod package
"""

import os.path
import types
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
                fname, chain_id = l.split()

                if not os.path.exists(fname):
                    fname = '%s/%s'%(T.testRoot(), fname)
                    
                if not len(fname) == 0:
                    result[ fname ] = chain_id
                    
            except:
                fname = l.strip()
                
                if not os.path.exists(fname):
                    fname = '%s/%s'%(T.testRoot(), fname)                
                    
                result[ fname ] = ''

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


def verify_fasta( target ):
    """
    Verify that a given file or string is in Fasta format.
    The definition used for a fasta file here is that:
     - first line starts with '>'
     - the following sequence lines are not longer that 80 characters
     - the characters has to belong to the standard amino acid codes

    @param target: name of fasta file OR file contents as list of strings
    @type  target: str OR [str]
    
    @return: conforms to the fsata format
    @rtype: True/False
    """
    if not type(target) == types.ListType:
        if os.path.exists( target ):
            f = open( target, 'r' )
            target = f.readlines()

    if not target[0][0] == '>':
        print 'Fasta format does not contain description line.'
        return False

    for i in range( 1, len(target) ):
        if len( target[i] ) >= 80:
            print 'Fasta sequence lines longer that 80 characters'
            return False

        for j in target[i]:
            aa_codes = MU.aaDicStandard.values() + [ '\n' ]
            if not j.upper() in aa_codes:
                print 'Invalid amino acid code: %s'%j.upper()
                return False

    return True




#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype:  int
        """
        vf = verify_fasta( T.testRoot() + '/Mod/project/target.fasta')

        if local:
            globals().update( locals() )

        return vf


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    

if __name__ == '__main__':

    test = Test()
    
    assert test.run( local=1 ) ==  test.expected_result()


