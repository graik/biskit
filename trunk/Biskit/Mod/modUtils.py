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
"""utility funtions for Mod package"""

import os.path
import types
import Biskit.molUtils as MU

def parse_tabbed_file( fname ):
    """
    fname - str, file name of key : value mapping
    -> { key : value, }
    """
    f = open( fname )

    result = {}
    for l in f:
        if not l[0] == '#':

            try:
                fname, chain_id = l.split()
                if not len(fname) == 0:
                    result[ fname ] = chain_id
            except:
                fname = l.strip()
                result[ fname ] = ''

    f.close()

    return result


def format_fasta(seq, width=60):
    """
    Transform a given sequence to fasta format
    seq -str, sequence
    -> str, string sequence in fasta format
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
    -> True/False
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
