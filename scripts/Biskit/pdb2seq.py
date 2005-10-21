#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

import sys

sys.path.insert(0, '/home/Bis/raik/data/tb/dock/scripts/modules')

from Biskit.tools import *
from Biskit import PDBModel

def _use():
    print """
    Extract AA sequence from PDB.

    Syntax pdb2seq.py |pdb_file|
    """

def printSequence( header, s ):

    n_chunks = len( s ) / 80

    print ">%s" % header 

    for i in range(0, n_chunks+1):

        if i * 80 + 80 < len( s ):
            chunk = s[i * 80 : i * 80 + 80]
        else:
            chunk = s[i * 80 :]

        print chunk


if __name__ == '__main__':

    if len( sys.argv ) < 1:
        _use()

    m = PDBModel( sys.argv[1] )

    seq = m.sequence()

    printSequence( m.pdbCode, seq )
