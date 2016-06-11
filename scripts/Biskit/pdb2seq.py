#!/usr/bin/env python
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
