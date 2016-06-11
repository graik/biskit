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


from Biskit import PDBModel
from Biskit.tools import *
import numpy as N
import Biskit.mathUtils as MU


def getSS( model, cutoff=4.0 ):

    cys_mask = model.mask( lambda a: a['residue_name'] in ['CYS', 'CYX']\
                           and a['name'] == 'SG')

    model = model.compress( cys_mask )

    diff = MU.pairwiseDistances( model.xyz, model.xyz )

    count = N.sum( N.sum( N.less( diff, cutoff ) ) ) - model.lenAtoms()

    return count / 2


if __name__ == '__main__':

    if len( sys.argv ) < 2:
        print """Count number of SS bonds in protein.
        Syntax: getSS.py |input1.pdb| |input_2.pdb| ..
        """
        sys.exit(0)

    for f in sys.argv[1:]:
        
        m = PDBModel( f )
        
        flushPrint( '%s : %2i\n' % (m.sourceFile(),getSS( m )) )
