#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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

## last $Author$
## last $Date$
## $Revision$

from Biskit import PDBModel
from Biskit.tools import *
from numpy.oldnumeric import *

def __pairwiseDistances( u, v):
    """
    pairwise distance between 2 3-D numpy arrays of atom coordinates.
    -> Numpy array len(u) x len(v)
    donated by Wolfgang Rieping.
    """
    diag1= diagonal(dot(u,transpose(u)))
    diag2= diagonal(dot(v,transpose(v)))
    dist= -dot(v,transpose(u))-transpose(dot(u,transpose(v)))
    dist= transpose(asarray(map(lambda column,a:column+a, \
                               transpose(dist), diag1)))

    return transpose(sqrt(asarray(map(lambda row,a: row+a, dist, diag2))))


def getSS( model, cutoff=4.0 ):

    cys_mask = model.mask( lambda a: a['residue_name'] in ['CYS', 'CYX']\
                           and a['name'] == 'SG')

    model.compress( cys_mask )

    diff = __pairwiseDistances( model.xyz, model.xyz )

    count = sum( sum( less( diff, cutoff ) ) ) - model.lenAtoms()

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
