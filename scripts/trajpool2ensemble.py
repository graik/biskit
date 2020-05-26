#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2019 Raik Gruenberg & Johan Leckner
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
## convert one normal Trajectory into EnsembleTraj
##

import os, sys

import biskit.tools as T
from biskit.md import Trajectory
from biskit.md import EnsembleTraj
from biskit.md.ensembleTraj import traj2ensemble

def use():
    if len( sys.argv ) < 2:
        print("""
Convert one normal Trajectory into EnsembleTraj. The input trajectory must
have frame names that allow sorting by time and ensemble member (see
EnsembleTraj.py for details).

traj2ensemble.py -i |in_traj| -n |n_members| -o |out_traj| -pdb |PDBCode| ]

    o        - out file name        (default: replace input file)
    n        - number of ensemble members to expect (default: 10)
    pdb      - PDB code to be stored in trajectory
""")
        sys.exit(0)

##########
## MAIN ##

use()

o = T.cmdDict( {'n':10} )

f_in  = T.absfile( o['i'] )
f_out = T.absfile( o.get('o', f_in) )
n = int( o['n'] )

T.flushPrint("Loading...")
t = T.load( f_in )

T.flushPrint("Converting %i frames..." % len(t) )

if isinstance(t, EnsembleTraj ):
    T.flushPrint( "Nothing to be done!\n")
    sys.exit(0)
    
t = traj2ensemble( t, n )
if 'pdb' in o:
    t.ref.pdbCode = o['pdb']

if f_in == f_out:
    os.rename( f_in, f_in + '_backup')

T.flushPrint("Saving...")
T.dump( t, f_out )

T.flushPrint("Done.\n")
