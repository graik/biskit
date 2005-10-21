#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## convert one normal Trajectory into EnsembleTraj
##
## last $Author$
## last $Date$
## $Revision$

import os, sys

import Biskit.tools as T
from Biskit import Trajectory
from Biskit import EnsembleTraj
from Biskit.EnsembleTraj import traj2ensemble

def use():
    if len( sys.argv ) < 2:
        print \
"""
Convert one normal Trajectory into EnsembleTraj. The input trajectory must
have frame names that allow sorting by time and ensemble member (see
EnsembleTraj.py for details).

traj2ensemble.py -i |in_traj| -n |n_members| -o |out_traj| -pdb |PDBCode| ]

    o        - out file name        (default: replace input file)
    n        - number of ensemble members to expect (default: 10)
    pdb      - PDB code to be stored in trajectory
"""
        sys.exit(0)

##########
## MAIN ##

use()

o = T.cmdDict( {'n':10} )

f_in  = T.absfile( o['i'] )
f_out = T.absfile( o.get('o', f_in) )
n = int( o['n'] )

T.flushPrint("Loading...")
t = T.Load( f_in )

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
T.Dump( t, f_out )

T.flushPrint("Done.\n")
