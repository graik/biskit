#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## trajFluct.py
## Calculate global and side chain fluctuation per atom
## last $Author$
## last $Date$

from Biskit.tools import *
from Biskit import Trajectory, PDBDope
from Biskit.ProfileCollection import ProfileError

def _use( options ):
    print """
trajFluct.py: Calculate global and side chain fluctuation per atom
              for a trajectory.

Syntax:	   trajFluct -i trajectory_file [-o result_trajectory]

                     
Options:   -i     pickled trajectory
           -o     file name for pickled result Trajectory
           
"""
    for key, value in options.items():
        print "\t-",key, "\t",value
    
    sys.exit(0)


### MAIN ###

default = {'o':'traj_fluct.dat' }

if len (sys.argv) < 2:
    _use( default )

options = cmdDict( default )

traj = Load( options['i'] )

traj.fit( prof='rms2avg', ref=traj.ref )
traj.fit( prof='rms2avg' )

f_global = traj.getFluct_global()

f_local = traj.getFluct_local()

ref = traj.getRef()

## disconnect ref from its source
ref.disconnect()

ref.setAtomProfile( 'fluct_global', f_global,
                comment='fluctuation around average position in Å' )

ref.setAtomProfile( 'fluct_local', f_local,
                comment='fluctuation after fitting to each residue backbone')

try:

    doper = PDBDope( ref )

    flushPrint('adding accessible surface profiles...')
    doper.addASA()

except ProfileError, why:
    errWriteln("Couldn't add ASA, ProfileError: ", why)

Dump( traj, options['o'] )
