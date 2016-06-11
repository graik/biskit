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

## trajFluct.py
## Calculate global and side chain fluctuation per atom

from Biskit.tools import *
from Biskit import Trajectory, PDBDope
from Biskit.ProfileCollection import ProfileError

def _use( options ):
    print """
trajFluct.py: Calculate global and side chain fluctuation per atom
              for a trajectory.

Syntax:  trajFluct -i trajectory_file [-o result_trajectory]

                     
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

traj = load( options['i'] )

traj.fit( prof='rms2avg', ref=traj.ref )
traj.fit( prof='rms2avg' )

f_global = traj.getFluct_global()

f_local = traj.getFluct_local()

ref = traj.getRef()

## disconnect ref from its source
ref.disconnect()

ref.atoms.set( 'fluct_global', f_global,
                comment='fluctuation around average position in A' )

ref.atoms.set( 'fluct_local', f_local,
                comment='fluctuation after fitting to each residue backbone')

try:

    doper = PDBDope( ref )

    flushPrint('adding accessible surface profiles...')
    doper.addSurfaceRacer()

except Exception, why:
    errWriteln("Couldn't add surface area profiles, ProfileError: ", why)

dump( traj, options['o'] )
