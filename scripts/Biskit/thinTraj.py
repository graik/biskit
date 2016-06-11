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
##

from Biskit.tools import *
from Biskit.EnsembleTraj import traj2ensemble

import os.path

o = {'step':'5'}

def _use():
    print """
    This script is used only for the test_multidock example and is
    used to remove frames from the test trajectory to speed up
    subsequent test steps. With he default setting of step=5 this
    will result in a 100 frame trajectory.

    thinTraj.py -i traj.dat [-step |int|] 

         i     - pickled Trajectory object
         step  - int, 1..keep all frames, 2..skip first and every second, ..
    
Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)



##########
## MAIN ##

if len( sys.argv) < 2:
    _use()

o = cmdDict( o )


print "Loading trajectory .. "
traj = load( absfile(o['i']) )
traj = traj2ensemble( traj )

print "Removing frames ..."
thin_traj = traj.thin( step=int(o['step']) )

print "Dumping thinned trajectory with %i frames (%i frames removed)..."\
      %(thin_traj.lenFrames(), traj.lenFrames()-thin_traj.lenFrames())
dump( thin_traj, absfile( o['i'] ) )

