#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## $Date$
## $Revision$

from Biskit.tools import *
from Biskit import EnsembleTraj

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
traj = Load( absfile(o['i']) )

print "Removing frames ..."
thin_traj = traj.thin( step=int(o['step']) )

print "Dumping thinned trajectory with %i frames (%i frames removed)..."\
      %(thin_traj.lenFrames(), traj.lenFrames()-thin_traj.lenFrames())
Dump( thin_traj, absfile( o['i'] ) )

