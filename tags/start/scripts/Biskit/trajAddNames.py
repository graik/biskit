#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## Add file names of frames to trajectory
## last $Author$
## last $Date$

from Biskit.tools import *
from Biskit import Trajectory

if __name__ == '__main__':

    if len( sys.argv ) < 2:
        print \
"""
trajAddNames.py -i |in_traj.| -o |out_traj.| -f |file1 file2 file..|

Add/replace file names of frames to existing (pickled) Trajectory.
"""
        sys.exit(0)
        
    options = cmdDict({})

    print "Loading ", options['i']
    
    traj = Load( options['i'] )

    if len( traj.frames ) != len( options['f']):
        print "ERROR: frame number != file name number"
        sys.exit(1)

    ## strip of folders etc.

    options['f'] = map( lambda s: stripFilename( s ), options['f'] )

    traj.frameNames = options['f']

    print "Dumping ", options['o']

    Dump( traj, options['o'] )

    print "Done"
    
