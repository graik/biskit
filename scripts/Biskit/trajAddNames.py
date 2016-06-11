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

## Add file names of frames to trajectory

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
    
    traj = load( options['i'] )

    if len( traj.frames ) != len( options['f']):
        print "ERROR: frame number != file name number"
        sys.exit(1)

    ## strip of folders etc.

    options['f'] = map( lambda s: stripFilename( s ), options['f'] )

    traj.frameNames = options['f']

    print "Dumping ", options['o']

    dump( traj, options['o'] )

    print "Done"
    
