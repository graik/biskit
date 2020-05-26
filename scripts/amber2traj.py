#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2019 Raik Gruenberg
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
## Convert single amber crd into Trajectory object

import sys
from biskit.tools import *
from biskit.md import Trajectory
import numpy as N

def _use():

    print("""
Convert single amber crd into Trajectory object

amber2traj.py -i sim.crd -o traj_0.dat -r ref.pdb [-b -wat -hyd -rnres
              -code PDBC -step 10 ]

    -i     input amber trajectory
    -o     output file with pickled biskit Trajectory object
    -r     reference PDB, must have identical atom content+order as sim.crd
    -b     traj has box info (3 additional coordinates per frame)
    -wat   delete WAT, Cl-, Na+ residues (during parsing)
    -hyd   delete all hydrogens (after parsing)
    -rnres rename amber residues HIE/HID/HIP, CYX to HIS and CYS
    -code  PDB code of molecule [else first 4 letters of ref file name]
    -step  only take every [step]th frame 
    """)
    sys.exit( 0 )


if __name__ == '__main__':

    if len( sys.argv ) < 2:
        _use()

    o = cmdDict( {'o':'traj.dat', 'i':'sim.crd'} )
    fcrd = o['i']
    fpdb = o['r']
    fout = o['o']
    box  = 'b' in o
    wat  = 'wat' in o
    hyd  = 'hyd' in o
    rnres  = 'rnres' in o
    code = o.get('code', None)
    step = int(o.get('step', 1))

    t = Trajectory( fcrd, fpdb, hasbox=box, rmwat=wat, verbose=True)
    
    if code:
        t.ref.pdbCode = code

    if hyd:
        t.ref.addChainId( keep_old=1 ) ## preserve chain-delimiters
        t.removeAtoms( lambda a: a['element'] == 'H' )

    if step != 1:
        t = t.takeFrames( N.arange(0,len(t),step) )

    print("Dumping result to ", fout)
    dump( t, absfile(fout) )

    print("Done")
