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
## Convert single amber crd into Trajectory object

import sys
from Biskit.tools import *
from Biskit import AmberCrdParser

def _use():

    print """
Convert single amber crd into Trajectory object

amber2traj.py -i sim.crd -o traj_0.dat -r ref.pdb [-b -wat -hyd -rnres
              -code PDBC ]

    -i     input amber trajectory
    -o     output file with pickled biskit Trajectory object
    -r     reference PDB, must have identical atom content+order as sim.crd
    -b     traj has box info (3 additional coordinates per frame)
    -wat   delete WAT, Cl-, Na+ residues (after parsing)
    -hyd   delete all hydrogens (after parsing)
    -rnres rename amber residues HIE/HID/HIP, CYX to HIS and CYS
    -code  PDB code of molecule [first 4 letters of ref file name]
    """
    sys.exit( 0 )


if __name__ == '__main__':

    if len( sys.argv ) < 2:
        _use()

    o = cmdDict( {'o':'traj_0.dat', 'i':'sim.crd'} )
    fcrd = o['i']
    fpdb = o['r']
    fout = o['o']
    box  = o.has_key( 'b' )
    wat  = o.has_key('wat')
    hyd  = o.has_key('hyd')
    rnres  = o.has_key('rnres')
    code = o.get('code', None)

    p = AmberCrdParser( fcrd, fpdb, box, rnres, pdbCode=code )
    t = p.crd2traj()

    if wat:
        t.removeAtoms( lambda a: a['residue_name'] in ['WAT', 'Na+', 'Cl-'] )

    if hyd:
        t.ref.addChainId( keep_old=1 ) ## preserve chain-delimiters
        t.removeAtoms( lambda a: a['element'] == 'H' )

    print "Dumping result to ", fout
    dump( t, absfile(fout) )

    print "Done"
