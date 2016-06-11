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

## pdb2traj.py
## Collect coordinate frames of a molecule

import Biskit.tools as tools
from Biskit import Trajectory, EnsembleTraj

import sys, os

def _use( options ):
    print """
pdb2traj.py: Collect many coordinate frames ( pdb or pickled PDBModel ) of one
             molecule. Write Trajectory object. Waters are removed.

Syntax:    pdb2traj -i pdb1 pdb2 ..  [ -e -r |ref_structure| -o |out_file| -f
                    -wat -c ]
OR         pdb2traj -d folder/with/pdbs/or/models [ -r ... ]

Options:   -i     input pdb files or pickled PDBModel objects
           -d     folder containing input pdb or pickled PDBModel files
           -e     create EnsembleTraj, input files must be ordered first
                  by time then by member; x_30_10.pdb sorts before x_100_09.pdb
           -r     reference structure with PDB records (incl. waters),
                  if not given, the first file from -i is used
           -wat   delete TIP3, HOH, Cl-, Na+ from ref and all frames
           -o     file name for resulting pickled Trajectory object
           -f     fit to reference (dry reference if given)
           -c     analyze atom content of all files seperately before casting
                  them to reference. Default: only analyze first file in -i.

Note about reference structure: The atom order and content of the files given
with -i is adapted to the order/content of the reference PDB but NOT
vice-versa. Snapshots can hence have additional atoms (which are removed) but
they must have, at least, all the atoms that are in the reference. 

Default options:
"""
    for key, value in options.items():
        print "\t-",key, "\t",value

    sys.exit(0)

### MAIN ###
############

if len (sys.argv) < 3:
    _use( {'o':'traj.dat'} )

options = tools.cmdDict( {'o':'traj.dat'} )

## get all PDBs and models directly from a directory (avoids shell limits)
if 'd' in options:
    d = tools.absfile( options['d'] )
    l = os.listdir( d )
    l = [ x for x in l if x[-4:].upper() == '.PDB' or x[-6:] == '.MODEL' \
          or x[-7:].upper() == '.PDB.GZ' ]
    l = [ os.path.join( d, f ) for f in l ]

    options['i'] = l

if 'e' in options:
    t_class = EnsembleTraj
else:
    t_class = Trajectory

traj = t_class( options['i'],
                options.get('r', None),
                rmwat   =options.has_key('wat'),
                castAll =options.has_key('c') )

## remove dependencies to source
traj.ref.disconnect()

if options.has_key('f'):
    traj.fit()

tools.dump( traj, options['o'] )
