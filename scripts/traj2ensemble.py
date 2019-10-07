#!/usr/bin/env python
## numpy-oldnumeric calls replaced by custom script; 09/06/2016
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
## pool several trajectory objects to one ensemble trajectory
##

import sys

import biskit.tools as T
from biskit.md import Trajectory
from biskit.md import EnsembleTraj
from biskit.md.ensembleTraj import traj2ensemble

import biskit.core.oldnumeric as N0

def use():
    if len( sys.argv ) < 2:
        print("""
Pool several trajectory objects to one ensemble trajectory.
Each sub-trajectory is considered as traj of one ensemble member.
This script is ignoring any profiles of the given trajectories and
re-assigns new frame names ( the Trajectory.concat() method is not
used to allow handling of larger trajectories )

traj2ensemble.py -i |in_traj1 in_traj2 ..| -o |out_traj|
                [-s |start_frame(incl.)| -e |end_frame(excl.)| -step |step|
                 -ref |ref_pdb_or_model| -pdb |PDBCode| -prot ]

    s,e,step - start, end position and stepping for each of the input traject.
    ref      - PDB or pickled PDBModel with reference coordinates, by default,
               the reference frame of the first trajectory is taken
    pdb      - PDB code to be stored in trajectory
    prot     - delete all non-protein atoms (not by default)
""")
        sys.exit(0)


def addFrameNames( traj, trajNumber ):

    traj.frameNames = []

    for i in range( len( traj ) ):
        traj.frameNames += ['t%04i_ens%02i' % (i, trajNumber ) ]


def loadTraj( f, trajIndex, start=0, end=None, step=1, prot=False ):
    """Load traj from file, add frame names, extract portion if requested"""
    
    t = T.load( T.absfile( f ) )
    addFrameNames( t, trajIndex )

    e = end or len( t )

    if start or end or (step != 1):
        t = t.takeFrames( list(range( start, e, step)) )

    if prot:
        t.keepAtoms( N0.nonzero(t.ref.maskProtein()) )
        
    return t


############
### MAIN ###

use()

o = T.cmdDict( {'o':'traj_ensemble.dat'} )

out = T.absfile( o['o'] )
inLst = T.toList( o['i'] )
start = int( o.get('s','0') )
end   = o.get( 'e', None )
if end:
    end = int( end )
step = int( o.get('step',1) )
prot = 'prot' in o

ref = o.get('ref',None)
if ref:
    ref = PDBModel( T.absfile( ref ) )
    if 'prot' in o:
        ref = ref.compress( ref.maskProtein() )


T.flushPrint("Loading and appending trajectories...")
r = loadTraj( inLst[0], 0, start, end, step, prot=prot )

for i in range( 1, len( inLst ) ):

    t = loadTraj( inLst[i], i, start, end, step )

    if t.ref.equals( r.ref ) != [1,1]:
        raise Exception( 'Incompatible reference structure.' )

    r = r.concat( t )
    
    T.flushPrint('#')
    
print(" Done")

r.ref.disconnect()

if 'pdb' in o:
    r.ref.pdbCode = o['pdb']


## T.flushPrint("Converting to EnsembleTraj...")
## r = traj2ensemble( r, len(inLst))

## T.flushPrint( ("Done\nDumping %i member ensemble traj to "% result.n_members) + o['o'] )
## T.dump( result, T.absfile( o['o'] ) )
