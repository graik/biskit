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
from Biskit import ReduceCoordinates, EnsembleTraj

import time
import os.path

o = {'step':'1'}

def _use():
    print """
    Reduce all-atom trajectory to trajectory with only one backbone and
    up to 2 side chain atoms per residue.
    The new atoms are the centers of mass of several atoms and carry the
    weight of the pooled atoms in an atom profile called 'mass'.

    reduceTraj.py -i traj.dat [-o traj_reduced.dat -amber -red |red_traj.dat|]

         i     - pickled Trajectory object (all atoms)
         o     - alternative output file name (default: 'reduced_' + input.dat)
         red   - pickled reduced traj, just update ref model in given traj.
         amber - rename amber HIE/HIP/HID, CYX -> HIS, CYS; unwrap atom names
    
Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)

def renameAmberRes( model ):

    for a in model.getAtoms():
        if a['residue_name'] == 'CYX':
            a['residue_name'] = 'CYS'
        if a['residue_name'] in ['HIE','HID','HIP']:
            a['residue_name'] = 'HIS'

def unwrapAmberAtoms( model ):
    """ptraj puts last letter/number of 4-letter atom names first. Undo.
    could be avoided if ptraj would be told: trajout |file| pdb nowrap
    """
    numbers = map( str, range(10) )

    for a in model.getAtoms():
        if len(a['name'])==4 and a['name'][0] in numbers:
            a['name'] = a['name'][1:] + a['name'][0]


##########
## MAIN ##

start_time = time.time()

if len( sys.argv) < 2:
    _use()

o = cmdDict( o )

## default out name
if not 'o' in o:
    o['o'] = 'reduced_' + stripFilename( absfile(o['i'] ) ) + '.dat'
    o['o'] = os.path.dirname( absfile(o['i'])) + '/'+o['o']

print "Loading trajectory .. "
t = Load( absfile(o['i']) )

if 'amber' in o:
    renameAmberRes( t.getRef() )
    unwrapAmberAtoms( t.getRef() )


red = ReduceCoordinates( t.getRef() )
ref = red.reduceToModel( t.ref.getXyz() )


if 'red' in o:
    ## just replace ref model in already reduced trajectory
    tred = Load( absfile(o['red']))
    tred.setRef( ref )
    Dump( tred, absfile( o['red']) )
    sys.exit(0)

print "Reducing .. "
frames = red.reduceXyz( t.frames, axis=1 )

t.frames = frames
t.ref = ref

print "Dumping.."
Dump( t, absfile( o['o'] ) )

print "Done in %3.2f min" % ((time.time() - start_time)/60.0)
