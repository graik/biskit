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
## Analyze entropy of a single trajectory with ptraj


import sys

import Biskit.tools as t
from Biskit.AmberEntropist import AmberEntropist
from Biskit import LogFile

def _use():
    print """
Analyze entropy of a single or one rec and one lig trajectory with ptraj.

Syntax:  a_trajEntropy.py -i |traj1.dat+traj2.dat| [ -o |result.dic|
                      -ref |ref_structure| -cast -chains |chain_indices|
                      -border |chain| -split -shift -shuffle
                      -s |startFrame| -e |endFrame| -step |frameOffset|
                      -ss |member_startFrame| -se |member_endFrame|
                      -ex_n |exclude_n_members|
                      -ex1 |ex_from_traj1| -ex2 |ex_from_traj2|
                      -ex3 |exclude_member_tripple|
                      -atoms |CA CB ..| -heavy -protein
                      -nice |level|
                      -parm |parm_file| -crd |crd_file| -f_out |ptraj_out|
                      -f_template |ptraj_template|
                      -log |log_file| -debug -verbose ]
Options:
        i          1 trajectory or 2 trajectories connected by '+'
        o          file name for pickled result dictionary
        s          skip first |s| frames (of complete trajectory)
        e          skip frames after |e| (of complete trajectory)
        ss         skip first |ss| frames of each member trajectory
        se         skip frames after |se| of each member trajectory
        atoms      considered atoms (default: all)
        heavy      remove hydrogens (default: don't)
        solvent    retain solvent atoms and ions (default: don't)
        protein    remove non-protein atoms (default: don't)
        ref        pickled PDBModel, Complex, or Trajectory
        cast       equalize atom content of traj and ref         [no]
        chains     list of integer chain indices e.g -chains 0 1 [all]
        border     1st chain of 2nd molecule for -split, -shift, -shuffle
        split      split complex trajectory and fit rec and lig separately
                   (requires -border with first lig chain)       [no]
        shuffle    shuffle the order of rec vs. lig frames
        thin       use randomly distributed fraction of frames, e.g. 0.2 [all]
        step       frame offset, use every step frame, e.g. 5            [all]
        ex1        exclude these members from 1st trajectory, e.g. 3 6
        ex2        exclude these members from 2nd trajectory (if given)
        ex_n       exclude first n members                       [0]
        ex3        exclude |ex3|rd tripple of members, e.g. 2 excludes 3,4,5
                   (0 excludes nothing)                          [0]

        f_template alternative ptraj input template [default template]
        f_out      target name for ptraj output file        [discard]
        nice       nice level                                     [0]
        log        file for program log                       [STOUT]
        debug      keep all temporary files
        verbose    print extended progress messages to log    [log != STDOUT]

Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def defOptions():
    o = {}
    o['step'] = 1
    o['o']    = 'entropy_result.dic'

    ## uncomment for testing
##     o['step'] = 2
##     o['start']= 100
##     o['stop'] = 400
##     o['ref']  = '~/interfaces/c15/com_wet/ref.complex'
##     o['i'] = '~/interfaces/c15/rec_pcr_00/traj.dat+\
##               ~/interfaces/c15/lig_pcr_00/traj.dat'
    return o

##########
## MAIN ##


options = t.cmdDict( defOptions() )

if not ('i' in options ):
    _use()

## convert command line parameters
for k in ['s','e','ss', 'se', 'step','fit_s', 'fit_e', 'nice','thin', 'ex_n',
          'ex3', 'shift', 'border']:
    if k in options:
        options[k] = int( options[k] )

for k in ['debug', 'verbose', 'cast', 'split', 'shuffle', 'shift', 'heavy',
          'solvent', 'protein']:
    if k in options:
        options[k] = 1
    else:
        options[k] = 0

for k in ['chains','ex1', 'ex2', 'ex']:
    if k in options:
        options[k] = t.toIntList( options[k] )

if 'atoms' in options:
    options['atoms'] = t.toList( options['atoms'] )

if 'ex1' in options and 'ex2' in options:
    options['ex'] = ( options['ex1'], options['ex2'] )
else:
    options['ex'] = options.get( 'ex', options.get('ex1', None) )

if 'log' in options:
    options['log'] = LogFile( options['log'] )

f_in = options['i']
del options['i']

a = AmberEntropist( f_in, **options )
a.run()

t.dump( a.result, options['o'] )

print "Dumped detailed result to %s. (for python unpickling)" % options['o']
print "Entropy in cal/mol-kelvin (total, vibrational): ",
print a.result['S_total'], ',', a.result['S_vibes']
print
