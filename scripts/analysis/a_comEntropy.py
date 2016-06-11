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
## Calculate vibrational entropy of binding from trajectory of rec, lig, and
## complex


import sys, time

import Biskit.tools as t
import Biskit.mathUtils as mu
import Biskit.PVM.hosts as hosts
from Biskit.AmberEntropyMaster import AmberEntropyMaster

def _use():
    print """
    Run many AmberEntropist calculations on many nodes. The Master has
    a standard set of 13 protocols to run on rec, lig, and com
    trajectories, as well as on every single member trajectory - in
    total 113.  It accepts one variable parameter, e.g. s(tart). Each
    protocol is then run for all values of the variable parameter.
    The script puts many temporary trajectories into the folder with the
    input trajectories -- consider creating a new folder for each trajectory!

Syntax:  a_comEntropy.py -rec |rec.traj| -lig |lig.traj| -com |com.traj|
                         -out |out.dat| [ -cr |rec_chains| -zfilter |cutoff|
                         -s |from| -e |to| -ss |from| -se |to|
                         -thin |fraction| -step |offset|
                         -var |option| -vrange |v1 v2..| -jack
                         -exrec |members| -exlig |members| -excom |members|
                         -hosts |name1 name2..| -clean  -single ]

Options:
        rec    - str, free rec trajectory
        lig    - str, free lig trajectory
        com    - str, complex trajectory
        out    - str, file name for pickled result
        cr     - [int], chains of receptor in complex trajectory [n_chains rec]
        var    - str, name of variable option [ s ]
        vrange - [any], set of values used for variable option
                 OR 'start:stop:step' i.e string convertable to arange() input
        jack   - set up leave-one-trajectory-out jackknife test [don't]
                 (replaces var with 'ex1' and vrange with range(1,n_members+1))

        zfilter- float, kick out outlyer trajectories using z-score threshold
                 [None->don't]
        exrec  - [int], exclude certain members of receptor ensemble    [[]]
        exlig  - [int], exclude certain members of ligand  ensemble     [[]]
        excom  - [int], exclude certain members of complex ensemble     [[]]

        clean  - remove pickled ref models and member trajectories [0]
        hosts  - [str], nodes to be used [all known]
        h      - int, number of nodes to be used from all known [all]
        single - run only one job on multi-processor nodes [0]
        mem    - float, run only on machines with more than |mem| GB RAM [0]
        debug  - don't delete output files [0]

        ... parameters for AmberEntropist -- can also be given as -var
        cast    - equalize free and bound atom content [1]
        s,e     - int, start and stop frame                 [0, to end]
        ss, se  - int, start and stop frame of single member trajectories
                  (only works with EnsembleTraj; overrides s,e)
        atoms   - [ str ], names of atoms to consider       [all]
        heavy   - delete all hydrogens                      [don't]
        protein - delete non-protein atoms                  [don't]
        step    - int, frame offset                         [no offset]
        thin    - float, use randomly distributed fraction of frames [all]
                  (similar to step but sometimes better)
        all     - only calculate with all members, no single member values
        ex      - [int] exclude same members from rec, lig and com 
        ex_n    - int, exclude last n members  OR...                 [0]
        ex3     - int, exclude |ex3|rd tripple of trajectories       [0]
                  (0 excludes nothing, 1 excludes [0,1,2] )
        ex1     - int, exclude ex1-th member remaining after applying ex [None]
                  (0 excludes nothing, 1 excludes [0] )

        ... parameters for AmberCrdEntropist, Executor, Master
        f_template - str, alternative ptraj input template  [default]
        verbose    - print progress messages to log     [log != STDOUT]
        w          - show x-windows   [no]
        a          - 0|1, add hosts to PVM [1]

Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def defOptions():
    o = {}
    return o

##########
## MAIN ##

options = t.cmdDict( defOptions() )

if len( options ) < 4:
    _use()

## convert command line parameters
for k in ['s','e','ss', 'se', 'step','a','h', 'ex_n', 'ex3', 'ex1']:
    if k in options:
        options[k] = int( options[k] )

for k in [ 'zfilter', 'thin', 'mem']:
    if k in options:
        options[k] = float( options[k] )

for k in ['debug', 'verbose', 'cast', 'clean', 'w', 'all', 'heavy', 'jack',
          'protein']:
    if k in options:
        options[k] = 1

for k in [ 'cr', 'exrec', 'exlig', 'excom']:
    if k in options:
        options[k] = t.toIntList( options[k] )

if 'atoms' in options:
    options['atoms'] = t.toList( options['atoms'] )

if 'single' in options:
    options['hosts'] = mu.nonredundant( options.get('hosts', hosts.cpus_all) )

if 'h' in options:
    options['hosts'] = options.get('hosts', hosts.cpus_all)[ : options['h'] ]

if 'mem' in options:
    ram = hosts.ram_dic
    cpus= options.get('hosts', hosts.cpus_all)
    options['hosts'] = []
    for h in cpus:
        if ram.get(h, ram['default']) >= options['mem']:
            options['hosts'] += [ h ]
            ram[h] = ram.get(h, 0 ) - options['mem']


## ## initialize
master = AmberEntropyMaster( **options )

## ## ## ## run
master.start()
