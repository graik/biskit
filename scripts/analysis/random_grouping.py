#!/usr/bin/env python
## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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
## select set of non-redundant, non-native random complexes
##

import time
import random
import Biskit.oldnumeric as N0
import os.path, os

from Biskit.tools import *
from Biskit.Dock import ComplexRandomizer, ComplexList
from Biskit.Dock import ContactMaster, ComplexGroups

from Biskit import Trajectory
from Biskit.PVM.hosts import nice_dic, cpus_own

default_options = {'o':'random_complexes_',
                   'ref':'../../com_wet/ref.complex',
                   'h':len(cpus_own),
                   'co':'./coms'}

def _use():
    print """
random_grouping.py -cl |complex_list|

Options:
      ref .. pickled native complex
      h   .. number of hosts
      co  .. folder name for result complexes
      o   .. base name for other result files
      a   .. add hosts to PVM before starting

Default options:
"""
    for key, value in default_options.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)


def contact( cl, refComplex, fout, add_hosts=1 ):
    """
    Contact the list of random complexes
    -> contacted ComplexList
    """
    flushPrint('Contacting...\n')

    master = ContactMaster(cl, 1,
                           cpus_own[:int(options['h'])],
                           refComplex = refComplex,
                           force = ['fnarc_10','c_ratom_10'],
                           niceness = nice_dic,
                           outFile = fout,
                           add_hosts = add_hosts,
                           show_output = 0 )

    return master.calculateResult()


def filter_zero_contacts( cl ):
    """
    """
    s = [ len( c['c_ratom_10']['nonzero'] ) for c in cl ]

    return cl.take( N0.nonzero( s ) )

def group( cl, pw, pwmin ):
    """
    Cluster the list.
    pw    - float, min overlap with center complex
    pwmin - float, min overlap between all cluster members
    -> ComplexGroups
    """
    flushPrint('grouping...')
    cg = ComplexGroups( cl )

    cg.group( pw, pwmin )

    cg.assignClusters()

    cg.reportOverlaps()

    flushPrint('done grouping\n')
    return cg


def selectClusters( cg, n ):
    """
    select n out of all cluster centers in cg
    cg - ComplexGroups
    n  - int, number of clusters to select
    -> ComplexList
    """
    ## centers only
    nr_cl = cg.centerComplexes()

    ## randomize order
    for c in nr_cl:
        c['r'] = random.random()
    nr_cl = nr_cl.sortBy('r')

    if n > len( nr_cl ):
        n = len( nr_cl )

    return nr_cl[:n]


def saveComplexes( cl, folder ):
    """
    """
    if not os.path.exists( folder ):
        os.mkdir( folder )
        
    for i in range(len(cl)):
        dump( cl[i], '%s/random_%02i.complex' % (folder,i) )

#######
# MAIN

if len(sys.argv) < 2:
    _use()
    
options = cmdDict( default_options )

## parameters
try:
    fout = absfile( options['o'] )
    cl   = load( options['cl'] )
    n_cl = len( cl )
    ref  = load( options['ref'] )
    com_out = absfile( options['co'] )
except:
    print "Missing or wrong option:"
    print lastError()
    _use()
    
add_hosts = options.has_key('a')

## contacting and very rough grouping
cl = contact( cl, ref, fout + 'cont.cl', add_hosts=add_hosts )

cl = filter_zero_contacts( cl )  ## happened in one case

## cl = cl.filter( 'fnac_10', (0., 0.) )
cl = cl.filter( 'fnarc_10', (0., 0.) )
flushPrint('%i orientations removed because of fnarc > 0.0\n' % (n_cl-len(cl)))

cg = group( cl, 0.001, 0.0001 )

dump( cg, fout+'grouped.cg' )

## final selection
nr_cl = selectClusters( cg, 10 )
flushPrint(
    'dumping %i non-redundant, non-overlapping, non-native complexes\n'\
        % len( nr_cl ) )
dump( nr_cl, fout+'nr.cl')

flushPrint('Dumping random complexes.')
saveComplexes( nr_cl, com_out )

## for visualisation
flushPrint('\nPreparing trajectory...')

t = Trajectory( [ c.model() for c in nr_cl ] )
t.ref.addChainId()
t.ref.writePdb( fout+'nr_traj_ref.pdb', ter=2 )
t.writeCrd( fout+'nr_traj.crd')

