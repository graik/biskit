#!/usr/bin/env python
## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## -*- coding: iso-8859-1 -*-
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
##

import sys
import os.path as osp

import biggles
import numpy as N

import biskit.tools as T
from biskit.md.ensembleTraj import EnsembleTraj, traj2ensemble

def syntax():

    if len( sys.argv ) < 2:
        
        print("""
    Analyze normal MD.
    Syntax:  a_trajRmsd.py -i traj.dat [ -step |offset|
                           -o |out.eps| -title |plot_title| -show ]
    Default options:
    """)
        o = options()
        for key, value in list(o.items()):
            print("\t-",key, "\t",value)

        sys.exit(0)
        

def options():
    default = { 'step':1 }
    return T.cmdDict( default )

def calcRmsd( traj ):
    mCA = traj.ref.maskCA()
    traj.fit( ref=traj.ref, prof='rmsCA_ref', comment='rmsCA to ref')
    traj.fit( mask=mCA, prof='rmsCA_av', comment='rmsCA to avg')
    traj.fit( ref=traj[-1], mask=mCA, prof='rmsCA_last',
                  comment='all CA to last member frame')

    traj.fit( mask=traj.ref.maskHeavy(), prof='rmsHeavy_av', comment='rms heavy atoms to avg')
    traj.fit( mask=traj.ref.maskHeavy(), ref=traj[-1], prof='rmsHeavy_last',
                  comment='all heavy atoms to last member frame')

    

def plot( traj, title ):
    p = t.profiles.plotArray('rmsCA_ref', 'rmsCA_av',('rmsHeavy_last', 'rmsCA_last' ))
    p.title = title
    p.xlabel= 'frame #'
    p.ylabel= 'RMSD [Å]'
    return p


########
## MAIN
########

syntax()

## get and clean up options
o = options()
o['step'] = int( o['step'] )
o['i'] = T.absfile( o['i'] )
o['o'] = o.get('o',
              '%s/%s_rms.eps' % (osp.dirname(o['i']), T.stripFilename(o['i'])))
o['show'] = 'show' in o

T.flushPrint( "Loading..." )
t = T.load( o['i'] )
T.flushPrint( "done loading trajectory with %i frames." % len(t) )

if o['step'] != 1:
    t = t.thin( o['step'] ) 

T.flushPrint( "Fitting ...")
calcRmsd( t )
T.flushPrint( "done." )

p = plot( t, o.get( 'title', T.stripFilename(o['i']) ) )

if o['show']:
    p.show()

T.flushPrint( "Saving plot to %s" % o['o'] )
p.write_eps( o['o'], width="18cm", height="29cm" )
T.flushPrint( "done." )
