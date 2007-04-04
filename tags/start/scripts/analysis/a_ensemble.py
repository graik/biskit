#!/usr/bin/env python
## -*- coding: iso-8859-1 -*-
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Date$
## last $Author$

import sys
import os.path as osp

import biggles
import Numeric as N

from Biskit.tools import *
from Biskit.EnsembleTraj import EnsembleTraj, traj2ensemble

def syntax():

    if len( sys.argv ) < 2:
        
        print """
    Analyze ensemble MD.
    Syntax:  a_ensemble.py -i traj.dat [ -z |outlier-z-score| -step |offset|
                           -o |out.eps| -title |plot_title| ]
    Default options:
    """
        o = options()
        for key, value in o.items():
            print "\t-",key, "\t",value

        sys.exit(0)
        

def options():
    default = { 'z':1.0, 'step':1 }
    return cmdDict( default )

def calcRmsd( traj ):
    mCA = traj.ref.maskCA()
    t.fit( ref=t.ref, prof='rmsCA_ref', comment='CA to ref')
    t.fitMembers( mask=mCA, prof='rmsCA_av', comment='CA to member avg')
    t.fitMembers( refIndex=-1, mask=mCA, prof='rmsCA_last',
                  comment='all CA to last member frame')

def plot( traj, title ):
    p = t.plotMemberProfiles('rmsCA_ref', 'rmsCA_av', 'rmsCA_last')
    p.title = title
    p.xlabel= 'frame #'
    p.ylabel= 'RMSD [�]'
    return p

def markOutliers( traj, z, page ):

    outliers = N.nonzero( traj.outliers( z=z, mask=traj.ref.maskCA() ) )

    for o in outliers:
        t = traj.takeMember( o )

        ## cross out outliers in plot
        prof = N.array( t.profiles['rmsCA_ref'] ).tolist()
        prof.extend( t.profiles['rmsCA_last'] )
        maxV = max( prof )
        
        line = biggles.Line( (0,0), (len(t),maxV) )

        page[ o / 2, o % 2 ].add( line )


########
## MAIN
########

syntax()

## get and clean up options
o = options()
o['z'] = float( o['z'] )
o['step'] = int( o['step'] )
o['i'] = absfile( o['i'] )
o['o'] = o.get('o',
               '%s/%s_rms.eps' % (osp.dirname(o['i']), stripFilename(o['i'])))

print "Loading...",
t = Load( o['i'] )
print "done."

if not isinstance( t, EnsembleTraj ):
    t = traj2ensemble( t )

if o['step'] != 1:
    t = t.thin( o['step'] ) 

print "Fitting...",
calcRmsd( t )
print "done."

p = plot( t, o.get( 'title', stripFilename(o['i']) ) )

print "Getting outliers ...",
markOutliers( t, o['z'], p )
print "done."

p.show()

print "Saving plot...",
p.write_eps( o['o'], width="18cm", height="29cm" )
print "done."