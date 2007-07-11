#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## last $Author$
## $Date$
## $Revision$

import Biskit.Mod.modUtils as modUtils

from Biskit.Mod import *
from Biskit.Mod import Modeller as M
from Biskit.Mod import TemplateCleaner as TS
from Biskit.Mod import TemplateFilter as TF
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

from Biskit.Pymoler import Pymoler
import glob
from Biskit.Trajectory import Trajectory

import sys, os.path
import string

def _use( o ):
    print """
Build model using Modeller.

Syntax: model.py [ -o |outFolder| -log |logFile| -h |host_computer|
                   -v |verbosity|
                   -zfilter |x.x| -idfilter |percent| ]

Options:
    -o         output folder for results      (default: .)
    -log       log file                       (default: STDOUT)
    -h         host computer for calculation  (default: local computer)
               -> must be accessible w/o password via ssh, check!
    -s         show structures on Pymol superimposed on average
    -dry       do not run Modeller; just set up and try post-processing
               existing results
    -zfilter   ignore templates that are more than z standard deviations below
               the average sequence identity to the target (0..switch off)
	       (default: see Biskit.Mod.TemplateFilter.py)
    -idfilter  ignore templates that have less than idfilter fraction sequence
               identity to target sequence (the best template is always kept)
	       (0..switch filtering off)
    -v         verbosity level (0 - 3)
    -? or help .. this help screen

input: ./templates/modeller/*.pdb
       ./t_coffee/final.pir_aln

output: ./modeller/modeller.log
                  /*.B9999000??   <- models
		  /model_??.pdb   <- processed PDBs ordered by modeller score

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
    sys.exit(0)


def defaultOptions():
    return {'o':'.',
            'log': None,
            'h':None,
	    'v':1,
	    'zfilter':TF.Z_CUTOFF,
	    'idfilter':TF.ID_CUTOFF,
            }

def convertOptions( o ):
    o['verbose'] = int( o.get('v',1) )
    del o['v']
    
    o['outFolder'] = tools.absfile( o.get('o','.') )
    del o['o']
    
    o['zfilter'] = float( o.get('zfilter',0) )
    o['idfilter']= float( o.get('idfilter',0))

    o['log'] = o.get('log', None)
    if o['log']:
	o['log'] = LogFile( o['o'] + '/' + options['log'], 'a' )

    o['host'] = o.get('h', None)
    del o['h']


    return o

def __printMatrix( matrix ):
    """
    Print the part right of the diagonal in a matrix
    """
    nr = len( matrix )
    for i in range(nr): print '%5i'%(i+1),
    for i in range(nr):
        print '\n%2i'%(i+1),
        for k in range(i):
            print ' '*5,
        for j in range(i, nr):
            print '%5.2f'%matrix[i,j],
            
### MAIN ###

options   = tools.cmdDict( defaultOptions() )

if '?' in options or 'help' in options:
    _use( defaultOptions() )

options = convertOptions( options )

## Test, comment in for testing
# outFolder = tools.testRoot() + '/Mod/project'

###################
## Modeller
##
## Build model using Modeller.

## input: templates/modeller/*.pdb
##        t_coffee/final.pir_aln
##
## output: modeller/modeller.log
##                 /*.B9999000??   <- models
try:

    if options['verbose'] > 0:
        print "\n"+\
              "Type model.py -? or model.py -help for a full list of options!"

    m8 = M( **options )
   
    if not 'dry' in options:
	r = m8.run( options['host'] )  ## comment out for testing

except:
    EHandler.error( 'Error while modelling.')


#####################
## Show output

## show result in PyMol
if options.has_key('s'):
    names = []

    ## fit backbone of all models to average until convergence
    models =  glob.glob( '%s/target.B*'%(outFolder + m8.F_RESULT_FOLDER) )
    traj = Trajectory( models )
    traj.blockFit2ref( mask=traj[0].maskBB() )

    ## calculate and print rmsd matrix
    rmsHeavy = traj.pairwiseRmsd()
    print '\nHEAVY ATOM RMSD BETWEEN MODELS::'
    __printMatrix( rmsHeavy )

    ## same thing for backbone atoms
    BBMask = traj[0].maskBB()
    traj.blockFit2ref( mask = BBMask )
    rmsBB = traj.pairwiseRmsd( aMask = BBMask  )
    print '\nBACKBONE RMSD BETWEEN MODELS:'
    __printMatrix( rmsBB )

    ## get all templates
    templates =  glob.glob( '%s*'%( outFolder + TS. F_MODELLER ))
    
    p=Pymoler()
    ## add models
    for i in range( len(traj) ):
	print 'Pymol.add model', i
        p.addPdb( traj[i], 'model_%i'%i ) 
        names  += [ 'model_%i'%i ]

    ## load all templates and align to first model
    for t in templates:
        name = string.split( string.split(t, '/')[-1], '.')[0]
        names  += [ name ]
        p.addPdb( t, name )
        p.add( 'align %s, model_0'%name )

    ## color things in different colors
    colors = p.addColors( len(names),
                          firstColor=[1.0, 0.0, 0.0],
                          lastColor=[0.0, 0.0, 1.0] )

    for i in range( len(names) ):
        p.add('color %s, %s'%(colors[i], names[i]) )
                
    p.add('hide all')
    p.add('show line')

    ## don't show sidechains and carboxyl oxygen - too messy
    p.add('select bb, not name c+ca+n')
    p.add('hide line, bb')
    
    p.add('select none')
    p.add('zoom all')

    print 'writing pymol script and models...'
    
    p.show()
