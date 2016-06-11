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
## create list of random complex orientations
##

from Biskit.tools import *
from Biskit.Dock import ComplexRandomizer, ComplexList
from Biskit import Trajectory

default_options = {'o':'random_complexes.cl', 'n':'25',
                   'ro':'rec_centered.model', 'lo':'lig_centered.model'}

def _use():
    print """
random_complexes.py -r |rec_model| -l |lig_model| [ -o |out_file| -n |number|
                    -traj |traj_out_name| -ro |rec_out| -lo |lig_out|
                    -debug copy_inp_file ]

Remark:
to create valid PCRModels for rec and lig use (in rec_wet/, lig_wet/)
        1pdb2model.py -i ????.pdb -psf ????.psf -o xplor.model
The waters in the PSF are deleted. They can be in the model but don't have to. 

Options:
     r     pickled PCRModel, receptor (psf file has to be valid)
     l     pickled PCRModel, ligand   (psf file has to be valid)
     ro    file name for rec copy (centered and no waters)
     lo    file name for lig copy (centered and no waters)
     o     file name for result ComplexList
     n     number of random complexes to generate
     traj  file name for optional amber crd and pdb (for visualisation)
     debug keep Xplor input file and don't delete temporary files

Default options:
"""
    for key, value in default_options.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)



#### MAIN #####

if len(sys.argv) < 2:
    _use()
    
options = cmdDict( default_options )

try:
    rec = load( options['r'] )
    lig = load( options['l'] )
    ro  = absfile( options['ro'] )
    lo  = absfile( options['lo'] )
    out = absfile( options['o'] )
    n   = int( options['n'] )
    traj= absfile( options.get('traj', None) )
    finp= absfile( options.get('debug', None))
    debug = ( finp is not None )
except:
    print "Missing or wrong option:"
    print lastError()
    _use()
    
cr = ComplexRandomizer( rec, lig, ro, lo, debug=debug )
result = ComplexList()

flushPrint('generating...')
for i in range( n ):
    result.append( cr.random_complex( inp_mirror=finp ) )

    flushPrint('#')

flushPrint('\nDumping result...')
dump( result, out )


if traj:
    flushPrint('\nWriting trajectory...')
    t = Trajectory( [ c.model() for c in result ] )
    t.ref.writePdb( traj + '.pdb' )
    t.writeCrd( traj + '.crd' )

flushPrint('done')
