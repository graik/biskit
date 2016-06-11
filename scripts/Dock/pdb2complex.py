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
##

import sys

import Biskit.tools as T
from Biskit import PDBModel
from Biskit.Dock import Complex as ProteinComplex

def _use( options ):
    print """
pdb2complex.py  - create a reference Complex (without waters)

Syntax:  pdb2complex.py  -c |complex pdb|
                         -r |chain index| -l |chain index|
                         -o |output name|
                     
Options:   -c     complex pdb file or pickled PDBModel object
           -r     receptor chain list (e.g. 0 1 ) 
           -l     ligand      ~       (e.g. 2 )
           -o     output file
           -lo,lr ligand, receptor model output file
           
Default options:
"""
    for key, value in options.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)

### MAIN ###
############

options = T.cmdDict( {'o':'ref.complex', 'lo':'lig.model', 'ro':'rec.model' } )

if len (sys.argv) < 3:
    _use( options )

## create a reference complex
print "Loading..."
ref_com =  PDBModel( options['c'] )
print "Removing water..."
ref_com.remove( lambda a: a['residue_name'] in ['TIP3','HOH','WAT'] )

## extract rec and lig chains
rec_chains = T.toIntList( options['r'] )
lig_chains = T.toIntList( options['l'] )

print "Extracting rec and lig..."
ref_rec = ref_com.takeChains( rec_chains )
ref_lig = ref_com.takeChains( lig_chains )

## create Protein complex
com = ProteinComplex( ref_rec, ref_lig )

print "Saving..."
ref_lig.saveAs( T.absfile( options['lo'] ) )

ref_rec.saveAs( T.absfile( options['ro'] ) )

T.dump( com, T.absfile( options['o']) )
