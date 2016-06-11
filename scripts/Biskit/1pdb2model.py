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

from Biskit.tools import *
from Biskit import PDBModel, PCRModel

def _use( o ):

    print """
Syntax: 1pdb2model.py -i |file1| [-o |outfile| -psf |psf_file| -wat -amber
                      -pdb |PDBCode| ]

Result: stand-alone pickled PDBModel or PCRModel
        
Options:
    -i      input PDB or pickled PDBModel
    -psf    psf file name -> will generate PCRModel instead
    -o      output file name (default: pdbfile.model)
    -wat    skip water residues (WAT TIP3 WWW H2O) and Cl-, Na+
    -amber  rename CYX -> CYS, HID/HIE/HIP -> HIS
    -pdb    pdb code to be stored in model.pdbCode

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


def go( fin, fout, fpsf=None, skipRes=None, amber=0, pdb=None ):

    if fpsf:
        m = PCRModel( absfile(fpsf), absfile(fin), skipRes=skipRes )
    else:
        m = PDBModel( absfile(fin), skipRes=skipRes )

    if amber:
        renameAmberRes( m )

    m.pdbCode = pdb or m.pdbCode

    m.saveAs( absfile(fout) )

##### MAIN ########

if __name__ == '__main__':

    options = cmdDict( {} )

    if len( options ) < 2:
        _use( options )

    skip = None
    if options.has_key('wat'):
        skip = ['TIP3', 'WAT', 'HOH', 'H2O', 'Na+', 'Cl-']

    if not options.has_key('o'):
        options['o'] = options['i'].split('.')[0]+'.model'

    amber = options.has_key('amber')
    
    go( options['i'], options['o'], options.get('psf',None), skip, amber,
        options.get('pdb',None) )
