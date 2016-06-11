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
import os.path as osp

def _use( o ):

    print """
Convert a pickled PDBModel into a PDB file.

Syntax: model2pdb.py -i |infile| -o |outfile| [ -wat -ter 0|1|2 -codeprefix ]

Options:
    -i      one or more pickled PDBModel(s)
    -o      output file name (default: infile.pdb ) (ignored if >1 input file)
    -wat    skip water residues (WAT TIP3 WWW H2O) and Cl-, Na+
    -ter 0  don't write any TER statements
    -ter 1  try restoring original TER statements
    -ter 2  put TER between all detected chains
    -codeprefix  add model's pdbCode entry as prefix to out file name

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)

##### MAIN ########

if __name__ == '__main__':

    options = cmdDict( {'ter':1} )

    if len( options ) < 2:
        _use( options )

    skip = None
    if options.has_key('wat'):
        skip = ['TIP3', 'WAT', 'HOH', 'H2O', 'Na+', 'Cl-']

    fin = [ absfile(f) for f in toList( options['i'] ) ]
    ter = int( options['ter'] )
    codeprefix = options.has_key( 'codeprefix' )

    for f in fin:

        fout= options.get('o', None) or stripSuffix( f ) + '.pdb'

        m = PDBModel( f )

        if codeprefix:
            code = m.pdbCode or ''
            fpath, fname = osp.split( fout )
            fout = osp.join( fpath, code + fname )

        m.writePdb( fout, ter=ter )

        errWriteln('dumped ' + fout )

