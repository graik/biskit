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

## Create an incomplete Amber topology from a PDB.
## The parm file will have the same atom content as the PDB -
## i.e. it will have missing atoms and cannot be used for simulation.
## 

import sys

from Biskit import AmberParmBuilder, PDBModel, LogFile
import Biskit.tools as t

def _use( options ):
    print """
Create amber topology and coordinate files from a PDB.

Syntax: am_pdb2parm.py -i |PDBfile| [-o |out.parm| ...any of options below ]
    OR: am_pdb2parm.py -i |PDBfile| -mirror [-o |out.parm| ]
    
Result: out.parm, out.crd, out.pdb, (and leap.log in current folder)

Special option -mirror: create a parm for exact atom content of input PDB
                      (no S-S bonds formed, atoms missing from residues..)
                      This parm can be used for ptraj but not for simulations!
Options:
        ocrd      - str, target file for crd (coordinates) [|f_out_base|.crd]
        opdb      - str, target file for pdb               [|f_out_base|.pdb]
        hetatm    - keep hetero atoms                                 [don't]
        cap       - put ACE and NME capping residue on chain breaks   [don't]
        capN      - int int, indices of chains that should get ACE cap     []
        capC      - int int, indices of chains that should get NME cap     []
        box       - float, minimal distance of solute from box edge    [10.0]
        nocenter  - do *not* re-center the input coordinates         [center] 
        fmod      - str str, list of files with amber parameter modifications
                    (to be loaded into leap with loadAmberParams)          []
        fprep     - str str, list of files with amber residue definitions
                    (to be loaded into leap with loadAmberPrep)            []

        leap_template - str, path to template file for leap input [use default]
        leaprc        - str, path to parameter file for leap [use default]
        leap_out      - str, target file for leap.log [default: discard]
        leap_in       - str, target file for leap.in script [default: discard]
        leap_pdb      - str, target file for cleaned input pdb      [discard]
        leap_bin      - str, path to tleap binary [use default]
        norun         - do not run leap, only prepare files
        debug         - keep all temporary files

        more -key value pairs for place holders in  leap input template
       
Comments:
- The protonation state of histidines is decided from the H-atoms that
  are found (HD, HE, or both). After that all H are removed to be added again
  by leap.
- Cleaning tries to convert non-standard residues to the closest standard one.
- Cleaning removes non-standard atoms (and atoms following them) from standard
  residues.
- Cleaning keeps the largest / first of multiple occupancies
- Ends of chains are assumed if the residue numbering jumps backward, if there
  is a TER record or chain ID or segid change, or if there is a chain break.
- A chain break is assumed if there is an untypical gap in the chain of back-
  bone atoms (see PDBModel.chainBreaks() ).
- The index of the first chain is 0.
- Original waters are deleted. 
- As usual, options can also be put into a file and loaded with the -x option

Default options:
"""
    for key, value in options.items():
        print "\t-",key, "\t",value
        
    #sys.exit(0)

options = t.cmdDict( {'o':'out.parm'} )

try:

    f_out = options['o']

    if 'ocrd' in options: 
        options['f_out_crd'] = options['ocrd']
    if 'opdb' in options:
        options['f_out_pdb'] = options['opdb']
    if 'box' in options:
        options['box'] = float( options['box'] )

    options['cap'] = 'cap' in options
    options['capN']= t.toIntList( options.get('capN',[]))
    options['capC']= t.toIntList( options.get('capC',[]))
    options['hetatm'] = 'hetatm' in options
    options['norun'] = 'norun' in options
    options['debug'] = 'debug' in options
    options['center'] = not 'nocenter' in options

    if 'log' in options:
        options['log'] = LogFile( options['log'] )

    if 'norun' in options:
        fbase = t.stripSuffix( t.absfile( options['i'] ) )
        options['leap_in'] = options.get('leap_in', fbase+'_leap.in')
        options['leap_pdb']= options.get('leap_pdb',fbase+'_forleap.pdb')

    a = AmberParmBuilder( options['i'], **options )
    del options['debug']

    if not 'mirror' in options:
        a.parmSolvated( f_out, **options )

    else:
        a.parmMirror( f_out, **options )

except KeyError, why:
    _use( options )

except Exception, why:
    print "There was an error..."
    print t.lastError()
    print t.lastErrorTrace()
