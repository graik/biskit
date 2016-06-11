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

## Contributions: Olivier PERIN, Raik Gruenberg


from Biskit.Mod.Analyse import Analyse as A
from Biskit.Mod.ValidationSetup import ValidationSetup as VS
from Biskit.Mod.Benchmark import Benchmark as B
import Biskit.tools as T
import sys, os
import Biskit.Pymoler as Pymoler

def _use( o ):

    print """
Syntax: analyse.py -d |main project folder| [-s |1||0] ]
                       
Result: Performing model analysis for each main project folder given.
        Outputs a folder 'analyse' containing:
        
        * analyse/global_results.out
          various data about the model, see file header.

        * analyse/local_results.out:
          residue rmsd profile to taget and mean rmsd to tagret

        * modeller/final.pdb:
        the 'best' model with the mean residue rmsd in the B-factor column


Options:
        -d          [str], list of project directory
        -s          show the structure final.pdb im PyMol
"""
    
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = T.cmdDict()
    f = os.getcwd()
    
    if '?' in options or 'help' in options:
        _use( options )

    if not os.path.exists(f + VS.F_RESULT_FOLDER) and not options.has_key('d'):
        print 'Current directory is not a valid modeling project folder.' 
        _use( options )
        
    ## Try to add project folders
    ## look for default cross-validation projects
    d = []
    if os.path.exists( f + VS.F_RESULT_FOLDER ):
        d = [ f ]

    if options.has_key('d'):
        folders = T.toList(options['d'])
    else: 
        folders = d

                       
    T.flushPrint("Starting job...\n")

    for f in folders:
            a = A(outFolder=f)
            a.go()

    T.flushPrint("Done.\n")
    
    ## show result in PyMol
    if options.has_key('s'):
        p=Pymoler()
        p.addPdb( folders[0] + a.F_FINAL_PDB )
        p.add('color_b')
        p.add('select na, b<0')
        p.add('color grey, na')
        p.add('hide all')
        p.add('show cartoon')
        p.add('show stick')
        p.add('hide stick, name o+c+n')
        p.add('select ca, /////ca')
        p.add('label ca,"%s-%s"%(resn, resi)')
        p.add('select none')
        p.show()


        
