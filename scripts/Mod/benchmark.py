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

## Contributions: Olivier PERIN

from Biskit.Mod.Benchmark import Benchmark as B
from Biskit.Mod.ValidationSetup import ValidationSetup as VS
import Biskit.tools as T
import sys, os
import glob
import os.path

def _use( o ):

    print """
Syntax: benchmark.py -d |list of folders| 
                     [ -modlist |model_list| -ref |reference|]
                                                            
Result: Performing various benchmark tasks for each folder given.
        A folder validation/benchmark containing:

        * validation/????/benchmark/Fitted_??.pdb:
        Benchmark model iteratively superimposed on its known structure.

        * validation/????/benchmark/rmsd_aa.out:
        All-atom rmsd of the benchamark modela. (1) without iterative fitting,
        (2) with iterative fitting and (3) the percentage of atoms that has
        been removed during the iterative fitting.

        * validation/????/benchmark/rmsd_ca.out:
        same as above, but only for C-alpha atoms

        * validation/????/benchmark/rmsd_res_??:
        gives the C-alpha rmsd for each residue.

        * validation/????/benchmark/PDBModels.list:
        pickled PYTHON list of PDBModels. Each model contains
        the benchmark information in the atom and residue profiles:
        'rmsd_aa', 'rmsd_ca', 'rmsd_res'. See PDBModel.profile()!

        
Options:
    -d          [str], list of project validation directories
    -modlist    str, the path to the 'PDBModels.list' from the
                  project directory
    -ref        str, the path to the 'reference.pdb' from
                  the project directory (known structure)
"""
    
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = T.cmdDict()
    f = os.getcwd()
    
    if '?' in options or 'help' in options or 'h' in options:
        _use( options )

    if not os.path.exists( f + B.F_INPUT_FOLDER ) and not options.has_key('d'):
        print 'Current directory is not a valid modeling project folder.' 
        _use( options )
        
    ## Try to add project folders
    ## look for default cross-validation projects
    d = []
    
    if os.path.exists( f + VS.F_RESULT_FOLDER ):
        d = glob.glob( f + VS.F_RESULT_FOLDER + '/*' )

    if options.has_key('d'):
        folders = T.toList(options['d'])
    else: 
        folders = d
    
    reference = options.get('ref', None)
        
    model_list = options.get('modlist', None)

                       
    T.flushPrint("Starting job...\n")

    for f in folders:
        T.flushPrint("\tWorking on %s\n"%os.path.split(f)[1])
        b = B( outFolder=f )
        b.go(model_list = model_list, reference = reference)

    T.flushPrint( "Done.\n")
