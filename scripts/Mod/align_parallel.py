#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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

## Contributions: Olivier PERIN, Raik Gruenberg
## last $Author$
## last $Date$
## $Revision$

from Biskit.Mod.AlignerMaster import AlignerMaster
import Biskit.tools as T
import Biskit.hosts as hosts
import sys

def _use( o ):

    print """
Syntax: align_parallel.py -d |list of folders| -h |host|
                       [-pdb |pdbFolder| -ft |fastaTemplates|
                       -fs |fastaSequences| -fta |fastaTarget|
                       -fe |ferror|]

    pvm must be running on the local machine!

Result: Parallel alignment for each project directory given
        
Options:
    -d    [str], list of project directory (full path)
    -h    int, number of hosts to be used
    -a    first add hosts to pvm
    -pdb  str, pdbFolder for the pdb *.alpha
    -ft   str, path to 'templates.fasta'
    -fs   str, path to 'nr.fasta'
    -fta  str, path to 'target.fasta'
    -fe   str, path to the error file for the AlignerMaster
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = T.cmdDict()

    if len( sys.argv ) < 4:
        _use( options )

                       
    folders = T.toList(options['d'])
    hostNumber = int(options['h'])

   
    pdbFolder = options.get('pdb', None)

    fastaTemplates = options.get('ft', None)

    fastaSequences = options.get('fs', None)

    fastaTarget = options.get('fta', None)

    ferror = options.get('fe',None)

    show_output = 'w' in options


    print "Initialize Job queue.."

    master = AlignerMaster(hosts=hosts.cpus_all[ : hostNumber ],
                           folders=folders,
                           pdbFolder = pdbFolder,
                           fastaTemplates=fastaTemplates,
                           fastaSequences=fastaSequences,
                           fastaTarget=fastaTarget,
                           ferror=ferror,
                           show_output=show_output)
    
    print "Start jobs .."
    master.calculateResult()
