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

from Biskit.StructureMaster import StructMaster

import Biskit.tools as T
import Biskit.PVM.hosts as hosts
import sys

def _use( o ):

    print """
Syntax: pdbs2model.py -i |file1 file2 ..| [-h |host| -c |chunk| -a -w
                      -o |other_outFolder| -wat -s]

    pvm must be running on the local machine!

Result: pickled PDBModel object for each pdb file with same file name
        but ending in '.model'
        
Options:
    -h    number of hosts to be used
    -a    first add hosts to pvm
    -c    chunk size, number of pdb's passed to each node at once
    -w    display a xterm window for each node
    -o    destination folder (default: same where pdb file comes from)
    -wat  skip water residues (WAT TIP3 WWW H2O)
    -amber  rename CYX -> CYS, HID/HIE/HIP -> HIS, unwrap atom names
          (this creates models with the same atom/res names as pdbs created
           with ambpdb -p top.parm -aatm -bres < some.crd > some.pdb )
    -s    sort atoms alphabetically within residues

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)
     

if __name__ == '__main__':

    options = T.cmdDict( {'h':10, 'c':5 } )

    if len( sys.argv ) < 2:
        _use( options )

    hostNumber = int( options['h'] )
    chunk = int( options['c'] )
    windows = options.has_key('w')
    amber = options.has_key('amber')
    sort = options.has_key('s')

    outFolder = None
    if options.has_key('o'):
        outFolder = T.absfile( options['o'] )

    ## get all PDBs and models directly from a directory (avoids shell limits)
    if 'd' in options:
        d = T.absfile( options['d'] )
        files = os.listdir( d )
        files = [ x for x in l if x[-4:].upper() == '.PDB' or \
                  x[-6:] == '.MODEL' or x[-7:].upper() == '.PDB.GZ' ]
        files = [ os.path.join( d, f ) for f in files ]
    else:
    ## extract file names
        files = T.toList( options['i'] )

    ## add needed amount of hosts 
    add_hosts = options.has_key('a')

    skipWat = options.has_key('wat')

    fileDic = {}
    for f in files:
        fileDic[T.absfile( f )] = ''

    print "Initialize Job queue.."

    master = StructMaster(fileDic, chunk , hosts.cpus_all[ : hostNumber ] ,
                          outFolder, skipWat=skipWat, amber=amber,
                          sort=sort, show_output=windows, add_hosts=add_hosts )

    print "Start jobs .."
    master.start()

