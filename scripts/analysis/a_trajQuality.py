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

from Biskit.QualMaster import QualMaster


def _use( options ):
    print """
Syntax: a_trajQuality -i |traj_1 traj_2 .. traj_n| [-a -h |n_hosts| -w]
    pvm must be running on the local machine!

Result: eps with quality plots in folder of traj files
        
Options:
    -h    number of hosts to be used
    -a    first add hosts to pvm
    -w    display a xterm window for each node

Default options:
"""
    for key, value in options.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = cmdDict( {'h':20} )

    if len( sys.argv ) < 2:
        _use( options )

    ## extract file names
    files = toList( options['i'] )

    master = QualMaster( files, int(options['h']),
                         show_output=options.has_key('w'),
                         add_hosts=options.has_key('a') )

    master.start()
