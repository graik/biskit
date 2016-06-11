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
from Biskit.PVM.TrackingJobMaster import restart

def use():
    if len( sys.argv ) < 2:
        print \
"""
Restart a distributed calculation.
Syntax:  restartPVM.py -i |rst_file| [-a]
Options:
         i  .. restart file containing result of TrackingJobMaster.getRst()
         a  .. add hosts to PVM
"""
        sys.exit(0)

## MAIN ##

use()

cmd = T.cmdDict()

T.flushPrint('Loading restart data...')
rst = T.load( cmd['i'] )

hosts = []
if 'a' in cmd:
    hosts = [ h['host'] for h in rst['hosts'] ]

master = restart( rst, hosts=hosts )

T.flushPrint('Master initialized for restart.')

master.start()





