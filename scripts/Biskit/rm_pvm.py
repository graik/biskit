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

## delete all /tmp/pvm* files on all hosts

import os
import sys, string
from Biskit.PVM.hosts import nodes_all

def cleanPvm( host ):
    cmd = "ssh %s 'rm /tmp/pvm*'"
##     cmd = "ssh %s 'rm -f /work/@*'"
    os.system( cmd % host )

def df( host ):
    cmd = "ssh %s 'cd /tmp; df -h .'"
    os.system( cmd % host )

s = "removing local /tmp/pvm* on ALL hosts! Continue? (Y/N)+ENTER: "
if string.upper( raw_input(s) ) != 'Y':
    print "Cancelled"
    for h in nodes_all:
        try:
            print h
            df( h )
        except:
            pass
    sys.exit(0)

for h in nodes_all:

    try:
        print h
        cleanPvm( h )
    except Exception, why:
        pass
##         print "Error cleaning %s: ", why
