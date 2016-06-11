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

## Contributions: Olivier PERIN

from Biskit.Mod.ValidationSetup import ValidationSetup as VS
import Biskit.tools as T
import sys, os, os.path


def _use( o ):

    print """
Setup the cross-validation folder for one or several projects
        
Syntax: setup_validation.py [ -o |project folder(s)| ]
                       
Options:
    -o          .. one or several project folders (default: current)
    -? or -help .. this help screen

Default options:
"""
    
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = T.cmdDict({'o':[ os.getcwd() ]})

    if '?' in options or 'help' in options:
        _use( options )
                       
    folders = T.toList( options['o'] )

    if not os.path.exists( folders[0] +'/templates'):
        print 'Current directory is not a valid modeling folder.' 
        _use( options )

    T.flushPrint( "Creating folders and links...\n" )
  
    for f in folders:
        sv = VS(outFolder=f)
        sv.go(f)

    T.flushPrint( "done\n" )
