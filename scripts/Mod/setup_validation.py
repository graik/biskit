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

## Contributions: Olivier PERIN
## last $Author$
## last $Date$
## $Revision$

from Biskit.Mod.ValidationSetup import ValidationSetup as VS
import Biskit.tools as T
import sys


def _use( o ):

    print """
Syntax: setup_validation.py -d |list of folders| 
                       
Result: Setup the validation folder for each project directory given
        
Options:
        -d    [str], list of project directories
"""
    
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    options = T.cmdDict()

    if len( sys.argv ) < 2:
        _use( options )

                       
    folders = T.toList(options['d'])

    print "Initialize Job queue..\n"
  
    for f in folders:
        sv = VS(outFolder=f)
        sv.go(f)

    print "Now it's done..\n"
