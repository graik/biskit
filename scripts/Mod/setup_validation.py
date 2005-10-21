#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: Olivier PERIN
##
##
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
