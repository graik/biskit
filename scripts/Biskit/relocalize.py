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
## Re-create LocalPath in PDBModel or Complex from the original file name
## todo: Add options to replace substitutions or fragments of the original
##       path

import Biskit.tools as tools
from Biskit import PDBModel, PCRModel, LocalPath
from Biskit.Dock import Complex

import os, sys

def localizeComplex( com, replace ):
    localizeModel( com.rec_model, replace )
    localizeModel( com.lig_model, replace )
    return 1

def localizeModel( m, replace, ftarget=None ):
    """
    m - PDBModel
    replace - [ [str, str], [..]]
    ftarget - str, target file name
    """
    result = 0
    print '\t', m.source.formatted(),
    f = m.source.original()
    s = m.source

    for repl in replace:
        f = f.replace( repl[0], repl[1] )
    
    m.source = LocalPath( f )
    print " -> ", m.source.formatted()

    ## don't slim if model is it's own source
    if ftarget and tools.absfile( ftarget ) == m.source.local():
        m.forcePickle = 1

    return 1


##########
## MAIN ##

options = tools.cmdDict()

if not 'i' in options:
    print """Re-create LocalPath in PDBModel or Complex from the original
file name (and the current environment variables).
relocalize.py -i |file1 file2 ..|
"""
    sys.exit(0)


files = tools.toList( options['i'] )

repl = []
if 'r' in options:
    relacements = tools.toList( options['r'] )
    repl = [ r.split( ':' ) for r in relacements ]

for f in files:
    try:
        print "Re-localizing ", f
        f = tools.absfile( f )
        o = tools.load( f )
        result = 0

        if o.__class__ in [ PCRModel, PDBModel ]:
            result = localizeModel( o, repl, f )

        else:
            if o.__class__ in [ Complex ]:
                result = localizeComplex( o, repl )

            else:
                print "unknown class", o.__class__

        if result:
            f_bak = f + '__old'
            os.rename( f, f_bak )
            tools.dump( o, f )
            print "..done"
        else:
            print "..skipped"


    except Exception, why:
        print "ERROR converting %s: %s" % (f, str(why))
            
