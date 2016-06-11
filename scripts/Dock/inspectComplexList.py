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
from Numeric import *

def _use( o ):

    print """Check info dict of ComplexList for missing values.
Syntax: checkComplexes.py |complex_cont.cl|

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)

def __valuesOf( l, k ):
    x = l.valuesOf( k, None )
    x = [ sum(i) for i in x ]
    return x

def check( l ):
    print '%-20s  %6s' % ('infokey', 'invalid')
    print 30*'-'
    keyList = l[0].keys()
    keyList.sort()
    for k in keyList:
        x = __valuesOf( l, k )
        try:
            print '%-20s: %6i' % (k, sum( equal( x, None ) ) )
        except:
            print lastError()
            print '%-20s: %6s' % (k, '?')
            pass


if len( sys.argv ) < 2:
    options = cmdDict( {} )

    if len( options ) < 1:
        _use( options )
    

input = absfile( sys.argv[1] )

flushPrint( 'Loading...' )
cl = load( input )
flushPrint( '%i complexes read\n\n' % len( cl ) )


check( cl )

