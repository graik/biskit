#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
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


def check( l ):
    print '%-20s  %6s' % ('infokey', 'invalid')
    print 30*'-'
    for k in l[0].keys():
        x = l.valuesOf( k, None )
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
cl = Load( input )
flushPrint( '%i complexes read\n\n' % len( cl ) )


check( cl )

