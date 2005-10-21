#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## Fix the import of Numeric.array_constructor in pickled Biskit modules
## Many Biskit pickles have e.g. array_constructor dumped as
## PDBModel.array_constructor

import shutil
import Biskit.tools as t
import sys, os

fin = t.absfile( sys.argv[1] )
fout = fin + '_backup'

exit = 0

try:
    o = t.Load( fin )
    print "%s looks alright.\nnothing to be done." % fin
    exit = 1
except ImportError:
    pass
except AttributeError:
    pass
except:
    print "Something else is wrong with %s:" % fin
    print t.lastErrorTrace()
    exit = 1

if exit:
    sys.exit()

print fin, '->', fout

s = open( fin ).read()

s2 = s.replace( 'cBiskit.PDBDope\narray_constructor',
                'cNumeric\narray_constructor' )

s2 = s2.replace( 'cBiskit.PDBModel\narray_constructor',
                 'cNumeric\narray_constructor' )

s2 = s2.replace( 'cPDBModel\narray_constructor',
                 'cNumeric\narray_constructor' )

s2 = s2.replace( 'cMathUtils\narray_constructor',
                 'cNumeric\narray_constructor' )

s2 = s2.replace( 'cfit\narray_constructor',
                 'cNumeric\narray_constructor' )

s2 = s2.replace( 'cTrajectory\narray_constructor',
                 'cNumeric\narray_constructor' )

s2 = s2.replace( 'cBiskit.Dock.ContactMaster\narray_constructor',
                 'cNumeric\narray_constructor' )

shutil.copy( fin, fout )

f = open( fin, 'w' )
f.write( s2 )
f.close()

try:
    print "checking...",
    o = t.Load( fin )
    print "Done"
except ImportError:
    print "could not fix %s" % fin
    os.rename( fout, fin )
except:
    print "There is still some other problem: "
    print t.lastError()
    print t.lastErrorTrace()
    
