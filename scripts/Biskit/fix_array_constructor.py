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
    o = t.load( fin )
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
    o = t.load( fin )
    print "Done"
except ImportError:
    print "could not fix %s" % fin
    os.rename( fout, fin )
except:
    print "There is still some other problem: "
    print t.lastError()
    print t.lastErrorTrace()
    
