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
import sys, pickle, cPickle

import Biskit.tools as T
from Biskit import PDBModel

def __use():
    if len( sys.argv ) < 2:
        print """
Update old python pickles. Unpickle some python objects and pickle them back
to the same filename.  Disable strict class checking if necessary -- this
allows to load classes that have meanwhile changed their base
class.

Usage:
redump.py  file1 file2 ..
"""
        sys.exit(0)


def resolveClass(module, name):
    exec( 'import %s' % module )
    return eval( '%s.%s' % ( module, name ) )

class PickleUpgrader(pickle.Unpickler):
    """
    Adapted from
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/286203
    """

    def find_class(self, module, cname):
        # Pickle tries to load a couple things like copy_reg and
        # __builtin__.object even though a pickle file doesn't
        # explicitly reference them (afaict): allow them to be loaded
        # normally.
        if module in ('copy_reg', '__builtin__'):
            thing = pickle.Unpickler.find_class(self, module, cname)
            return thing
        return resolveClass(module, cname)


def sloppyload( f ):
    """
    f - str, file name
    -> any, unpickled object
    """
    try:
        T.flushPrint( "Loading " + str(f) + '\n' )

        return T.load( T.absfile( f ) )

    except cPickle.UnpicklingError:

        print "Trying to load %s in sloppy mode..." % f
        return PickleUpgrader(open(T.absfile(f))).load()
    

#########
## Main
#########

__use()

fs = sys.argv[1:]

for f in fs:

    try:
        o = sloppyload( f )

        ## don't slim PDBModels that are their own source
        if isinstance( o, PDBModel ) and str( o.source ) == T.absfile( f ):
            o.forcePickle = 1
            
        T.flushPrint('Dumping %s\n' % f )
        T.dump( o, T.absfile( f ) )

    except:
        print "Error with ", f
        print T.lastError()
