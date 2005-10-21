#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

from Biskit.tools import *
from Biskit import PDBModel

def _use():
    print """
    castPdbs: Convert two similar PDBs in two PDBs with equal atom content.
    PDBs must not have any HETATOMs. TIP3 residues are removed.

    Syntax: castPdbs.py -i1 |pdb1| -i2 |pdb2| -o1 |outFile1| -o2 |outFile2|
                        [ -c1 |0 1 ..| -c2 |2 3 ..| ]

    i1, i2   file names of PDBs to be compared
    o1, o2   file names for result pdbs
    c1, c2   chain numbers (starting 0) to take from i1 and i2 (default: all)
    
Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


if __name__ == '__main__':

    if len( sys.argv ) < 4:
        _use()

    options = cmdDict( {} )

    f1 = absfile( options['i1'] )
    f2 = absfile( options['i2'] )
    o1 = absfile( options['o1'] )
    o2 = absfile( options['o2'] )

    errWriteln("loading pdbs...")
    
    m1 = PDBModel( f1 )
    m2 = PDBModel( f2 )

    if options.has_key('c1'):
        m1 = m1.takeChains( toIntList( options['c1'] ) )

    if options.has_key('c2'):
        m2 = m2.takeChains( toIntList( options['c2'] ) )

    m1.removeRes( 'TIP3' )
    m2.removeRes( 'TIP3' )

    m1.sort()
    m2.sort()

    errWriteln("compare atoms of pdbs...")
    mask1, mask2 = m1.equalAtoms( m2 )

    errWriteln("removing %i atoms from %s" % (sum( logical_not( mask1 ) ), f1))
    m1.remove( logical_not( mask1 ) )
    
    errWriteln("removing %i atoms from %s" % (sum( logical_not( mask2 ) ), f2))
    m2.remove( logical_not( mask2 ) )

    errWriteln("writing new pdbs..." )
    m1.writePdb( o1 )
    m2.writePdb( o2 )
