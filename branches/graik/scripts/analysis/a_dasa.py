#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
import os, sys
from Numeric import *
from Biskit.tools import *
from Biskit import *

def use():
    print """Calculate change in access. and molecular surface upon binding.
    Syntax:   a_dasa.py -r |rec_model| -l |lig_model| -c |com_model|
    """
    sys.exit()

def area( m ):
    m.keep( nonzero( m.maskProtein()) )
    d = PDBDope( m )
    d.addSurfaceRacer( probe=1.4, vdw_set=2 )

    return sum( m.profile( 'AS' ) )

def MS( m ):
    return sum( m.profile('MS') )

## MAIN ##

if len( sys.argv ) < 2:
    use()

o = cmdDict({})

rec = Load(o['r'])
lig = Load(o['l'])
com = Load(o['c'])

diff_ASA = area( rec ) + area( lig ) - area( com )
diff_MS  = MS( rec ) + MS( lig ) - MS( com )

gamma_ASA = 46 ## cal/mol/A²; Sharp et al 1991, Noskov & Lim 2001
gamma_MS  = 69 ## cal/mol/A²; Jackson & Sternberg 1995
gamma_MS  = 45
T = 298        ## K

print "change in accessible surface: ", diff_ASA
print "change in molecular  surface: ", diff_MS

print 'hydrophobic "entropy" in cal/mol/K'

print 'based on ASA: ', diff_ASA * gamma_ASA / T
print 'based on MS : ', diff_MS  * gamma_MS / T
