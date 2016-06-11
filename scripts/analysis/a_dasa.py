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

rec = load(o['r'])
lig = load(o['l'])
com = load(o['c'])

diff_ASA = area( rec ) + area( lig ) - area( com )
diff_MS  = MS( rec ) + MS( lig ) - MS( com )

gamma_ASA = 46 ## cal/mol/A2; Sharp et al 1991, Noskov & Lim 2001
gamma_MS  = 69 ## cal/mol/A2; Jackson & Sternberg 1995
gamma_MS  = 45
T = 298        ## K

print "change in accessible surface: ", diff_ASA
print "change in molecular  surface: ", diff_MS

print 'hydrophobic "entropy" in cal/mol/K'

print 'based on ASA: ', diff_ASA * gamma_ASA / T
print 'based on MS : ', diff_MS  * gamma_MS / T
