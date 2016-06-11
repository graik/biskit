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

## 
## 

from Numeric import *
from Biskit.tools import *
from Biskit.Dock import ComplexList

def _use():
    print """
Reduce macro docked complex list (rec * lib * 5120)
to a normal complex list (rec * lig * 512)

   i - complexList to be reduced
   o - name of reduced list
 
    """
    default = _defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def _defOptions():
    return { 'o':'complexes_512.cl'}


def reduceList( c, rm=range(1,12), lm=range(1,12), soln=512 ):
    """
    Create an matrix: len(rec_model) * len((lig_model) * solutions
    containing the values of the info dic with given key.
    """
    c_new = ComplexList()

    flushPrint('Removing HEX solutions greater than 512.\n')
    c = c.filter( 'soln', (1, soln) )
    
    flushPrint('Sorting list ')
    for r in rm:
        flushPrint('|')
        cr = c.filter( 'model1', r )
        for l in lm:
            flushPrint('.')
            cl = cr.filter( 'model2', l )
            cl.sortBy('soln')
            c_new += cl
            
    return c_new


def test():
    options = _defOptions()
    dir = '/home/Bis/raik/data/tb/interfaces/c13/dock_pcr_multi_0919/hex0919'
    options['i'] = dir + '/complexes.cl'
    options['o'] = dir + '/complexes_reduced_512.cl'
    return options

###############
# main function

if __name__ == '__main__':
    if len(sys.argv) < 2:
        _use()
    
#    options = test()
    options = cmdDict( _defOptions() )

cFile = os.path.abspath( options['i'] )

flushPrint('Loading complex list \n')
cList = load( options['i'] )

flushPrint('Reducing complex list \n')
new_cList = reduceList( cList )

dump( new_cList, absfile(options['o']) )

