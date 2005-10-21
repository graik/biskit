#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## 
## 
## last $Author$
## last $Date$
## $Revision$

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
cList = Load( options['i'] )

flushPrint('Reducing complex list \n')
new_cList = reduceList( cList )

Dump( new_cList, absfile(options['o']) )

