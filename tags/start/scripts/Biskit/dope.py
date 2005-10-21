#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Author$
## last $Date$


import os.path
import sys

import Biskit.tools as T
from Biskit.PDBModel import PDBModel
from Biskit.PDBDope import PDBDope

import Biskit.LocalPath

def _use( options ):
    print """

Syntax:	   dope.py -s sourceModel -i otherModels
                  [-so sourceOut -o othersPrefix -dic old_model_dic ]

Add conservation, accessibility profiles and foldX energies to a reference
model and models linking to this reference.

1) if sourceOut is given: Remove waters from source, put conservation score
   into profiles, saveAs sourceOut
2) update each of |otherModels| from their respective source, make |sourceOut|
   their new source, remove atoms (shouldn't be changed from |sourceOut|) and
   pickle them down to same file name plus |othersPrefix| if given.
3) update old model dic if given

Example 1:
dope.py -s ../../rec_wet/1B39.pdb -so ../../rec_wet/dry.model \
-i *.model -dic 1B39_model.dic
-> create source and update model.dic

Example 2:
dope.py -s ../../rec_wet/dry.model \
-i *.model -dic 1B39_model.dic
-> source already there, update model.dic


                     
Options:   
           
"""
    for key, value in options.items():
        print "\t-",key, "\t",value
    
    sys.exit(0)

class ConvertError(Exception):
    pass


def prepareSource( inFile, outFile ):
    """
    Strip waters, add profiles and save as doped source model.
    """
    
    source = PDBModel( inFile )
    
    source.remove( lambda a: a['residue_name'] in ['HOH','WAT','TIP3'] )
    
    source = source.sort()
    
    doper = PDBDope( source )

    doper.addASA()
    doper.addSurfaceMask()
    doper.addFoldX()
    doper.addSurfaceRacer( probe=1.4 )
    doper.addDensity()
    
    try:
        doper.addConservation( )
    except:
        errWriteln('\n ERROR: Conservation profile could not be added to '\
                   + str(sourceOut) + '\n' )
       
    source.saveAs( outFile )

    return source


def changeModel( inFile, prefix, sourceModel ):

    print '\nget ' + os.path.basename( inFile ) + '..',

    model = PDBModel( inFile )

    model.update()

    model = model.sort()

    eq = model.equals( sourceModel )
    if not eq[0] and eq[1]:
        raise ConvertError('source and other models are not equal: ' + str(eq))

#    model.validSource()
    model.setSource( sourceModel.validSource() )

    model.atomsChanged = 0
    model.xyzChanged = 1

    model.update( lookHarder=1 )

    ## accessib. could have changed from source
    doper = PDBDope( model )
    doper.addASA()
    doper.addSurfaceMask()
    doper.addSurfaceRacer( probe=1.4 )
    doper.addDensity()
    
    ## foldX energies change with coordinates, too
    doper.addFoldX()

    outFile = os.path.dirname( inFile ) + '/' + prefix +\
              T.stripFilename( inFile ) + '.model' 

    T.Dump( model, outFile )

    print '-> ' + os.path.basename( outFile )


def updateModelDic( f ):
    """ Call update() on all models in a dict to make them aware of the
    new profiles."""

    print 'Updating ', f

    d = T.Load( T.absfile( f ) )

    for m in d.values():
        m.update( lookHarder=1 )

    T.Dump( d, T.absfile( f ) )

##########
## MAIN ##

default = {'o':'' }

if len (sys.argv) < 2:
    _use( default )

options = T.cmdDict( default )

if options.has_key('i'):
    inLst = T.toList( options['i'] )
    inLst = [ T.absfile( f ) for f in inLst ]
else:
    inLst = []

sourceIn = T.absfile( options['s'] )

sourceOut = None
if options.has_key('so'):
    sourceOut = T.absfile( options['so'] )

prefix = options['o']

print 'Preparing source ' + str(os.path.basename(sourceIn))\
      + ' -> ' + str(sourceOut)

if sourceOut:
    source = prepareSource( sourceIn, sourceOut )
else:
    try:
        source = T.Load( sourceIn )
        if not source.profile('cons_ent', None):
            raise Exception()
    except:
        errWriteln('No -so given: -s must be PDBModel with cons. profiles')

print 'Changing other models '

for f in inLst:

    changeModel( f, prefix, source ) 

if options.has_key('dic'):

    print 'Updating model.dic'
    updateModelDic( options['dic'] )
