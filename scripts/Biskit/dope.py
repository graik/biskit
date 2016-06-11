#!/usr/bin/env python
## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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


import os.path
import sys
import Biskit.oldnumeric as N0

import Biskit.tools as T
from Biskit.PDBModel import PDBModel
from Biskit.PDBDope import PDBDope

import Biskit.LocalPath

def _use( options ):
    print """

Syntax:    dope.py -s sourceModel -i otherModels [-p [fx surf dens cons]]
                  [-so sourceOut -o othersPrefix -dic old_model_dic ]
                  [-nosort -wat]

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

   -s      source input PDB or pickled PDBModel
   -p      profiles to be calculated:
             fx   ... foldx energies from FoldX (not a real profile)
             surf ... surfrace accessibilities and curvature fom surfrace
             dens ... atomic densities
             cons ... sequence conservation from HMM
             dssp ... secondary structure from DSSP
             delphi.. info record (no profile): DelPhi electrostatic potential
   -so     filename of updated (source) model pickle
   -i      PDBModels that should be linked to updated source
   -o      pickle updated -i models with this file name prefix
   -dic    file name of pickled model dict to be updated
   -nosort do not sort atoms within residues
   -wat    keep waters

Default options:   

"""
    for key, value in options.items():
        print "\t-",key, "\t",value

    sys.exit(0)


class ConvertError(Exception):
    pass


def prepareSource( inFile, outFile, wat=1, sort=1,
                   foldx=1, surf=1, dens=1, cons=1, dssp=1, delphi=0 ):
    """
    Strip waters, add profiles and save as doped source model.
    """

    source = PDBModel( inFile )

    if wat:
        source.remove( lambda a: a['residue_name'] in ['HOH','WAT','TIP3'] )

    if sort:
        source = source.sort()

    doper = PDBDope( source )

    if surf:
##         doper.addASA()
##         doper.addSurfaceMask()
        doper.addSurfaceRacer( probe=1.4 )

    if foldx:
        doper.addFoldX()

    if dens:
        doper.addDensity()

    if dssp:
        doper.addSecondaryStructure()
        
    if delphi:
        doper.addDelphi()

    try:
        if cons:
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

    #model.atomsChanged = 0
    for k in model.atoms:
        model.atoms[k,'changed'] = N0.all( model[k] == sourceModel[k] )

    model.xyzChanged = ( 0 != N0.sum( N0.ravel( model.xyz - sourceModel.xyz)) )

    model.update( updateMissing=1 )

    if model.xyzChanged:

        doper = PDBDope( model )

        if 'MS' in sourceModel.atoms.keys():
            doper.addSurfaceRacer( probe=1.4 )

        if 'density' in sourceModel.atoms.keys():
            doper.addDensity()

        if 'foldX' in sourceModel.info.keys():
            doper.addFoldX()
            
        if 'delphi' in sourceModel.info.keys():
            doper.addDelphi()

    outFile = os.path.dirname( inFile ) + '/' + prefix +\
            T.stripFilename( inFile ) + '.model' 

    T.dump( model, outFile )

    print '-> ' + os.path.basename( outFile )


def updateModelDic( f ):
    """ Call update() on all models in a dict to make them aware of the
    new profiles."""

    print 'Updating ', f

    d = T.load( T.absfile( f ) )

    for m in d.values():
        m.update( updateMissing=1 )

    T.dump( d, T.absfile( f ) )

##########
## MAIN ##

default = {'o':'', 'p':'surf dens'}

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
    source = prepareSource( sourceIn, sourceOut,
                            wat  =('wat'  not in options),
                            sort =('nosort' not in options),
                            foldx=('fx'   in options['p']),
                            surf =('surf' in options['p']),
                            dens =('dens' in options['p']),
                            cons =('cons' in options['p']),
                            dssp =('dssp' in options['p']),
                            delphi=('delphi' in options['p']))
else:
    try:
        source = T.load( sourceIn )
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
