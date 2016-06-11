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
# last $Author$
# last $Date$

from Biskit import PCRModel
from Biskit.tools import *
from Biskit.Dock.hexTools import *

import os
import sys
from types import *
import tempfile

def _use():
    print """
pcr2hex  pool many pdb's into one seperated by MODEL/ENDMDL to be used by hex.

Syntax:  pcr2hex -psf |in.psf| -pdb |in1.pdb| |in2.pdb| ... [-s |modelFolder| ]

         -psf     psf file name
         -pdb     list of pdb file names
         -nos     don't pickle each PDB as pickled w/o waters to this folder

Result:  -pdb file, pdbCode_hex.pdb (first 4 characters of the first pdb file
          are taken as pdbCode)
         -model dictionary, pdbCode_models.dic
         -modelFolder/in1.model, in2.model, unless -nos has been given
"""
    sys.exit(0)


def loadModels( psfName, pdbNames ):
    """
    psfName - str, psf file name, same for all pdbs
    filenames - [ str, str, ..], set of pdb filenames
    filenames with ~, ~user, ../ are expanded to absolute names

    -> {1:PCRModel, 2:PCRModel, 3:..}, model dictionary
    """
    modelDic = {}
    
    ## if only one pdb file, instead of list
    if type(pdbNames) is StringType:
        pdbNames = [pdbNames]
    
    counter = 0
    for f in  pdbNames:

        counter += 1
        model = PCRModel( absfile(psfName), absfile(f) )   

        # chain id removed by xplor, must be identical to ref complex
        model.addChainFromSegid() 

        model.removeRes( 'TIP3' )

        modelDic[ counter ] = model

    return modelDic


## def writeHexPdb( modelDic, pdbName):
##     """
##     write pdb for hex with models separated by MODEL%i/ENDMODEL
##     """
##     ## open new file
##     out = open(pdbName, 'w')
##     ## fetch name for temporary file
##     tmpFile = tempfile.mktemp('_pcr2hex.pdb')  
    
##     numbers = modelDic.keys()
##     numbers.sort()
##     for n in numbers:

##         ## write temporary pdb and open it
##         modelDic[ n ].writePdb( tmpFile )
        
##         pdb_temp = open( tmpFile )
        
##         out.write("MODEL%6i\n" % n)
##         lines = pdb_temp.readlines()    # get all lines
##         for line in lines:
##             if line[:3] <> 'END':       # filter out END
##                 out.write(line)
##         out.write("ENDMDL\n")

##         ## remove temporary file
##         pdb_temp.close()
##         os.remove( tmpFile )

##     out.write("END")
##     out.close()


def writeModelDic( modelDic, pdbcode, noSave=0 ):
    """
    !! PDBModelError if there is no source file
    """

    for m in modelDic.values():

        ## save if requested
        if not noSave:
            fname = stripFilename( m.sourceFile() ) + '.model'
            m.saveAs( absfile(fname) )
                     
    ## pickle dictionary
    dumpName =  pdbcode + '_model.dic'    # model dump file name
    dump(modelDic, dumpName)    # dump all models to file
    

############################################
# MAIN

def test():
    """Insert for testing"""
    s = '~johan/dock/scripts/test_complex11/'
    options = get_cmdDict( ['-psf',s+'rec/1A2P.psf',
                            '-pdb',s+'pcr_rec/1A2Pfit_1.pdb',
                            s+'pcr_rec/1A2Pfit_2.pdb'], {})
    return options

if __name__ == '__main__':
    if len(sys.argv) < 3:
        _use()
    
    import time
    start = time.time()
    
    options = cmdDict({})
#    options = test()
    
    ## create PCRModels and put them into dictionary
    modelDic = loadModels( options.get('psf',None), options['pdb'] )
    
    ## fetch PDB code (by default first 4 letters of file name)
    pdbCode = modelDic[1].getPdbCode()

    ## write hex pdb file
    createHexPdb( modelDic, pdbCode + '_hex.pdb')

    noSave = options.has_key( 'nos' )

    writeModelDic( modelDic, pdbCode, noSave )

    print "done in " + str( time.time() - start ) +"s"
