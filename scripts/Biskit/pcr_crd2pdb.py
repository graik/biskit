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

## extract pdb files from xplor PCR trajectories

from Biskit.tools import *
from Biskit import ExeConfig
import Biskit.settings as settings

import os
from os.path import *
import tempfile
import shutil
import commands

options = {'i':'.',
           'o':'.',
           'n_iter':'50',
           'skip':'250',
          }

def _use():
    print """
Extract PDB files from Xplor PCR trajectories.

Syntax: pcr_crd2pdb -i |pcrFolder| -t |psfFolder| -o |outFolder|
                    [ -n_iter |n_iterations| -skip |stepping| -z ]

\t-z       gzip crd files
\t-skip    MD-step intervall for PDBs (500 = 1/ps)
\t-n_iter  number of iterations per ensemble member (50 = 50ps)

Default options:"""
    for key in options.keys():
        print "\t-",key, "\t",options[key]

    sys.exit(0)



class Runner:

    def __init__(self, fin, ftopo, fout, iterations=50, skip=250, zip=0):

        self.pcr_folder = '"%s/"' % absfile( fin )
        self.topo_folder =  '"%s/"' % absfile( ftopo )
        self.out_folder =  '"%s/"' % absfile( fout )

        self.n_iteration = iterations
        self.skip = skip

        if (skip % 10 != 0 ):
            print "WARNING: skip should be multiple of 10 (assuming that" +\
                  " the crd files contain each 10th MD step)."

        self.pdbCode = '"%s"' % self.getPdbCode()

        ## XPlor binary
        self.xbin = ExeConfigCache.get( 'xplor' ).bin
##        self.xbin  = settings.xplorpcr_bin

        ## template for xplor input
        self.ftemplate = absfile( dataRoot() + '/xplor/pcr_crd2pdb.inp')

        self.afterzip = zip


    def run( self ):

        fout = tempfile.mktemp('crd2pdb_inp')
        self.fillTemplate( self.ftemplate, fout, self.__dict__ )

        print "unzipping crd ..."
        self.unzip()

        cmd = "%s < %s > %s " % (self.xbin, fout, 'xplor.out')

        try:
            print "executing: ", cmd, "\n"
            os.system( cmd )

        except:
            errWriteln("Error executing command " + cmd)
            errWriteln("Problem: " + str( lastError() ) )

        if self.afterzip:
            print "zipping crd ..."
            self.zip()

        os.remove( fout )
            

    def getPdbCode( self ):
        """Hack, assume psf folder contains CODE.psf"""

        files = os.listdir( self.topo_folder.replace('"','') )

        for f in files:
            if ( len( f ) == 8) and f[-4:].upper() == '.PSF':
                return f[:4]

    def unzip( self ):
        """Unzip files needed by xplor script."""

        f = self.pcr_folder.replace('"','')
        try:
            cmd = 'gunzip %s*crd.gz; gunzip %satomselection.gz' % (f, f )
            os.system( cmd )
        except:
            errWriteln('Error unzipping: '+ lastError())


    def zip( self ):
        f = self.pcr_folder.replace('"','')
        try:
            cmd = 'gzip %s*crd' % ( f )
            os.system( cmd )
        except:
            errWriteln('Error zipping: ' + lastError())
        

    def fillTemplate(self, fTemplate, fout, valueDic ):
        """
        Read template input file with formatstr placeholders, insert values
        from valueDic. Write it back.
        fTemplate - String, filename for template
        valueDic  - Dictionary, {placeHolder:value}
        """
        try:

            result = []
            line = None
            for line in open( fTemplate ):
                result += line % valueDic

            f = open( fout, 'w' )
            f.writelines( result )
            f.close()

        except KeyError, why:
            s =  "Unknown option in template file."
            s += "\n  template file: " + fTemplate
            s += "\n  Template asked for a option called " + str( why[0] )
            s += "\n  template line:\n  " + str(line)
            s += "\n  Please give a value for this option at the command line."
            s += "\n  E.g: runAmber.py -parm ... -%s some_value -out..." %\
                 str( why[0] )

            raise BiskitError, s
        
        except:
            s =  "Error while adding template file."
            s += "\n  template file: " + fTemplate
            s += "\n  output file: " + fout
            s += "\n  template line:\n  " + str(line)
            s += "\n  available arguments:\n"

            for i in valueDic.keys():
                s += "\t%25s\t%s\n" % (i, str( valueDic[i] ) )
                
            s += "\n  Error:\n  " + lastError()

            raise BiskitError, s

###########################
# MAIN
###########################

if len(sys.argv) < 2:
    _use()

op = cmdDict( options )

r = Runner( op['i'], op['t'], op['o'], int( op.get('n_iter', 50)),
            int( op.get('skip', 250)), zip=op.has_key('z') )

r.run()

