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

## Start XPlor MD for interface project

from os.path import *
import os
import shutil
import commands
import time
import sys
from Biskit.tools import *

modules = os.path.dirname(testRoot()) +'/'
sys.path.insert(0, modules)

#from Biskit.tools import *
from Biskit import BiskitError
import Biskit.settings as settings
from Biskit.Errors import BiskitError


options = {'f':0,
           'r':'.',
           'n':0,
           'h':'localhost',
           'i':dataRoot() + '/xplor/inpPcr',
           'parm':dataRoot() + '/xplor/toppar',
           'nstep':'500',
          }

def _use():
    print """
Start single Xplor PCR job on remote host.

Syntax: runPCR -t |psfFolder| -h |host|
\t\t[-f |Force| -r |resultFolder| -n |nice| -i |inpFolder| ]

Options:
\t-f\tforce constant for PCR restraint
\t-r\tbase folder for result (sub-folder will be created)
\t-t\tfolder with topology (psf, pdb)
\t-n\tnice value
\t-h\thost computer (accessible via ssh)
\t-i\tfolder with all input files, must contain restart_h2o.inp'
\t-parm\tfolder with param19.* files

A MD folder called pcr_<PDBCode> is created. The force constant is
written to a file 'oldenergy' in this folder. The topology folder is
copied to the new md folder and renamed to the PDB code
(which is taken from the first 4 letters of the psf file name).
The job is started via ssh on the remote host.
A summary of all used parameters is written to runReport.out.

NOTE: The pdb/psf file has to be 4 characters long and start with
      a number (the segId has to conform to the same format).

By using python -i runPCR.py .. a python shell remains open and
the job can be killed with the command r.kill()

Default options:"""
    for key in options.keys():
        print "\t-",key, "\t",options[key]

    sys.exit(0)


class Runner:
    """
    Prepare, run and manage single Xplor-PCR job
    """

    def __init__(self, force, resultFolder, psfFolder,
                 host='localhost', nice=0,
                 inpFolder=dataRoot() + '/xplor/inpPcr',
                 **options ):
        """
        force - int, force constant for PCR-Runs
        resultFolder - string, folder in wich a sub folder named
                       pcr_00 or pcr_|force| will be created
        psfFolder    - string, folder with pdb and psf of the system
        host         - string, host computer to be used
        nice         - int, nice level
        inpFolder    - string, folder with inp files needed
        parm         - str, path to toppar folder with param19*
        """
        self.force = force
        self.fResult = absfile( resultFolder )
        self.fTopo = absfile( psfFolder )
        self.fInp  = absfile( inpFolder )
        self.host  = host
        self.nice  = nice
        self.mdFolder = ''    ## folder to start XPlor in

        ## take additional parameters into self
        self.__dict__.update( options )

        ## XPlor binary
        self.xbin  = settings.xplorpcr_bin

        ## will hold PID of xplor job
        self.pid = 0

        self._prepareFolder()
        

    def _getPdbCode( self ):
        """Look in folder for file 1ABCxyz.psf, return first 4 letters.
        -> string
        """
        files = os.listdir( self.fTopo )

        for f in files:
            if f[-4:] == '.psf':
                return f[:4]

        raise ProgError('no psf file in topology folder!')
    

    def _prepareFolder(self):
        """
        Create folder for run, copy over inp files and topology.
        """

        pdbCode = self._getPdbCode( )
        
        newTopoFolder = self.fResult + '/' + pdbCode
        self.mdFolder = self.fResult + '/pcr_%02i' % self.force

        try:

            ## create folder for trajectory
            status = self.mdFolder
            os.mkdir( self.mdFolder )

            ## copy over topology
            status = newTopoFolder
            if not exists( newTopoFolder ):
                shutil.copytree( self.fTopo, newTopoFolder )

            ## copy over input files
            status = self.fInp
            inpFiles = os.listdir( self.fInp )

            try: inpFiles.remove('CVS')
            except: pass

            for f in inpFiles:
                target = self.mdFolder + '/' + f
                shutil.copy( self.fInp + '/' + f, self.mdFolder )
                self.__fillTemplate( target, target, self.__dict__ )

            ## create force constant inp file
            f = open( self.mdFolder + '/initialpcr','w')
            f.write( str(self.force) )
            f.close()

            ## create force constant inp file for restarts
            ## (not needed for normal start..)
            f = open( self.mdFolder + '/oldenergy', 'w')
            f.write( str(self.force) )
            f.close()

            ## write PDB code to file
            f = open( self.mdFolder + '/segid', 'w' )
            f.write( '"%s"\n' % self._getPdbCode()  )
            f.close()

        except:
            errWriteln('Error while preparing directory for MD.')
            errWriteln('Working on: ' + status)
            errWrite( lastError() )


    def _setPID( self ):

        time.sleep(2)

        try:
            cmd = "ssh %s '/bin/ps -C %s -o pid='" % \
                  (self.host, stripFilename( self.xbin ) )

            lines = commands.getoutput( cmd )

            lines = lines.split()

            self.pid = int( lines[-1] )
            
        except:
            errWriteln("Couldn't get PID with %s" % cmd )
            errWriteln("Problem: " + str( lastError()) )


    def runXplor( self ):
        xinp = 'restart_h2o.inp'
        xout = 'restart_h2o.out'

        try:
            f = open( self.mdFolder + '/host.out' , 'w')
            f.write( self.host )
            f.close()
        except:
            errWriteln("Error writing host name to %s." % self.mdFolder )
            errWriteln("Problem: " + str( lastError() ) )

##         cmd = "ssh %s 'cd %s; /bin/nice -%i %s ' & "\
##               %(self.host, self.mdFolder, self.nice, self.xbin)
        if self.nice > 0:
            nice_cmd = '/bin/nice -n %i' % self.nice
        else:
            nice_cmd = ''
        
        cmd = "ssh %s 'cd %s; %s %s < %s > %s &' & "\
              %(self.host, self.mdFolder, nice_cmd, self.xbin, xinp, xout)
        
        try:
            ## Problem: job gets automatically niced down to 19, zsh-related
            print "executing: ", cmd, "\n"
            os.system( cmd )

        except:
            errWriteln("Error executing command " + cmd)
            errWriteln("Problem: " + str( lastError() ) )

        self._setPID()
        

    def __fillTemplate(self, fTemplate, fout, valueDic):
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


    def kill( self ):

        try:
            cmd = "/usr/bin/ssh %s '/bin/kill -9 %i'" % (self.host, self.pid )
            
            if os.system( cmd ) <> 0:
                raise ProgError("non-0 exit status. \nCommand: "+cmd)

        except:
            print lastError()


    def report( self ):

        s = 'host: %s\n' % self.host
        s += 'pid:  %s\n' % self.pid
        s += 'nice: %i\n' % self.nice
        s += 'binary:        %s\n' % self.xbin
        s += 'result folder: %s\n' % self.mdFolder
        s += 'topology:      %s\n' % self.fTopo
        s += 'force:         %i\n' % self.force
        
        return s
            

        
###########################
# MAIN
###########################

if len(sys.argv) < 3:
    _use()

op = cmdDict( options )

r = Runner( int( op['f'] ), op['r'], op['t'], host=op['h'], nice=int(op['n']),
            inpFolder=op['i'], nstep=op['nstep'], parm=op['parm'] )

r.runXplor()

print r.report()

try:
    fn = r.mdFolder + '/runReport.out'
    f = open(fn, 'w') 
    f.write( r.report() )
    f.close()
except:
    errWriteln('Error creating ' + fn )
    errWriteln( str( lastError() ) )
