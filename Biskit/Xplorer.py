##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## $Revision$
## last $Date$
## last $Author$

import tempfile
from time import time
import os

import settings
import Biskit.tools as t
from Biskit import BiskitError
from Biskit import StdLog

class XplorerError( BiskitError ):
    pass

class TemplateError( XplorerError ):
    pass

class RunError( XplorerError ):
    pass

class Xplorer:
    """Prepare and run a xplor job
    @todo:  Most of this has been abstracted into Executor
            -> subclass Executor and delete redundant parts
    """

    def __init__( self, inpFile, xout=None, bin=settings.xplor_bin,
                  node=None, nice=0, log=None, debug=0, verbose=None,
                  **params ):
        """
        inpFile - str, file name of inp template
        xout    - str, file name of xplor log file (None->discard) [None]
        bin     - str, file name of xplor binary     [settings.xplor_bin]
        node    - str, host for calculation (None->local)          [None]
        nice    - int, nice level                                     [0]
        log     - Biskit.LogFile, program log (None->STOUT)        [None]
        debug   - 0|1, keep all temporary files                       [0]
        verbose - 0|1, print xplor command to log         [log != STDOUT]
        [ additional key=value pairs for inp file ]
        """

        self.__dict__.update( params )

        self.bin = t.absfile( bin ) or settings.xplor_bin
        self.inp = t.absfile( inpFile )
        self.xout = t.absfile( xout ) or tempfile.mktemp('_xplor.log')
        self.finp = None
        self.node = node or os.uname()[1]   ##.split('.')[0]
        self.nice = nice
        self.keepLog = xout is not None
        self.debug = debug

        ## Log object for own program messages
        self.log = log or StdLog()
        self.verbose = verbose
        if self.verbose is None:
            self.verbose = (log is not None)

        ## time needed for last run
        self.runTime = 0

        ## version as of creation of this object
        self.initVersion = self.version()

        ## will contain Xplor log file after run finished
        self.logLines = None


    def version( self ):
        return 'Xplorer $Revision$'


    def runXplor( self, finp ):
        """
        Run Xplor with given inp file and wait until it is finished.
        Called by run().
        finp - str, name of existing input file (w/o place holders)
        -> float, run time in s
        """
        start_time = time()

        str_nice = "%s -%i" % (settings.nice_bin, self.nice)

        lead_cmd = settings.ssh_bin
        str_ssh = "%s %s" % ( settings.ssh_bin, self.node )
        
        cmd = "%s %s %s < %s > %s"\
              %(str_ssh, str_nice, self.bin, finp, self.xout)

        if self.verbose: self.log.add_nobreak("executes: "+cmd+"..")

        if os.spawnvp(os.P_WAIT, lead_cmd, cmd.split() ) <> 0:
            raise RunError("non-0 exit status. \nCommand: "+cmd)

        if self.verbose: self.log.add(".. finished.")

        return time() - start_time


    def run( self, inp_mirror=None ):
        """
        calls (in that order): prepare(), runXplor(), parseLog(),
        finish()/failed(), cleanup()
        inp_mirror - str, file name for formatted copy of inp file [None]
        """

        self.prepare()

        try:
            inp = open( self.inp, 'r' ).read()

            self.finp = inp_mirror or tempfile.mktemp('_xplor.inp')

            self.generateInp( inp, self.finp )

            self.runTime = self.runXplor( self.finp )

            self.parseLog()

        except RunError, why:
            try:
                self.failed()
            finally:
                self.cleanup()
            raise RunError, why

        try:
            if self.isFailed():
                self.failed()
            else:
                self.finish()
        finally:
            self.cleanup()
        

    def cleanup( self ):
        """Remove temporary files, override it but call it in the child"""
        if not self.keepLog and not self.debug:
            t.tryRemove( self.xout )
            
        if not self.debug:
            t.tryRemove( self.finp )


    def parseLog( self ):
        """Parse xplor log into self.logLines"""
        self.logLines = open( self.xout ).readlines()


    def isFailed( self ):
        """-> 1, if xplor died"""
        return self.logLines == None or \
               self.logLines[-1].find("X-PLOR: exit time") == -1

    def prepare( self ):
        """called before running xplor, override!"""
        pass

    def finish(self ):
        """called after normal termination, override!"""
        pass

    def failed( self ):
        """called after aborted calculation, override!"""
        s = ''
        try:
            s += ''.join( self.logLines[-7:] )
        except:
            pass
                
        raise RunError, 'Xplor terminated with an error:\n' + s

    def generateInp(self, inp, fout):
        """
        Replace formatstr place holders in inp by fields of this class.
        inp - str, content of the input file with place holders
        -> str, complete inp file
        !! raise TemplateError
        """
        try:

            s = inp % self.__dict__
            f = open( fout, 'w')
            f.write(s)
            f.close()
            
        except KeyError, why:
            s =  "Unknown option/place holder in template file."
            s += "\n  template file: " + self.inp
            s += "\n  Template asked for a option called " + str( why[0] )
            raise TemplateError, s
        
        except Exception, why:
            s =  "Error while creating template file."
            s += "\n  template file: " + self.inp
            s += "\n  why: " + str( why )
            s += "\n  Error:\n  " + t.lastError()
            raise TemplateError, s


    def saveXLog( self, fname ):
        """
        Save the content of the XPlor log file to a new file.
        fname - str, target file name
        !! XPlorerError, if logLines is empty (Xplor hasn't been run)
        """
        if not self.logLines:
            raise XPlorerError('Xplor log file is (yet) empty.')

        f = open( t.absfile(fname), 'w')
        f.writelines( self.logLines )
        f.close()

        

### TEST #####

if __name__ == '__main__':

    pass
