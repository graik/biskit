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

"""
Prepare and run a xplor job.
"""

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
    """
    Prepare and run a xplor job

    @todo:  Most of this has been abstracted into Executor next thing
            to do is to subclass Executor and delete redundant parts 
    """

    def __init__( self, inpFile, xout=None, bin=settings.xplor_bin,
                  node=None, nice=0, log=None, debug=0, verbose=None,
                  **params ):
        """
        @param inpFile: file name of inp template
        @type  inpFile: str
        @param xout: file name of xplor log file (None->discard)
                     (default: None)
        @type  xout: str
        @param bin: file name of xplor binary (default: settings.xplor_bin)
        @type  bin: str
        @param node: host for calculation (None->local) (default: None)
        @type  node: str
        @param nice: nice level (default: 0)
        @type  nice: int
        @param log: Biskit.LogFile, program log (None->STOUT) (default: None)
        @type  log: str OR None
        @param debug: keep all temporary files (default: 0)
        @type  debug: 0|1
        @param verbose: print xplor command to log         [log != STDOUT]
        @type  verbose: 0|1
        @param params: additional key=value pairs for inp file
        @type  params: key=value
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
        """
        Version of class.
        
        @return: version
        @rtype: str
        """       
        return 'Xplorer $Revision$'


    def runXplor( self, finp ):
        """
        Run Xplor with given inp file and wait until it is finished.
        Called by run().
        
        @param finp: name of existing input file (w/o place holders)
        @type  finp: str
        
        @return: run time in s
        @rtype: float
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
        Run the callculation. This calls (in that order):
          - L{ prepare() },
          - L{ runXplor() },
          - L{ parseLog() },
          - L{ finish() }/ L{ failed() },
          - L{ cleanup() }
          
        @param inp_mirror: file name for formatted copy of inp file
                           (default: None)
        @type  inp_mirror: str
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
        """
        Remove temporary files, override it but call it in the child
        """
        if not self.keepLog and not self.debug:
            t.tryRemove( self.xout )

        if not self.debug:
            t.tryRemove( self.finp )


    def parseLog( self ):
        """
        Parse xplor log into self.logLines
        """
        self.logLines = open( self.xout ).readlines()


    def isFailed( self ):
        """
        Calculation completion check.
        
        @return: if xplor died
        @rtype: 1
        """
        return self.logLines == None or \
               self.logLines[-1].find("X-PLOR: exit time") == -1


    def prepare( self ):
        """
        called before running xplor, override!
        """
        pass


    def finish(self ):
        """
        called after normal termination, override!
        """
        pass


    def failed( self ):
        """
        called after aborted calculation, override!

        @raise RunError: if Xplor terminated with an error
        """
        s = ''
        try:
            s += ''.join( self.logLines[-7:] )
        except:
            pass

        raise RunError, 'Xplor terminated with an error:\n' + s


    def generateInp(self, inp, fout):
        """
        Replace formatstr place holders in inp by fields of this class.
        
        @param inp: content of the input file with place holders
        @type  inp: str
        @param fout: output file mane (xplor input file)
        @type  fout: str
        
        @return: complete inp file
        @rtype: str
        
        @raise TemplateError: if unknown option/place holder in template file
        @raise TemplateError: if error while creating template file
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
        
        @param fname: target file name
        @type  fname: str
        
        @raise XPlorerError: if logLines is empty (Xplor hasn't been run)
        """
        if not self.logLines:
            raise XPlorerError('Xplor log file is (yet) empty.')

        f = open( t.absfile(fname), 'w')
        f.writelines( self.logLines )
        f.close()



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """

    inp = """
!! short test  minimzation 

! ------------------------------------------------------------
! Place holders to be inserted by Python script
! ------------------------------------------------------------
evaluate ($ligandpsf = "%(lig_psf)s")
evaluate ($ligandpdb = "%(lig_pdb)s")
evaluate ($param19 = "%(param19)s" )
evaluate ($lig_out = "%(lig_out)s" )

! -------------------------------------
! read psf and pdb files for the ligand
! -------------------------------------
structure @@$ligandpsf end
coor @@$ligandpdb
delete selection= (resname TIP3) end

! ------------------
! set toplogies etc.
! ------------------
parameter
  	reset
        @@$param19
end

! ------------------------------------------------------------
! minimize
! ------------------------------------------------------------
flags exclude * include bond angle impr elec end
minimize powell nstep=10            end

! ---------------------------------------------------
! write minimized PDBs
! ---------------------------------------------------
write coor output= $lig_out end
stop
"""
    
    def run( self, verbose=0 ):
        """
        run function test

        @return: 
        @rtype:  
        """
        dir_out = tempfile.mkdtemp( '_test_Xplorer' )

        ## write a test template inp file to disc
        f = open( dir_out +'/test.inp', 'w')
        f.writelines( self.inp )
        f.close()

        ## input template variables
        param = t.projectRoot() + '/external/xplor/toppar/param19.pro'
        pdb_out = dir_out +'/lig.pdb'
        log_out = dir_out +'/test.out'
        pdb_in = t.testRoot() + '/lig/1A19.pdb' 
        psf_in = t.testRoot() + '/lig/1A19.psf'
        
        x= Xplorer( dir_out +'/test.inp',
                    xout = log_out,
                    verbose = verbose,
                    lig_psf = psf_in,
                    lig_pdb = pdb_in,
                    param19 = param,
                    lig_out = pdb_out )

        x.run()
        
        print 'The minimized structure and the X-Plor log file has been written to %s and %s, respectively'%(pdb_out, log_out)
        
        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 
        @rtype: 
        """
        return 1

        

if __name__ == '__main__':

    test = Test()

    assert test.run( verbose=1 ) == test.expected_result()



