##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
"""
Interface to Xplor
"""
## allow relative imports when calling module as main script for testing https://www.python.org/dev/peps/pep-0366/
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import tempfile
from time import time
import os

import biskit.settings as settings
import biskit.tools as T
from biskit.errors import BiskitError

from .executor import Executor, RunError

class XplorerError( BiskitError ):
    pass

class Xplorer(Executor):
    """
    Prepare and run a xplor job.
    
    Last tested with Xplor-NIH 2.45

    @see: L{Executor} for detailed documentation!

    @note: Command configuration in: biskit/Biskit/data/defaults/exe_xplor.dat
    """

    def __init__( self, template, pipes=None, **params):
        """
        Create an Executor instance for XPlor jobs. The only argument required
        is a template for the control script passed to Xplor. As usual, any
        keyword=value pairs given to this constructor, can be used to
        substitute place holders in the template (see L{Executor}!).
        
        @param template: template for Xplor input file. This can be the
                         template itself or the path to a file containing it.
        @type  template: str
        @param pipes: push input script directly through a pipe rather
                      than writing it to disc; If is not None, this option
                      overrides the default setting from exe_xplor.dat.
                      The default in exe_xplor is pipes = 1 (True).
        @type  pipes: 1|0 or bool

        @param kw: additional key=value parameters for Executor
        for example::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local)
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        Executor.__init__( self, 'xplor', template=template, **params )

        if pipes is not None:
            self.exe.pipes = pipes

        #: will contain copy of Xplor log file after run finished
        self.logLines = None


    def postProcess( self ):
        """
        Parse xplor log into self.logLines
        """
        self.logLines = open( self.f_out ).readlines()


    def isFailed( self ):
        """
        Calculation completion check.
        
        @return: if xplor died
        @rtype: 1
        """
        return self.logLines is None or \
               self.logLines[-1].find("X-PLOR: exit time") == -1


    def fail( self ):
        """
        called after aborted calculation, override!

        @raise RunError: if Xplor terminated with an error
        """
        s = ''
        try:
            s += ''.join( self.logLines[-7:] )
        except:
            pass

        raise RunError('Xplor terminated with an error:\n' + s)


    def saveXLog( self, fname ):
        """
        Save the content of the XPlor log file to a new file.
        
        @param fname: target file name
        @type  fname: str
        
        @raise XPlorerError: if logLines is empty (Xplor hasn't been run)
        """
        if not self.logLines:
            raise XplorerError('Xplor log file is (yet) empty.')

        f = open( T.absfile(fname), 'w')
        f.writelines( self.logLines )
        f.close()



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.EXE ]

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
    def prepare(self):
        self.dir_out = tempfile.mkdtemp( '_test_Xplorer' )

        ## write a test template inp file to disc
        f = open( self.dir_out +'/test.inp', 'w')
        f.writelines( self.inp )
        f.close()


    def test_Xplorer( self ):
        """Xplorer test"""

        ## input template variables
        param = T.dataRoot() + '/xplor/toppar/param19.pro'
        pdb_out = self.dir_out +'/lig.pdb'
        log_out = self.dir_out +'/test.out'
        pdb_in = T.testRoot() + '/lig/1A19.pdb' 
        psf_in = T.testRoot() + '/lig/1A19.psf'
        
        self.x = Xplorer( self.dir_out +'/test.inp',
                     xout = log_out,
                     verbose = self.local,
                     lig_psf = psf_in,
                     lig_pdb = pdb_in,
                     param19 = param,
                     lig_out = pdb_out )

        self.x.run()

        if self.local:
            print("""
The minimized structure and the X-Plor log file
have been written to %s and %s, respectively
            """ % (pdb_out, log_out))

            print("""See x.logLines for the complete xplor log file!""")
        

    def cleanUp(self):
        T.tryRemove( self.dir_out, tree=1 )

if __name__ == '__main__':

    BT.localTest()


