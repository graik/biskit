## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
Default Error Handler for Biskit classes.
"""

import Biskit.tools as T
from Biskit.LogFile import ErrLog
from Biskit.Errors import HandledError, NormalError, FatalError

class ErrorHandler( object ):
    """
    Default Error Handler for Biskit classes.
    
    Note: When the parameter "fails" is set to False then no exception is raised. 
    Its purpose is only for error hunting in test cases, which means that you want
    to see if an error occurred and this is not treated as a program failure.
    "fails" parameter musn't be used by a final user, only for those contributors
    who want to do their own test cases.
    
    """

    def __init__( self, log=None, verbose=True ):
        """
        @param log: target of error messages, None->StdErr (default: None)
        @type  log: LogFile
        @param verbose: report errors to log,if False, errors are only silently 
                        written to ErrorHandler.lastError and .lastWarning
                        (default: True)
        @type verbose: bool
        """
        self.log = log or ErrLog()
        self.lastError =  ""
        self.lastWarning = ""
        self.verbose = verbose   # report errors to log; set to False for silent mode
        self.fails = True

    def __reportException( self, error=True, trace=True ):
        """
        Re-format and report last exception in backtrace.
        @param error: report Exception with line (default: 1)
        @type  error: 1||0
        @param trace: report error trace of last exception (if any)
        @type  trace: bool
        @return: str, Line of error
        """
        r = ''
        try:
            if error:
                r += '\n\t' + T.lastError() + '\n'
            if trace:
                r += 'TraceBack: \n' + T.lastErrorTrace() + '\n'
        except AttributeError, err:
            ## no exception in stack
            pass
        return r
        
    def fatal( self, message, error=True, trace=True ):
        """
        Handle a fatal error (likely a bug), stop program execution.

        @param message: message to be given to user
        @type  message: str
        @param error: report Exception with line (default: 1)
        @type  error: 1||0
        @param trace: report error trace of last exception (if any)
        @type  trace: bool

        @raise FatalError: 
        """
        s = '\nFatal Error: '+str(message)
        s += self.__reportException( error=error, trace=trace )
        self.lastError = s
                
        if self.verbose:
            self.log.add(s)
        
        if fails:
            raise FatalError, s


    def error( self, message, error=True, trace=True ):
        """
        Handle a normal error (like non-existing file) that is not
        necessarily a bug.

        @param message: message to be given to user
        @type  message: str
        @param error: report Exception with line (default: 1)
        @type  error: 1||0
        @param trace: report error trace of last exception, if any (default: 1)
        @type  trace: bool

        @raise NormalError: 
        """
        s = '\nError: '+str(message)
        s += self.__reportException( error=error, trace=trace )
        self.lastError = s

        if self.verbose:
            self.log.add(s)
            
        if fails:
            raise NormalError, s


    def warning( self, message, error=1, trace=0 ):
        """
        Issue a warning. No exception is raised.

        @param message: message to be given to user
        @type  message: str
        @param error: report Exception with line (default: 1)
        @type  error: 1||0
        @param trace: report full back trace to exception (default: 0)
        @type  trace: 1||0
        """
        s = '\nWarning (ignored): '+str(message)
        s += self.__reportException( error=error, trace=trace )
        self.lastWarning = s

        if self.verbose:
            self.log.add(s)



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """ErrorHandler test"""

    def prepare(self):
        import tempfile
        self.f_out = tempfile.mktemp( '_test_ErrorHandler' )

    def cleanUp(self):
        T.tryRemove( self.f_out )

    def test_ErrorHandler( self ):
        """ErrorHandler test"""
        from Biskit.LogFile import LogFile

        self.err_log = LogFile( self.f_out )

        self.e = ErrorHandler( log=self.err_log )
        self.e.warning( 'A warning' )

        if self.local:
            print 'An error log file was written to %s'%self.f_out

        lines = open(self.err_log.fname).readlines()
        self.assertEquals(lines[-1],'Warning (ignored): A warning\n')

if __name__ == '__main__':

    BT.localTest(debug=0)


