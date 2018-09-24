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
Default Error Handler for Biskit classes.
"""

import biskit.tools as T
from biskit.logFile import ErrLog
from biskit.errors import HandledError, NormalError, FatalError

class ErrorHandler( object ):
    """
    Default Error Handler for Biskit classes.
    """

    def __init__( self, log=None ):
        """
        :param log: target of error messages, None->StdErr (default: None)
        :type  log: LogFile
        """
        self.log = log or ErrLog()


    def fatal( self, message ):
        """
        Handle a fatal error (likely a bug), stop program execution.

        :param message: message to be given to user
        :type  message: str

        :raise FatalError: 
        """
        s = '\nFatal Error: '+str(message)
        s += '\n\t' + T.lastError() + '\n'
        s += 'TraceBack: \n' + T.lastErrorTrace() + '\n'

        self.log.add(s)
        raise FatalError


    def error( self, message ):
        """
        Handle a normal error (like non-existing file) that is not
        necessarily a bug.

        :param message: message to be given to user
        :type  message: str

        :raise NormalError: 
        """
        s = '\nError: '+str(message)
        s += '\n\t' + T.lastError()
        s += '\nTraceBack: \n' + T.lastErrorTrace() + '\n'

        self.log.add(s)
        raise NormalError


    def warning( self, message, error=1, trace=0 ):
        """
        Issue a warning. No exception is raised.

        :param message: message to be given to user
        :type  message: str
        :param error: report Exception with line (default: 1)
        :type  error: 1||0
        :param trace: report full back trace to exception (default: 0)
        :type  trace: 1||0
        """
        s = '\nWarning (ignored): '+str(message)
        try:
            if trace: error = 1
            if error:
                s += '\n\t' + T.lastError() + '\n'
            if trace:
                s += '\nTraceBack: \n' + T.lastErrorTrace() + '\n'
        except:
            pass

        self.log.add(s)



#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """ErrorHandler test"""

    def prepare(self):
        import tempfile
        self.f_out = tempfile.mktemp( '_test_ErrorHandler' )

    def cleanUp(self):
        T.tryRemove( self.f_out )

    def test_ErrorHandler( self ):
        """ErrorHandler test"""
        from biskit import LogFile

        self.err_log = LogFile( self.f_out )

        self.e = ErrorHandler( log=self.err_log )
        self.e.warning( 'A warning' )

        if self.local:
            print('An error log file was written to %s'%self.f_out)

        lines = open(self.err_log.fname).readlines()
        self.assertEqual(lines[-1],'Warning (ignored): A warning\n')

if __name__ == '__main__':

    BT.localTest(debug=0)

