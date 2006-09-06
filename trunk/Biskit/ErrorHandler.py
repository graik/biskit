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
Default Error Handler for Biskit classes.
"""
    
import Biskit.tools as tools
from Biskit.LogFile import ErrLog
from Biskit.Errors import HandledError, NormalError, FatalError

class ErrorHandler( object ):
    """
    Default Error Handler for Biskit classes.
    """

    def __init__( self, log=None ):
        """
        @param log: target of error messages, None->StdErr (default: None)
        @type  log: LogFile
        """
        self.log = log or ErrLog()


    def fatal( self, message ):
        """
        Handle a fatal error (likely a bug), stop program execution.
        
        @param message: message to be given to user
        @type  message: str
        
        @raise FatalError: 
        """
        s = '\nFatal Error: '+str(message)
        s += '\n\t' + tools.lastError() + '\n'
        s += 'TraceBack: \n' + tools.lastErrorTrace() + '\n'

        self.log.add(s)
        raise FatalError


    def error( self, message ):
        """
        Handle a normal error (like non-existing file) that is not
        necessarily a bug.
        
        @param message: message to be given to user
        @type  message: str
        
        @raise NormalError: 
        """
        s = '\nError: '+str(message)
        s += '\n\t' + tools.lastError()
        s += '\nTraceBack: \n' + tools.lastErrorTrace() + '\n'

        self.log.add(s)
        raise NormalError


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
        try:
            if trace: error = 1
            if error:
                s += '\n\t' + tools.lastError() + '\n'
            if trace:
                s += '\nTraceBack: \n' + tools.lastErrorTrace() + '\n'
        except:
            pass

        self.log.add(s)



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0  ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: something
        @rtype:  float
        """
        import tempfile
        from Biskit.LogFile import LogFile
        
        f_out = tempfile.mktemp( '_test_ErrorHandler' )
        err_log = LogFile( f_out )
        
        
        self.e = ErrorHandler( log=err_log )

        self.e.warning( 'A warning' )

        if local:
            print 'An error log file was written to %s'%f_out
            globals().update( locals() )
            
        return  1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: something
        @rtype:  float
        """
        return 1

        
if __name__ == '__main__':

    test = Test()

    assert  test.run( local=1 ) == test.expected_result()
    

