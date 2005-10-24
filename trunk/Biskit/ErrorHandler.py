## last $Author$
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
## last $Date$
## $Revision$

import Biskit.tools as tools
from Biskit.LogFile import ErrLog
from Biskit.Errors import HandledError, NormalError, FatalError

class ErrorHandler( object ):
    """
    Default Error Handler for Biskit classes.
    """

    def __init__( self, log=None ):
        """
        log - LogFile, target of error messages, None->StdErr [None]
        """
        self.log = log or ErrLog()


    def fatal( self, message ):
        """
        Handle a fatal error (likely a bug), stop program execution.
        message - str, message to be given to user
        !! FatalError
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
        message - str, message to be given to user
        !! NormalError
        """
        s = '\nError: '+str(message)
        s += '\n\t' + tools.lastError()
        s += '\nTraceBack: \n' + tools.lastErrorTrace() + '\n'

        self.log.add(s)
        raise NormalError


    def warning( self, message, error=1, trace=0 ):
        """
        Issue a warning. No exception is raised.
        message - str, message to be given to user
        error   - 1|0, report Exception with line [1]
        trace   - 1|0, report full back trace to exception [0]
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
