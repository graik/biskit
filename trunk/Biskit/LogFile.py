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
Simple log file.
"""

import Biskit.tools as T
import sys

class LogFile:
    """
    Simple log file that can be passed between objects.
    It is flushed after each writing and should hence always be
    up2date.
    """

    def __init__(self, fname, mode='w'):
        """
        @param fname: name of log file
        @type  fname: str
        @param mode: mode (default: w)
        @type  mode: str
        """
        self.fname = T.absfile( fname )
        self.mode  = mode
        self._f  = None


    def f( self ):
        """
        Open file only when needed for first time.
        
        @return: open file handle
        @rtype: object
        """
        if self._f is None:
            self._f = open( self.fname, self.mode )

        return self._f


    def add(self, s):
        """
        Add new line to logfile and flush

        @param s: line
        @type  s: str        
        """
        self.f().writelines(s+"\n")
        self.f().flush()


    def add_nobreak(self, s):
        """
        Add new line without creating a linebreak to logfile and flush

        @param s: line
        @type  s: str
        """
        self.f().writelines(s)
        self.f().flush()


    def __del__(self):
        """
        Close file handle.
        """
        if self._f is not None:
            self.f().close()
            self._f = None


class ErrLog( LogFile ):
    """
    Print to stderr instead.
    """

    def __init__(self):
        self.fname = None
        self.mode  = None
        self._f   = sys.stderr

    def __del__(self):
        pass


class StdLog( LogFile ):
    """
    Print to std out.
    """
    def __init__(self):
        self.fname = None
        self.mode  = None
        self._f   = sys.stdout

    def __del__(self):
        pass



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: something
        @rtype:  float
        """
        import tempfile

        f_out = tempfile.mktemp( '_test_LogFile' )
        
        l = LogFile( f_out, mode='w')

        l.add('A line')

        l.add_nobreak('A nonbreakling line')

        l.add('... more text.')

        if local:
            print 'A log file was written to %s'%f_out
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


