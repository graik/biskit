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
Simple log file.
"""

import biskit.tools as T
import sys

class LogFile:
    """
    Simple log file that can be passed between objects.
    It is flushed after each writing and should hence always be
    up2date.
    """

    def __init__(self, fname, mode='w'):
        """
        :param fname: name of log file
        :type  fname: str
        :param mode: mode (default: w)
        :type  mode: str
        """
        self.fname = T.absfile( fname )
        self.mode  = mode
        self._f  = None


    def f( self ):
        """
        Open file only when needed for first time.

        :return: open file handle
        :rtype: object
        """
        if self._f is None:
            self._f = open( self.fname, self.mode )

        return self._f


    def write(self, s):
        """Synonym for add_nobreak.
        :param s: line
        :type  s: str
        """
        self.add_nobreak( s )

    def writeln(self, s):
        """Synonym for add.
        :param s: line
        :type  s: str
        """
        self.add(s)

    def add(self, s):
        """
        Add new line to logfile and flush

        :param s: line
        :type  s: str        
        """
        if not type(s) is list:
            s = [ s ]

        for i in s:
            self.f().write(i)
            self.f().write('\n')

        self.f().flush()


    def add_nobreak(self, s):
        """
        Add new line without creating a linebreak to logfile and flush

        :param s: line
        :type  s: str
        """
        self.f().write(s)
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
import biskit.test as BT

class Test(BT.BiskitTest):

    def prepare(self):
        import tempfile
        self.f_out = tempfile.mktemp( '_test_LogFile' )

    def cleanUp(self):
        if T.tryRemove( self.f_out ) and self.local:
            print('log file removed.')

    def test_LogFile( self ):
        """LogFile test """

        self.l = LogFile( self.f_out, mode='w')

        self.l.writeln('1')
        self.l.write('2')
        self.l.writeln('3')

        if self.local:
            print('log file written to %s'%self.f_out)

        lines = open(self.f_out).readlines()
        self.assertEqual( lines, ['1\n','23\n'] )

    def test_StdLog(self):
        """StdLog test (only in interactive mode)"""
        if self.local:
            self.l = StdLog()
            self.l.write(' test message ... ')

    def test_ErrLog(self):
        """ErrLog test"""
        self.l = ErrLog()
        self.l.write(' test message ... ')


if __name__ == '__main__':

    BT.localTest()
