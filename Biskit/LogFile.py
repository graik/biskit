##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

import Biskit.tools as T
import sys

class LogFile:
    """
    Simple log file that can be passed between objects.
    It is flushed after each writing and should hence always be
    up2date.
    """

    def __init__(self, fname, mode='w'):
        self.fname = T.absfile( fname )
        self.mode  = mode
        self._f  = None

    def f( self ):
        """
        Open file only when needed for first time.
        -> file object
        """
        if self._f is None:
            self._f = open( self.fname, self.mode )

        return self._f

    def add(self, s):
        self.f().writelines(s+"\n")
        self.f().flush()

    def add_nobreak(self, s):
        self.f().writelines(s)
        self.f().flush()

    def __del__(self):
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
