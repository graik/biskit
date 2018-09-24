## Errors for all modules
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
Errors raised by Biskit scripts.
"""

class BiskitError( Exception ):
    """
    Root for all Errors raised by Biskit scripts.
    """
    pass


class HandledError( BiskitError ):
    """
    Error raised by the ErrorHandler after an Error has been reported.
    """
    pass


class FatalError( HandledError ):
    """
    Error raised by the ErrorHandler after a fatal Error has been reported.
    """
    pass


class NormalError( HandledError ):
    """
    Error raised by the ErrorHandler after a normal Error has been reported.
    """
    pass


class XplorInputError( BiskitError ):
    """
    Errors raised while generating xplor input script
    """
    pass

##################
##  (FAKE) TESTING        
###################
import biskit.test as BT

class Test(BT.BiskitTest):
    """Error test"""

    def test_Errors( self ):
        """Errors test (empty test)"""
        pass
