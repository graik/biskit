## Errors for all modules
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## last $Date$

import exceptions

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
