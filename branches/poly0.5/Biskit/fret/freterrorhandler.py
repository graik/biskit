"""
Comment from Raik:

This really duplicates Biskit.ErrorHandler without reporting Error traces from
exceptions. -- better use the modified ErrorHandler (see last checkin).

Testing: there is a unittest.TestCase.assertRaises method that you can use
instead of the assertEqual( 1, 2 ).
"""

from Biskit.ErrorHandler import ErrorHandler
from Biskit.Errors import BiskitError

class FRETError( BiskitError ):
    pass

class FRETFatal( BiskitError ):
    pass


class FRETErrorHandler :
    """
    Small class for error handling of the FRET package.
    """

    def __init__(self,verbose = False,log= None):
        """
        Creates an instance and defines variables. If "verbose" is then set to true, it will print in 
        screen the errors. It allvays stores the last error in "lastError" .
        """
        self.e = ErrorHandler(log)
        self.lastError =  ""
        self.lastWarning = ""
        self.verbose = verbose


    def fatal( self, message ):
        """
        Raises an error exception.
	"""
        self.lastError =message
        if self.verbose :
            raise FRETFatal(message)


    def error( self, message ):
        """
	Raises an error exception
	"""
        self.lastError =message
        if self.verbose :
            raise FRETError(message)


    def warning( self, message, error=1, trace=0 ):
        """
	See ErrorHandler::warning.
	"""
        self.lastWarning =message
        if self.verbose :
            self.e.warning(self.lastWarning)


##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for Error Handling"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass

    def test_Except(self):
        """FRET exception catching  test """

        f = FRETErrorHandler(verbose = True)

        try:
            f.error("lol!")
        except FRETError:
            pass
        else:
            self.assertEqual (1,2)

        try:
            f.fatal("lol!")
        except FRETFatal:
            pass
        else:
            self.assertEqual (1,2)

        f.verbose = False
        f.warning("lol!")
        self.assertEqual ( f.lastWarning, "lol!")

if __name__ == '__main__':
    BT.localTest()
