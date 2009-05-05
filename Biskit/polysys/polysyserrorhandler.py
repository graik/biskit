from Biskit.ErrorHandler import ErrorHandler
from Biskit.Errors import BiskitError


class PolySysError( BiskitError ):
    pass
    
class PolySysFatal( BiskitError ):
    pass
    
    
class PolySysErrorHandler :
	"""
	Small class for error handling of the PolySys package.
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
			raise PolySysFatal(message)


	def error( self, message ):
		"""
		Raises an error exception
		"""
		self.lastError =message
		if self.verbose :
			raise PolySysError(message)
		

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
		"""PolySys exception catching  test """
		
		p = PolySysErrorHandler(verbose = True)

		try:
			p.error("lol!")
		except PolySysError:
			pass
		else:
			self.assertEqual (1,2)
			
		try:
			p.fatal("lol!")
		except PolySysFatal:
			pass
		else:
			self.assertEqual (1,2)
		
		p.verbose = False
		p.warning("lol!")
		self.assertEqual ( p.lastWarning, "lol!")
		
if __name__ == '__main__':
	BT.localTest()
