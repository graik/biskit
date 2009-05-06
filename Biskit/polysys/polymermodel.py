
class PolymerModel:
	"""
	Base class for polymer physcs.
	"""
	
	def __init__ ( self ):
		"""
		Instantiation function. It defaults some values.
		"""
		
		self.x  = 0.
		self.T = 298.
		self.F = 0.

	def getContourLengthFromX ( self ) :
		pass
	
	
	def getForceFromX( self ):
		pass
