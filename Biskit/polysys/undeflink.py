from block import BlockEntity
from wormlikechain import 

class UndefLink (BlockEntity):
	"""
	Class for defining and creating peptide links with some characteristics in terms of force and
	length.
	"""
	
	def __init__(self, length=0., distance = 0.):
		"""
		Instantiation function.
		
		@param length: Contour length of the link.
		@type length: float
		@param distance: Start to end distance.
		@type distance: float
		"""
		
		BlockEntity.__init__(self,"Undef Link")
		self.length = length
		self.distance = distance
		self.polymodel = WormLikeChainModel
	
	def defSimModel(self , model):
		assert not isinstance(model,PolymerModel),"First argument may be a PolymerModel class."
		pass
	
	def __str__(self):
		return BlockEntity.__str__(self)
		
		
		
		
