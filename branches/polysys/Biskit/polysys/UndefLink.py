
from FRETProtein import BlockEntity

class UndefLink (BlockEntity):
	def __init__(self, length=0., distance = 0.):
		BlockEntity.__init__(self,"Undef Link")
		self.length = length
		self.distance = distance
		self.polymodel = None
	
	def __str__(self):
		return BlockEntity.__str__(self)