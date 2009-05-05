from FRETProtein import BlockEntity

class UndefLink (BlockEntity):
	def __init__(self, length=0., distance = 0.):
		BlockEntity.__init__(self,"Undef Link")
		self.length = length
		self.distance = distance
		self.polymodel = None
	
	def defSimModel(self , model):
		if not isinstance(model,PolymerModel):
			print "[ERROR UndefLink defSimModel] First parameter may be a PolymerModel class."
			return
	
	def __str__(self):
		return BlockEntity.__str__(self)
		
		
		
		
#~ class A:
	#~ def __init__ (self):
		#~ pass
	
	
#~ class B(A):
	#~ def __init__ (self):
		#~ A.__init__(self)

#~ a = A()
#~ b = B()

#~ print isinstance(a,A),isinstance(a,B),isinstance(b,A),isinstance(b,B)