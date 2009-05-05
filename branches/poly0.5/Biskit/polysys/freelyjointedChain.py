
import math
from PolymerModel import PolymerModel


class FreelyJointedChain (PolymerModel):
	
	def __init__ ( self ):
		PolymerModel.__init__(self)
		self.n = 0. 
		self.b = 0.
	
	
	def getContourLengthFromX (self) :
		"""
		Self explanatory name. It returns polymer length for a given start to end separation.
	
		Input:
			n = number of segments
			x = start to end distance (A)
			T = temperature, defaulted to room temperature (K)
	
		Output: 
			Contour Length (A)
		"""
		PolymerModel.getContourLengthFromX(self)
		
		k = 0.138 #pN A
				
		return n*math.sqrt(3*k*T*x/(n*F));
	
	
	
	def getForceFromX( self):
		"""
		Worm Like Chain basic equation.
		
		Input:
			p = persistence length (A)
			Lc = Contour Length (A)
			x = start to end distance (A)
			T = temperature, defaulted to rommtemperature (K)
		
		-Output: 
			Force (pN)
		"""
		
		PolymerModel.getForceFromX(self)
		
		k = 0.138 #pN A		

		
		Ec = (3.*k * T )/(float(n)*float(b)*float(b))
		
		return Ec*x

	

#~ n = 1400.;  b = 4. ;Lc = n*b; x = 0.75*Lc 

#~ fjc = FreelyJointedChain()


#~ F=  fjc.getForceFromX(n,b,x)

#~ Lc2=  fjc.getContourLengthFromX(n,x,F)


#~ #print F, Lc, Lc2 

#~ for i in range(0, 4000):
	#~ print i, fjc.getForceFromX(n,b,i)
	