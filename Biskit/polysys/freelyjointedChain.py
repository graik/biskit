
import math
from polymermodel import PolymerModel


class FreelyJointedChain (PolymerModel):
	
	def __init__ ( self ):
		PolymerModel.__init__(self)
		self.n = 0. 
		self.b = 0.
		self.Lc = 0.
	
	def getContourLengthFromX (self) :
		"""
		Self explanatory name. It returns polymer length for a given start to end separation.
	
		- Input:
			n = number of segments
			x = start to end distance (A)
			T = temperature, defaulted to room temperature (K)
	
		- Output: 
			Contour Length (A)
		"""
		
		
		PolymerModel.getContourLengthFromX(self)
		
		k = 0.138 #pN A
		
		n=self.n;T=self.T;x=self.x;F=self.F
		
		return n*math.sqrt(3*k*T*x/(n*F));
	
	
	
	def getForceFromX( self):
		"""
		Worm Like Chain basic equation.
		
		- Input:
			p = persistence length (A)
			Lc = Contour Length (A)
			x = start to end distance (A)
			T = temperature, defaulted to rommtemperature (K)
		
		- Output: 
			Force (pN)
		"""
		
		PolymerModel.getForceFromX(self)
		
		k = 0.138 #pN A		
		
		n=self.n;T=self.T;x=self.x;b = self.b
		
		Ec = (3.*k * T )/(float(n)*float(b)*float(b))
		
		return Ec*x

	
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
	""" Test cases for Polyfret"""
	
	def prepare(self):
		pass

	def cleanUp( self ):
		pass
	
	def test_FJC(self):
		"""FJC test cases"""
		n = 1400.;  b = 4. ;Lc = n*b; x = 0.75*Lc 

		fjc = FreelyJointedChain()
		fjc.x = x
		fjc.n = n
		fjc.b = b
		fjc.Lc = Lc
		fjc.T = 298.
		
		F=  fjc.getForceFromX()
		
		fjc.F = F
		
		Lc2=  fjc.getContourLengthFromX()

		if self.local: 
			print F, Lc, Lc2 
		
		self.assertEqual(Lc2 - Lc <1.,True)
		#~ for i in range(0, 4000):
			#~ print i, fjc.getForceFromX(n,b,i)
		
		
		
if __name__ == '__main__':
    BT.localTest()	


	