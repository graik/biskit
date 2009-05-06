
import math
from numpy import *
from emath import *
from polymermodel import PolymerModel



class WormLikeChainModel (PolymerModel):
	
	def __init__ ( self ):
		"""
		Instantiation function.
		"""
		PolymerModel.__init__(self)
		self.Lc = 0.
		self.p = 0.
	
	
	def getContourLengthFromX (self) :
		"""
		Self explanatory name. It returns polymer length for a given start to end separation.
	
		- Input:
			p = persistence length (A)
			x = start to end distance (A)
			T = temperature, defaulted to rommtemperature (K)
	
		- Output: 
			Contour Length (A)
	
		Very accurate solver for large values of x.
		"""
		
		PolymerModel.getContourLengthFromX(self)
		
		x = self.x; T = self.T;p=self.p; F= self.F
		
		k = 0.138 #K Boltzmann in pN  * A
		
		Ec = (k * T) /(4.0*p)
		
		## x^3+ax^2+bx+c =0
		a = -2.*x*(F+3.*Ec)/F
		b = x*x*(9.*Ec+F)/F
		c = -4.*x*x*x*Ec/F
		
		return cubic(a,b,c);
	
	
	def getContourLengthFromX2 (self) :
		"""
		See getContourLengthFromX for description.
		It seems that this third grade solver is more accurate for small values of x
		"""
		
		PolymerModel.getContourLengthFromX(self)
		
		x = self.x; T = self.T;p=self.p; F= self.F
		
		k = 0.138 #pN A
		Ec = (k * T) /(4.*p)
		
		## x^3+ax^2+bx+c =0
		a = -2.*x*(F+3*Ec)/F
		b = x*x*(9.*Ec+F)/F
		c = -4.*x*x*x*Ec/F
		
		## Q and R are always real
		Q= (a*a -3.*b)/9.
		R = (2.*a*a*a - 9.*a*b +27.*c) / 54.
		
		if R*R < Q*Q*Q : 
			## It has 3 real roots
			theta = math.acos(R/math.sqrt(Q*Q*Q))
			sqrQ =-2.* math.sqrt(Q)
			athi = a / 3.
			L1 = sqrQ*math.cos(theta/3)-athi
			L2 = sqrQ*math.cos((theta+6.28)/3.)-athi
			L3 = sqrQ*math.cos((theta-6.28)/3.)-athi
			return L1 , L2 , L3
			
		else:
			A = -sign(R)*(math.fabs(R) + math.sqrt(R*R-Q*Q*Q))**(1./3.)
			
			if A == 0:
				B = 0
			else:
				B = Q/A
			
			return (A+B)-(a/3.)

	
	
	def getForceFromX( self, p,Lc,x,T = 298.0):
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
		
		k = 0.138 #pN * A

		Ec = (k * T )/(4.0*p)
		
		return Ec*(((1.-(x/Lc))**-2) - 1+(4*x/Lc))

##############
## Test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
	""" Test cases for Polyfret"""
	
	def prepare(self):
		pass

	def cleanUp( self ):
		pass
	
	def test_WLC(self):
		"""WLC test cases"""
		p = 4.;Lc = 5280.; x = 0.75*5280 
		
		wlc = WormLikeChainModel()
		
		wlc.p = p;wlc.Lc=Lc;wlc.x=x;wlc.T = 298.
		
		F =  wlc.getForceFromX(p,Lc,x)
		
		wlc.F = F
		
		Lc2 =  wlc.getContourLengthFromX()
		Lc3 =  wlc.getContourLengthFromX2()

		self.assertEqual(Lc2[0] - Lc <2.,True)
		self.assertEqual(Lc3[1] - Lc <2.,True)
		
		
		
if __name__ == '__main__':
    BT.localTest()	


	
