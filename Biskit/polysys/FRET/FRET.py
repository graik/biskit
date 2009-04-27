from math import cos
from emath import vectorangle,norm

class FRET:
	
	"""
	Refraction coefficient of water at room temperature.
	"""
	WATER_REFRACTION_COEF = 1.3342
	
	def __init__(self, overlap, qyD,epsilon,kappa2 = 2./3, tauA = 0., tauD =1., tauDA = 0., refr_index=None):
		
		self.tauA = tauA
		self.tauD = tauD
		self.tauDA = tauDA
		self.overlap = overlap
		self.epsilon = epsilon
		self.qyD = qyD
		self.k2 = kappa2
		
		if refr_index == None:
			self.refr_index = self.WATER_REFRACTION_COEF
		else:
			self.refr_index = refr_index
	
	def energyTransferRate(self, thetaA,thetaD,thetaT,r):
		
		self.R0 = self.calcR0()	
		
		return ( 1.0/self.tauD) * (self.R0/r)**6

	def energyTransferEfficiency(self , r,kappa2 = None):
		
		kappa2 = self.k2 or kappa2

		R0 = self.calcR0()		
		
		return R0**6 / (R0**6 +r**6)
	
	
	def expEnergyTransferEfficiency(self, tauD =1., tauDA = 0.):
		return  1 - (self.tauDA/self.tauD)
	
	
	def calcK2( self,donor, acceptor, distance):
		#donor-acceptor angle
		thetaT = vectorangle(donor, acceptor)
		
		#donor-distance angle, distance vector is along acceptor y axis (0,y,0)
		thetaD = vectorangle(donor, distance)

		#acceptor-distance angle
		thetaA = vectorangle(acceptor,distance)
		
		self.k2 =  ( cos(thetaT) - 3.*cos(thetaD)*cos(thetaA))**2.
		
		return self.k2
		
	def calcR0 (self ): ## in A
		n4 = self.refr_index **4
		
		if self.k2 == None:
			print "[ERROR FRET calcR0] You need to calculate K2."
			
		self.R0 = (8.79e-5*self.k2*self.qyD*self.overlap/n4) ** (1./6) 
		
		return self.R0
		
	