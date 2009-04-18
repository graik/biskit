from math import cos


class FRET:
	
	def __init__(self, tauA, tauD, tauDA, overlap, qyD,epsilon,refr_index=None):
		self.tauA = tauA
		self.tauD = tauD
		self.tauDA = tauDA
		self.overlap = overlap
		self.epsilon = epsilon
		self.qyD = qyD
		self.WaterRefractionRoom = 1.3342
		if refr_index == None:
			self.refr_index = self.WaterRefractionRoom
		else:
			self.refr_index = refr_index
	
	def energyTransferRate(self, thetaA,thetaD,thetaT,r):
		self.kappa2 = ( cos(thetaT) - 3.*cos(thetaD)*cos(thetaA))**2. 
		self.efficiency = 1. - (self.tauDA/self.tauD)
		self.R0 = 2.11e-2 *(self.kappa2*self.overlap*self.refr_index**(-4)*self.qyD)**(1./6.)
		return ( 1.0/self.tauD) * (self.R0/r)**6

	def energyTransferEfficiency(self, thetaA,thetaD,thetaT,r):
		self.kappa2 = ( cos(thetaT) - 3.*cos(thetaD)*cos(thetaA))**2.
		self.R0 = 2.11e17 *(self.kappa2*self.overlap*self.epsilon*self.refr_index**(-4)*self.qyD)**(1./6.)		
		return 1/(1+(r/self.R0)**6)

	