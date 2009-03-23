
from math import *


class FRET:
	
	def __init__(self, tauA, tauD, tauDA, overlap, refr_index, qyD):
		self.tauA = tauA
		self.tauD = tauD
		self.tauDA = tauDA
		self.overlap = overlap
		self.refr_index = refr_index
		self.qyD = qyD
		self.WaterRefractionRoom = 1.3342
	
	def energyTransferRate(self, thetaA,thetaD,thetaT,r):
		self.kappa2 = ( cos(self.thetaT) - 3.*cos(self.thetaD)*cos(self.thetaA))**2. 
		self.efficiency = 1. - (self.tauDA/self.tauD)
		self.R0 = 2.11e-2 *(self.kappa2*self.overlap*(self.refr_index**(-4))*self.qyD)**(1./6.)
		return ( 1.0/self.tauD) * (self.R0/r)**6

