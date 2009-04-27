#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
## last $Author: graik $
## $Date: 2009-03-08 19:10:27 +0100 (Sun, 08 Mar 2009) $
## $Revision: 658 $

# MAIN AUTHOR: Victor Gil Sepulveda
# STATUS: Work in progress
# TODO: Test cases

from math import cos
from emath import vectorangle,norm

"""
The FRET Class encapsulates all the needed operations for FRET calculation. 
Methods of this class are though to be used independently or in a pipeline fashion (for instance calling first calcK2,
then calcR0 and then energyTransferRate). 
"""
class FRET:
	
	"""
	Refraction coefficient of water at room temperature.
	"""
	WATER_REFRACTION_COEF = 1.3342
	
	def __init__(self, overlap, qyD,kappa2 = 2./3, tauA = 0., tauD =1., tauDA = 0., refr_index=None):
		"""
		It initializes the FRET instance with some mandatory parameters for analytic results and then
		some optional experimental values.
		
		@param overlap: Integral overlap (nm^6).
		@type: float
		@param qyD: Donor quantum yield.
		@type: float
		@param kappa2: Orientation factor for the chromophores.
		@type kappa2: float
		@param tauA: Acceptor lifetime in absence of donor.
		@type tauA: float
		@param tauD: Donor lifetime in absence of acceptor.
		@type tauD: float
		@param tauDA: Donor lifetime in presence of the acceptor.
		@type tauDA: float
		@param refr_index: Refraction index of the medium
		@type refr_index: float
		"""
		self.tauA = tauA
		self.tauD = tauD
		self.tauDA = tauDA
		self.overlap = overlap
		self.qyD = qyD
		self.k2 = kappa2
		
		if refr_index == None:
			self.refr_index = self.WATER_REFRACTION_COEF
		else:
			self.refr_index = refr_index
	
	def energyTransferRate(self, r):
		"""
		Calculates the energy transfer rate for a pair of chromophores.
		
		@param r: Distance (A) between those molecules. Note that r is a scalar, not a vector.
		@type r: float
		
		@return: Energy transfer rate.
		@rtype: float
		"""
		self.R0 = self.calcR0()	
		
		return ( 1.0/self.tauD) * (self.R0/r)**6

	def energyTransferEfficiency(self , r, kappa2 = None):
		"""
		Calculates the energy transfer efficiency for a pair of chromophores.
	
		@param r: Distance (A) between those molecules. Note that r is a scalar, not a vector.
		@type r: float
		@param kappa2: Orientation factor for the chromophores. If this value is provided it overrides
				the last calculated one. 
		@type kappa2: float
		
		@return: Energy transfer efficiency.
		@rtype: float
		"""
		self.k2 = self.k2 or kappa2

		R0 = self.calcR0()		
		
		return R0**6 / (R0**6 +r**6)
	
	
	def expEnergyTransferEfficiency(self, tauD =1., tauDA = 0.):
		"""
		
		@param r: Distance (A) between those molecules. Note that r is a scalar, not a vector.
		@type r: float
		@param tauD: Donor lifetime in absence of acceptor.
		@type tauD: float
		@param tauDA: Donor lifetime in presence of the acceptor.
		@type tauDA: float
		
		"""
		return  1 - (self.tauDA/self.tauD)
	
	
	def calcK2( self, donor, acceptor, distance):
		"""
		Calculates the orientation factor of two chromophores. It overwrites the last value.
		
		@param r: Distance (A) between those molecules. Note that r is a scalar, not a vector.
		@type r: float
		@param donor: Transition dipole moment of the donor chromophore.
		@type donor: float iterable
		@param acceptor: Transition dipole moment of the acceptor chromophore.
		@type acceptor: float iterable
		@param distance: Distance from the start of one chromophore dipole to the start of the dipole moment
				of the other chromophore.
		@type distance: float tuple
		
		@return: Orientation factor.
		@rtype: float
		"""
		#donor-acceptor angle
		thetaT = vectorangle(donor, acceptor)
		
		#donor-distance angle, distance vector is along acceptor y axis (0,y,0)
		thetaD = vectorangle(donor, distance)

		#acceptor-distance angle
		thetaA = vectorangle(acceptor,distance)
		
		self.k2 =  ( cos(thetaT) - 3.*cos(thetaD)*cos(thetaA))**2.
		
		return self.k2
		
	def calcR0 (self ): 
		"""
		Calculates R0 (Forster distance) for the pair of chromophores defined in the FRET instance.
		
		@return: Forster distance (A).
		@rtype: float
		
		"""
		n4 = self.refr_index **4
		
		self.R0 = (8.79e-5*self.k2*self.qyD*self.overlap/n4) ** (1./6) 
		
		return self.R0
		
	