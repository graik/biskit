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
## last $Author$
## $Date: 2009-03-08 19:10:27 +0100 (Sun, 08 Mar 2009) $
## $Revision$

# MAIN AUTHOR: Victor Gil Sepulveda
# STATUS: Work in progress
# TODO: Test cases

from math import cos
from emath import vectorangle,norm


class FRET:
    """
    The FRET Class encapsulates all the needed operations for FRET calculation. 
    Methods of this class are though to be used independently or in a pipeline fashion (for instance calling first calcK2,
    then calcR0 and then energyTransferRate). 
    """
    
    WATER_REFRACTION_COEF = 1.3342
    """
    Refraction coefficient of water at room temperature.
    """
    
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
        
        @return: Energy transfer rate (ns^-1).
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
        if kappa2 != None:
            R0 = self.calcR0(kappa2)
        else:
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
        
    def calcR0 (self ,kappa2 = None): 
        """
        Calculates R0 (Forster distance) for the pair of chromophores defined in the FRET instance.
        
        @return: Forster distance (A).
        @rtype: float
        
        """
        n4 = self.refr_index **4
        
        if kappa2 == None:
            k2 = kappa2
        else:
            k2 = self.k2
        
        self.R0 = (8.79e-5*k2*self.qyD*self.overlap/n4) ** (1./6) 
        
        return self.R0

		
##############
## Test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
	""" Test cases for FRET """
	
	def prepare(self):
		self.fret = None

	def cleanUp( self ):
		pass
	
	
	def test_Efficiency(self):
		"""FRET efficiency test (mCerulean- mCitrine)"""
		
		self.fret = FRET( 3.5, 959732990869504)
		self.assertAlmostEqual( 0.799022 , self.fret.energyTransferEfficiency(50.),3)
		
		#efficiency at Forster distance is 0.5
		self.assertAlmostEqual( 0.5 , self.fret.energyTransferEfficiency(self.fret.calcR0()),2)
		
		#k2 may be overwritten
		self.fret.energyTransferEfficiency(50.,3.)
		self.assertEqual(self.fret.k2,3.)
		
	def test_K2(self):
		"""K2 calculation test (mCerulean- mCitrine)"""
		
		self.fret = FRET( 3.5, 959732990869504)
		
		#random test
		self.assertAlmostEqual(self.fret.calcK2((0.1,0.3,0.3),(4.,1.,0),(10.,0.,0.)) ,0.077399380805,5)
		
		#parallel test
		self.assertAlmostEqual(self.fret.calcK2((1.,0.,0.),(1.,0.,0.),(1.,0.,0.)),4.,5)
		
		#perpendicular test
		self.assertAlmostEqual(self.fret.calcK2((1.,0.,0.),(0.,1.,0.),(1.,0.,0.)),0.)
	
	def test_R0(self):
		"""R0 calculation test (mCerulean- mCitrine)"""
		
		self.fret = FRET( 3.5, 959732990869504)
		# k2 is 2/3
		self.assertEqual(self.fret.k2,2./3.)
		
		self.assertAlmostEqual(self.fret.calcR0(),62.9320757163)
		
		# k2 must be overwritten
		self.fret.calcK2((1.,0.,0.),(1.,0.,0.),(1.,0.,0.))
		
		# so R0 may be different
		self.assertAlmostEqual(self.fret.calcR0(),84.8328253872)
		
		
if __name__ == '__main__':

    BT.localTest()