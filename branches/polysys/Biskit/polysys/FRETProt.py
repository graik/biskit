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
##
## last $Author$
## last $Date: 2009-03-09 09:10:46 +0100 (Mon, 09 Mar 2009) $
## $Revision$

# MAIN AUTHOR: Victor Gil Sepulveda
# STATUS: Work in Progress
# TODO: Test cases 

from FRETUtils import dbPre
from Chromophore import Chromophore
from Protein import Protein


"""
Class for FRET Entity definition. It loads its parameters from a database of FRET proteins.
Note that a FRET can be used for proteic or not proteic FRET entities definition.
"""
class FRETEntity (Protein):
	"""
	Creates an instance of a FRET Protein ( or just a FRET entity if id parameter is defaulted).
	Please see FRETProtein::_loadFromDB for DB data storage information .
	Automatically creates its own Chromophore instance.
	
	@param name: ID  for the database search. If defaulted, an empty FRETProtein is created which can be filled with 
				any kind of chromophore or data by hand.
	@type name: string
	"""
	def __init__ (self,name="FRETProtein"):
		
		self.name = name
		
		if self.name != "FRETProtein":
			
			self._loadFromDB(name)
			
			self.chromophore = Chromophore(self.name+'_Cromophore',self._source,self)
			
			Protein.__init__(self,name,self._source+'.pdb')
			
			
			if not self.chromophore.defined:
				print "[WARNING FRETProtein __init__] Chromophore not defined. Defaulting." 
				self.name = "mCitrine"
				self.chromophore = Chromophore(self.name+'_Cromophore_Default',self._source,self)
				self.name = name
	"""
	@return: String representation of a FRET protein.
	@rtype: string
	"""
	def __str__(self):
		desc = Protein.__str__(self)
		desc += "\n["+str(self.name)+", lifetime: "+str(self.lifetime)+", abso. wavelenth: "+str(self.abswl)+", emis. wavelenth: "+str(self.emwl)+", quantum yield: "+str(self.quantumYield)+", PDB source file: "+self._source+"]"
		return desc
	
	"""
	Method for Assembly insertion. See BlockEntity::onInsertion for further information.
	This method has no interest for the coder if FRET module is used as an standalone module.
	
	@param myassembly: Assembly where the FRETProtein instance is going to be added.
	@type myassembly: Assembly
	"""
	def onInsertion(self,myassembly=None):
		if self.chromophore.defined:
			chromo = Block(self.chromophore.name,self.chromophore,self)
			chromo.addInterval(self.chromophore.atomrange[0],self.chromophore.atomrange[1])
			myassembly.addBlock(chromo)
	
	"""
	Loads FRET protein parameters from DB.
	
	A line of the database file ( 'fret_prots_db/single_parameters.db' ) contains a total of 7 entries separated by tabs plus an
	eigth comments entry.
	
	One example line is:	
	
	mCitrine2	3.6	516	529	90e3	0.62	1HUY	mCitrine V.2, changes in epsilon (latest)
	
	Where:
	- mCitrine2
	- 3.6 is the lifetime of the protein in absence of an acceptor.
	- 516 and 529 are the maximum absortion and emission wavelengths (nm).
	- 90e3 is the extinction coefficient at maximum absorvance (M^-1 cm^-1).
	- 0.62 is he quantum yield.
	- 1HUY is de PDB id of the structure model.
	
	And the rest of the line is for comments.
	
	Sections containing 'X' will be defaulted. An Error warning will raise if the missing parameter is mandatory.
	
	@param name: is the ID of the protein. As the method of search returns the last hit, if one searches for 'mCitrine' 
				in the DB and it contains 'mCitrine1' and 'mCitrine2', it will read the last.
	@type name: string
	
	@return: @True if there were no problems during loading.
	@rtype: bool
	"""
	def _loadFromDB(self,name):
		f = open ('./fret_prots_db/single_parameters.db',"r")
		lineas = f.readlines()
		f.close()
		line = -1
		
		if len(lineas)<=1:
			print "[WARNING FRETProtein _loadFromDB (__init__)] "+name+" is not an available name in database (empty or bad format database?)." 
			return False
		else:
			for i in range(len(lineas)):
				if name in lineas[i]:
					line = i 
		if line == -1 :
			print "[WARNING FRETProtein _loadFromDB (__init__)] "+name+" is not an available name in database (name not found)." 
			return False
		
		parameters = lineas[line].split('\t',7)
		
		self.lifetime = dbPre(parameters[1],'float')
		self.abswl = dbPre(parameters[2],'float',-1)
		self.emwl = dbPre(parameters[3],'float',-1)
		self.epsilon = dbPre(parameters[4],'float')
		self.quantumYield = dbPre(parameters[5],'float')
		
		self._source = parameters[6]
		if self._source =='X' or self.abswl == -1 or self.emwl == -1:
			print "[ERROR FRETProtein _loadFromDB (__init__)] some mandatory parameters are not defined (X) in DB." 
			return False
		
		if len(parameters)>6:
			self.notes = parameters[7]
		
		return True
		
##############
## Test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
	""" """
	
	def prepare(self):
		self.fret = FRET( 3.5, 959732990869504)
		self.efficiency = 0

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