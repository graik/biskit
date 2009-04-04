

import numpy as N
from Biskit.PDBModel import *
from  FRETProtein import *


class Chromophore(BlockEntity):
	
	def __init__(self,name):
		BlockEntity.__init__(self,name)
		self.indexes = None
		self.name = name
	
	def getModel (self,p,file = None) :
		"""
		Gets chromophore model structure from atoms
		"""
		if (self.indexes == None):
			print "[Polysys.Chromophore ERROR] No atomic info available."
		else:
			#get the atom info from model
			self.atoms_pos = p.indicesFrom('serial_number', a)
			mask = N.zeros( len(p) )
			N.put( mask, self.atoms_pos, 1 )
			self.atoms = p.compress (mask)
			if(file != None):
				self.atoms.writePdb(file)
			return self.atoms
	
		
	def getTransMoment(self):
		pass
	
	def __str__(self):
		return BlockEntity.__str__(self)

