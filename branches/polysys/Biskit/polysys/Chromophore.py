

import numpy as N
from Biskit.PDBModel import *
from  FRETProtein import *


class Chromophore(BlockEntity):
	
	def __init__(self,name,source,myFRETProtein):
		BlockEntity.__init__(self,name,myFRETProtein)
		self.indexes = None
		_loadFromDB(source)
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
		points = self.father.take(indexes).xyz
		return points[1] - points[0]
	
	def __str__(self):
		return BlockEntity.__str__(self)

	def _loadFromDB(self,source):
		f = open ('./fret_prots_db/Chromophores.db',"r")
		lineas = f.readlines()
		f.close()
		line = -1
		
		if len(lineas)<=1:
			print "[WARNING Chromophore _loadFromDB (__init__)] "+name+" is not an available name in database (empty or bad format database?)." 
			return False
		else:
			for i in range(len(lineas)):
				if source in lineas[i]:
					line = i 
		if line == -1 :
			print "[WARNING Chromophore _loadFromDB (__init__)] "+name+" is not an available name in database (name not found)." 
			return False
		
		parameters = lineas[line].split('\t')
		
		if len(parameters)<2:
			return False
			
			
		vectors = parameters[1:] 
		
		print vectors
		
		return True