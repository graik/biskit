

import numpy as N
from Biskit.PDBModel import *
from  FRETProtein import *
from FRETUtils import dbPre

# This class may not be accessible to the end user
class Chromophore(BlockEntity):
	
	def __init__(self,name,source,myFRETProtein):
		BlockEntity.__init__(self,name,myFRETProtein)
		
		
		self.name = name
		self.apppoint = 0
		self.vector = array([0,0,0])
		self.points = []
		self.modifiers = []
		self.atomrange = [0,0]
		
		if not self._loadFromDB(source):
			print "[WARNING Chromophore __init__] Couldn't create the chromophore." 
			self.defined = False
		else:
			self.defined = True
			
		
	
	def getModel (self,file = None) :
		"""
		Gets chromophore model structure from atoms
		"""
		if (self.indexes == None):
			print "[ERROR Chromophore getModel] No atomic info available."
		else:
			#get the atom info from model
			p = self.father
			self.atoms_pos = p.indicesFrom('serial_number', self.atomrange)
			mask = N.zeros( len(p) )
			N.put( mask, self.atoms_pos, 1 )
			self.atoms = p.compress (mask)
			if(file != None):
				self.atoms.writePdb(file)
			return self.atoms
	
		
	def getTransMoment(self):
		
		if not self.defined:
			print "[ERROR FRETProtein getTransMoment] Transition moment is not defined."
			return [0,0,0]
		
		k = 0
		print self.points
		print self.modifiers
		for i in self.points:
			print "*"
			_points = self.father.take(i).xyz
			p = (_points[1] - _points[0])*self.modifiers[k]
			k = k+1
			self.vector = self.vector+p
		return self.vector
	
	def __str__(self):
		return BlockEntity.__str__(self)

	def _loadFromDB(self,source):
		f = open ('./fret_prots_db/chromophores.db',"r")
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
		
		parameters = lineas[line].split()
		
		if len(parameters)<2:
			return False
			
		
		vectors = parameters[4:]
		
		self.apppoint = dbPre(parameters[3],'int')
		
		
		self.atomrange = dbPre( parameters[2],'int_range',[0,0])
		if 'X' in parameters[2]:
			print "[ERROR Chromophore _loadFromDB (__init__)] some mandatory parameters are not defined (X) in DB, Chromophore will not be defined." 
			return False
		
		if len(vectors) %2 != 0:
			print "[ERROR Chromophore _loadFromDB (__init__)] Error in vector description in "+source+" record" 
			return False
		
		self.modifiers = [0.]*(len(vectors)/2)
		self.points = [[0,0]]*(len(vectors)/2)
		
		for i in range(len(vectors)):
			j = i/2
			if i%2 == 0:
				self.modifiers[j] = dbPre(vectors[i],'float')
			else:
				self.points[j] = dbPre (vectors[i],'int_range',[0,0])
		
		print self.points
		print self.modifiers
		return True