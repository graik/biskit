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
## last $Author: graik $
## last $Date: 2009-03-09 09:10:46 +0100 (Mon, 09 Mar 2009) $
## $Revision: 660 $

# MAIN AUTHOR: Victor Gil Sepulveda
# STATUS: Work in Progress
# TODO: Test cases 

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