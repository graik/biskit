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

from Biskit.fret import FRET
from Biskit.fret.fretprotein import FRETProtein
from Biskit.fret.fretutils import dbPre


class FRETPair (object) :
	"""
	
	"""
	def __init__(self, donor,acceptor):
		self._donor = donor
		self._acceptor = acceptor
		self._loadFromDB(self._donor,self._acceptor)
		
	"""
	
	"""
	def __getDonor(self):
		assert isinstance( _donor,FRETProtein),\
			'Donor is not well defined'
		return self._donor
	"""
	
	"""
	def __setDonor(self,donor):
		assert isinstance( _donor,FRETProtein),\
			'Donor must be a FRET protein'
		self._donor = donor
		self._loadFromDB(self._donor,self._acceptor)
	
	"""
	
	"""
	donor = property(fget=__getDonor,fset=__setDonor)
	
	def __getAcceptor(self):
		assert isinstance( _acceptor,FRETProtein),\
			'Acceptor is not well defined'
		return self._acceptor

	def __setAcceptor(self,acceptor):
		assert isinstance( _acceptor,FRETProtein),\
			'Acceptor must be a FRET protein'
		self._acceptor = acceptor
		self._loadFromDB(self._donor,self._acceptor)

	acceptor = property(fget=__getAcceptor,fset=__setAcceptor)

	def getChromoDistance(self):
		assert( donor != None and acceptor != None) 
		return 0

	def getEnergyTransferRate(self):
		
		fretCalculator= FRET ()
	
	def _loadFromDB(self,donor,acceptor):
		f = open ('./fret_prots_db/pair_parameters.db',"r")
		lineas = f.readlines()
		#~ print self._acceptor.name, self._donor.name
		for l in lineas:
			parameters = l.split('\t',6)
			#~ print parameters
			if self._acceptor.name == parameters[0] and self._donor.name == parameters[1] :
				f.close()
				self._lifetime = dbPre(parameters[2],"float",2.5)
				
				self._overlap = dbPre(parameters[3],"float",-1)
				
				if self._overlap == -1 :
					print "[ERROR _loadFromDB] Mandatory parameter undefined in database (overlap)."
				return
		print "[ERROR FRETPair _loadFromDB] This pair isn't defined in database."
				
		f.close()
		
d = FRETProtein('mCerulean')
a = FRETProtein('mCitrine')
p = FRETPair(d,a)