from FRET import FRET
from Protein import FRETProtein
from FRETUtils import dbPre


class FRETPair (object) :
	
	def __init__(self, donor,acceptor):
		self._donor = donor
		self._acceptor = acceptor
		self._loadFromDB(self._donor,self._acceptor)
		

	def __getDonor(self):
		assert isinstance( _donor,FRETProtein),\
			'Donor is not well defined'
		return self._donor

	def __setDonor(self,donor):
		assert isinstance( _donor,FRETProtein),\
			'Donor must be a FRET protein'
		self._donor = donor
		self._loadFromDB(self._donor,self._acceptor)

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