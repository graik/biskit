from FRET import FRET
from FRETProtein import FRETProtein

class FRETPair (object) :
	
	def __init__(self, qy=0.,lt=0.):
		self.quantumYield = qy
		self.lifetime = lt
		self._donor = None
		self._acceptor = None
		self.overlap = 0

	def __getDonor(self):
		assert isinstance( _donor,FRETProtein),\
			'Donor is not well defined'
		return self._donor

	def __setDonor(self,donor):
		assert isinstance( _donor,FRETProtein),\
			'Donor must be a FRET protein'
		self._donor = donor

	donor = property(fget=__getDonor,fset=__setDonor)
	
	def __getAcceptor(self):
		assert isinstance( _acceptor,FRETProtein),\
			'Acceptor is not well defined'
		return self._acceptor

	def __setAcceptor(self,acceptor):
		assert isinstance( _acceptor,FRETProtein),\
			'Acceptor must be a FRET protein'
		self._acceptor = acceptor

	acceptor = property(fget=__getAcceptor,fset=__setAcceptor)

	def getChromoDistance(self):
		assert( donor != None and acceptor != None) 
		return 0

	def getEnergyTransferRate(self):
		
		fretCalculator= FRET ()
	
	def _loadFromDB(self,name):
		f = open ('./fret_prots_db/single_parameters.db',"r")
		lineas = f.readlines()
		f.close()
		line = -1
	