
from Chromophore import Chromophore
from Biskit.PDBModel import *
from numpy import array

#mCitrine : PDB - 1HUY
mCitrineChromA = [507,520]


class mCitrineChromophore ( Chromophore) :
	
	def __init__ (self):
		Chromophore.__init__(self,"mCitrineChromophore")
		self.indexes =  mCitrineChromA
	
	
	def getTransMoment ( self , mCitrine):
		Chromophore.getTransMoment(self)
		if (self.indexes == None):
			print "[PolySys.Chromophore ERROR]: You must define atom numbers before asking for the transistion moment"
			return 0,0,0
		else:
			#Calculate  vector
			# in H1UY, the transition dipole is defined roughly from atom 507 to 520
			points = mCitrine.take([507,520]).xyz
			return points[1] - points[0]
	

#~ protein = PDBModel( '1huy.pdb' )

#~ crom = Chromophore (  )
#~ crom.getVector()
#~ crom.setAList ( mCitrineChromAList)
#~ crom.getVector()

#~ crom2 = mCitrineChromophore (  )
#~ crom2.getVector()
#~ crom2.setAList ( mCitrineChromAList)
#~ print crom2.getTransMoment(protein)
