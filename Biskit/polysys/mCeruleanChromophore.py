from Chromophore import Chromophore
from Biskit.PDBModel import *
from numpy import array

#mCerulean : PDB - 2Q57
mCerulanChromAList = range(478,501)


class mCeruleanChromophore ( Chromophore) :
	
	def __init__ (self):
		Chromophore.__init__(self)
		self.indexes =  mCerulanChromAList
	
	
	def getTransMoment ( self , mCerulean):
		Chromophore.getTransMoment(self)
		if (self.indexes == None):
			print "[PolySys.Chromophore ERROR]: You must define atom numbers before asking for the transistion moment"
			return 0,0,0
		else:
			#Calculate  vector
			# in 2Q57, the transition dipole is defined roughly from atom 507 to 520
			points = mCerulean.take([507,520]).xyz
			return points[1] - points[0]
		

protein = PDBModel( '1huy.pdb' )

#~ crom = Chromophore (  )
#~ crom.getVector()
#~ crom.setAList ( mCitrineChromAList)
#~ crom.getVector()

crom2 = mCeruleanChromophore (  )
#~ crom2.getVector()
#~ crom2.setAList ( mCitrineChromAList)
print crom2.getTransMoment(protein)