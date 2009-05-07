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

"""
Comments from Raik:

* looks good!
* I converted tabs to spaces (doc strings are still a little messed up)
* docu: the get* methods are not telling which units they return
"""

from fret import FRET
from fretentity import FRETEntity
from fretutils import dbPre
from freterrorhandler import FRETErrorHandler, FRETFatal, FRETError
from emath import norm
import Biskit.tools as T

class FRETPair (object) :
    """
	A FRETPair is a container for the two FRET Proteins involved in the energy transfer. It tries to load all parameters needed for FRET
	calculations from the database.
	Mandatory paramaters for a FRETPair definition in the DB are the Overlap Integral and the IDs of the pair (their names). For more info.
	about the database see FRETPair::_loadFromDB doc.
	"""

    DEFAULT_DB = T.dataRoot()+'/fret/fret_prots_db/pair_parameters.db' # default value for database file


    def __init__(self, donor,acceptor,database="", verbose = True):
        """
		Creates an instance of a FRET Pair  

		@param donor:
		@type donor:
		@param acceptor:
		@type acceptor:
		@param database: If defined it will load data from this file instead. 
		@type database: string
		@param verbose: If True shows screen messages and throws exceptions.
		@type verbose: bool
		"""
        self._donor = donor
        self._acceptor = acceptor
        self.ehandler = FRETErrorHandler(verbose)
        self._loadFromDB(self._donor,self._acceptor,database)


    def __getDonor(self):
        """
		@return: Donor of the FRET pair.
		@rtype: FRETProtein
		"""
        assert isinstance( self._donor,FRETEntity),\
               'Donor is not well defined'
        return self._donor


    def __setDonor(self,donor):
        """
		@param donor: Donor of the FRET pair.
		@type donor: FRETProtein
		"""
        assert isinstance( self._donor,FRETEntity),\
               'Donor must be a FRET protein'
        self._donor = donor
        self._loadFromDB(self._donor,self._acceptor)


    """Set or get the donor protein of the FRET pair.	"""
    donor = property(fget=__getDonor,fset=__setDonor)



    def __getAcceptor(self):
        """
		@return: Acceptor of the FRET pair.
		@rtype: FRETProtein
		"""
        assert isinstance( self._acceptor,FRETEntity),\
               'Acceptor is not well defined'
        return self._acceptor


    def __setAcceptor(self,acceptor):
        """
		@param acceptor: Acceptor of the FRET pair.
		@type acceptor: FRETProtein
		"""
        assert isinstance( self._acceptor,FRETEntity),\
               'Acceptor must be a FRET protein'
        self._acceptor = acceptor
        self._loadFromDB(self._donor,self._acceptor)

    """Set or get the acceptor protein of the FRET pair.	"""
    acceptor = property(fget=__getAcceptor,fset=__setAcceptor)


    def getChromoDistance(self):
        """
		Calculates the distance between the starting points of the two dipole moments of the chromophores involved in 
		the energy transfer.
		"""
        if  self.donor == None or self.acceptor== None :
            self.ehandler.error ("Donor and acceptor may be defined before calculations")

        return self.acceptor.getAppPoint() - self.donor.getAppPoint()


    def getEnergyTransferRate(self):
        """
		Calculates the energy transfer rate.
		"""

        distance = self.getChromoDistance()

        fretCalculator= FRET (self.overlap, self.donor.quantumYield,tauD = self.donor.lifetime)

        k2 = fretCalculator.calcK2( self.donor.getChromoTransMoment(),self.acceptor.getChromoTransMoment(),distance)

        return fretCalculator.energyTransferRate(norm(distance))

    def getEnergyTransferEfficiency(self):
        """
		Calculates the energy transfer efficiency.
		"""
        distance = self.getChromoDistance()

        fretCalculator= FRET (self.overlap, self.donor.quantumYield)

        k2 = fretCalculator.calcK2( self.donor.getChromoTransMoment(),self.acceptor.getChromoTransMoment(),distance)

        return fretCalculator.energyTransferEfficiency(norm(distance))

    def getR0(self):
        """
		Calculates Forster's distance
		"""
        distance = self.getChromoDistance()

        fretCalculator= FRET (self.overlap, self.donor.quantumYield)

        k2 = fretCalculator.calcK2( self.donor.getChromoTransMoment(),self.acceptor.getChromoTransMoment(),distance)

        return fretCalculator.calcR0()

    def _loadFromDB(self,donor,acceptor,database= ""):
        """
		Loads the pair parameters from the database. Mandatory paramaters for a FRETPair definition in the DB ( './fret_prots_db/pair_parameters.db')
		are the Overlap Integral and the IDs of the pair (their names). 

		A line of the database file contains a total of 6 entries separated by tabs.
		One example line is:
		mCitrine	mCerulean	2.89	959732990869504	1.33	Using Short Linker (incorrect F distance?)

		Where 'mCitrine' and 'mCerulean' are the donor and acceptor. 2.89 (ns) is the lifetime of the donor in presence of the acceptor.
		9.59e14 (nm^6) is the overlap integral. 1.33 is an experimental calculation of Forster's distance and the rest of the line is 
		a comment.

		While comments are optional, all the other parameters must to be defined for each line. For unknown values one can define them
		as 'X' and they will be defaulted if needed. However, as said before, Overlap Integral and pair IDs must to be defined for the 
		FRETPair to be created.

		@param donor: Donor of the FRET pair.
		@type donor: FRETProtein
		@param acceptor: Acceptor of the FRET pair.
		@type acceptor: FRETProtein
		@param database: If defined it will load data from this file instead. 
		@type database: string
		"""

        if database == "" :
            database = self.DEFAULT_DB

        f = open (database,"r")
        lineas = f.readlines()
        f.close()

        if len(lineas)<=1:
            self.ehandler.fatal( "Empty or badly formatted database" )
            return False
        else:
            for l in lineas:
                parameters = l.split('\t',6)

                if self._acceptor.name == parameters[0] and self._donor.name == parameters[1] :

                    self.lifetime = dbPre(parameters[2],"float",2.5,self.ehandler)

                    self.overlap = dbPre(parameters[3],"float",-1,self.ehandler)

                    if self.overlap == -1 :
                        self.ehandler.error( "Mandatory parameter undefined in database (overlap).")

                    return

        self.ehandler.fatal( "This pair isn't defined in database.")



##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit import PDBModel

class Test(BT.BiskitTest):
    """ Test cases for FRETEntity"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass

    def test_Creation(self):
        """Instantiation test """
        d = FRETEntity('mCerulean',T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)
        a = FRETEntity('mCitrine',T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)

        try:
            p = FRETPair(d,a,verbose = True,database= T.testRoot() + "/fret/pair.db")
        except:
            self.assertEqual (1,2)

    def test_Loading(self):
        """Loading from database test """

        d = FRETEntity("LoadTest",T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)
        a = FRETEntity("LoadTest",T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)
        d2 = FRETEntity("mCitrine",T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)


        # Not available db
        p = FRETPair(d,a, database= T.testRoot() + "/fret/lol.db",verbose =False)
        self.assertEqual ( "Empty or badly formatted database" in p.ehandler.lastError,True) 

        # Not defined pair
        p = FRETPair(d,a, database= T.testRoot() + "/fret/pair.db",verbose =False)
        self.assertEqual ( "This pair isn't defined" in p.ehandler.lastError,True) 

        # No mandatory
        p = FRETPair(d2,a, database= T.testRoot() + "/fret/pair.db",verbose =False)
        self.assertEqual ( "Mandatory parameter undefined" in p.ehandler.lastError,True) 

    def test_FretCalcs(self):
        """FRET calculations test """
        d = FRETEntity('mCerulean',T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)
        a = FRETEntity('mCitrine',T.testRoot() + "/fret/single.db",chromo_autodef = False,verbose=True)

        try:
            p = FRETPair(d,a,verbose = True,database= T.testRoot() + "/fret/pair.db")
        except:
            self.assertEqual (1,2)

        structure = PDBModel(T.testRoot() + "/fret/2Q57.pdb")
        d.structure = structure
        structure = PDBModel(T.testRoot() + "/fret/1HUY.pdb")
        a.structure = structure

        self.assertAlmostEqual (norm(p.getChromoDistance()),51.877,3)
        self.assertAlmostEqual (p.getEnergyTransferRate(),0.268,2)
        self.assertAlmostEqual ( p.getEnergyTransferEfficiency(),0.484,2)
        self.assertAlmostEqual ( p.getR0(),51.350,2)

if __name__ == '__main__':

    BT.localTest()

