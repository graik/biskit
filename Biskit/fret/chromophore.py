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

"""
comments from Raik:
* converted tabs to space

* __str__: won't work any longer, will it?

* conceptually: The current way to anchor the chromophore into the protein 
structure is not very flexible. The exact atom numbers can easily change 
between different structures of the same molecule (or even different 
versions of the same PDB). Something based on residue and atom names 
would be probably more stable.

Mhm, and I still think it would be more logical to have all the fluorescence
parameters in a Chromophore entry and then have a way to anchor this 
chromophore into different structures (the job of FretEntity). 
So, Chromophore would provide all the parameters, FretEntity would map the 
Chromophore into a structure, and FretPair performs all the calculations.
"""

from numpy import zeros,array
from fretutils import dbPre
from freterrorhandler import FRETErrorHandler
import Biskit.tools as T


class Chromophore:
    """
	Class for chromophore data storage, mainly the transition dipole moment references. It tries to load is parameters from DB.
	"""

    DEFAULT_DB = T.dataRoot()+'/fret/fret_prots_db/chromophores.db' # default value for database file

    def __init__(self,name,source,database="", verbose = True):
        """
		Creates an instance of a given chromophore. Loads its parameters from the DB. Please see Chromophore::_loadFromDB .

		@param name : Name of the FRETEntity it belongs in databases.
		@type name: string
		@param source : Source of the FRETEntity it belongs in databases.
		@type source: string
		@param database : Database file.
		@type database: string
		@param verbose: If True shows screen messages and throws exceptions.
		@type verbose: bool
		"""

        self.name = name
        self.apppoint = 0
        self.vector = array([0,0,0])
        self.points = []
        self.modifiers = []
        self.atomrange = [0,0]

        self.ehandler = FRETErrorHandler(verbose)

        if not self._loadFromDB(source,database):
            self.ehandler.warning( "Couldn't create the chromophore ("+self.name+","+source+")." )


    def __str__(self):
        """
		@return: String representation of a chromophre.
		@rtype: string
		"""
        return BlockEntity.__str__(self)


    def _loadFromDB(self,source,database = ""):
        """
		Loads chromophore parameters from DB.

		A line of the database file ( 'fret_prots_db/chromophores.db' ) contains a total of 4 entries separated by tabs plus 2*n 
		for transition dipole moment definition.

		One example line is:

		2Q57	mCerulean	478,501	478		1	478,501	0.5	500,501

		Which means:
		- Model used for atom description is 2Q57 (PDB)
		- Id of the FRET Protein is mCerulean
		- The chromophore is defined in the range 478,501
		- Transition dipole moment starting point is atom 478

		And then the chromophore calculation parameters defined here as two pairs of values (a b) where a is a float 
		which multiplies b, a vector defined by the positions of two atoms.

		Here dipole moment would be a vector applied at atom 478 and defined as the sum of the vector starting in atom 
		478 to atom 501 and 0.5 times the vector starting at atom 500 and ending at atom 501.

		Sections containing 'X' will be defaulted. An Error warning will raise if the missing parameter is mandatory.

		@param source: PDB id of the model.
		@type source: string
		@param database: If defined it will load data from this file instead. 
		@type database: string
		"""

        if database == "" :
            database = self.DEFAULT_DB

        f = open (database,"r")

        lineas = f.readlines()
        f.close()

        line = -1
        if len(lineas)<=1 :
            self.ehandler.warning( self.name+" is not an available name in database (empty or bad format database?)." )
            return False
        else:
            for i in range(len(lineas)):
                if source in lineas[i] and self.name in lineas[i]:
                    line = i 

        if line == -1 :
            self.ehandler.error(  "There is no coincidence for ("+self.name+","+source+") (pair of values doesn't exist)." )
            return False

        parameters = lineas[line].split()

        if len(parameters)<2:
            return False

        vectors = parameters[4:]

        self.apppoint = dbPre(parameters[3],'int',ehandler = self.ehandler) -1 # as atoms start at 0

        self.atomrange = dbPre( parameters[2],'int_range',[0,0],ehandler = self.ehandler)
        self.atomrange[0] -= 1 # as atoms start at 0
        self.atomrange[1] -= 1 # as atoms start at 0

        if 'X' in parameters[2]:
            self.ehandler.error(  "Some mandatory parameters are not defined (X) in DB, Chromophore will not be defined." )
            return False

        if len(vectors) %2 != 0 or  len(vectors) ==0:
            self.ehandler.error( "Error in vector description in "+source+" record" )
            return False

        self.modifiers = [0.]*(len(vectors)/2)
        self.points = [[0,0]]*(len(vectors)/2)

        for i in range(len(vectors)):
            j = i/2
            if i%2 == 0:
                self.modifiers[j] = dbPre(vectors[i],'float',ehandler = self.ehandler)
            else:
                self.points[j] = dbPre (vectors[i],'int_range',ehandler = self.ehandler)
                self.points[j][0] -= 1 # as atoms start at 0
                self.points[j][1] -= 1 # as atoms start at 0

        return True

##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for Chromophore"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass

    def test_Creation(self):
        """Instantiation test """

        try:
            c=Chromophore("mCitrine","1HUY",T.testRoot() + "/fret/chromo.db",verbose = True)		
            self.assertEqual(c.ehandler.lastError,"")
            self.assertEqual(c.ehandler.lastWarning,"")
        except:
            self.assertEqual(1,2)



    def test_Loading(self):
        """Loading from database test """


        # Not available db
        c=Chromophore("mCitrine","1HUY",T.testRoot() + "/fret/lol.db",verbose = False)		
        self.assertEqual ( "Couldn't create the chromophore" in c.ehandler.lastWarning,True) 

        # Forcing errors

        # Non-existing pair
        c=Chromophore("mCitrininine","1HUY",database = T.testRoot() + "/fret/chromo.db",verbose = False)
        self.assertEqual ( "There is no coincidence for" in c.ehandler.lastError,True)

        # No mandatory
        c=Chromophore("Mandatory1","FAKE1",T.testRoot() + "/fret/chromo.db",verbose = False)
        self.assertEqual ( "mandatory parameters are not defined" in c.ehandler.lastError,True)

        # No vector
        c=Chromophore("Vector","FAKE2",T.testRoot() + "/fret/chromo.db",verbose = False)
        self.assertEqual ( "Error in vector description" in c.ehandler.lastError,True)


if __name__ == '__main__':

    BT.localTest()