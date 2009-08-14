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

from fretutils import dbPre
from Biskit import EHandler
from chromophore import Chromophore
from numpy import array
import Biskit.tools as T

class FRETEntity:
    """
    Class for FRET Entity definition. It loads its parameters from a database of FRET proteins.
    Note that a FRET can be used for proteic or not proteic FRET entities definition.
    """

    DEFAULT_DB = T.dataRoot()+'./fret/fret_prots_db/single_parameters.db' # default value for database file

    def __init__ (self,name="FRETEntity",database="", structure = None,chromo_autodef = False):

        """
        Creates an instance of a FRET Protein ( or just a FRET entity if id parameter is defaulted).
        Please see FRETProtein::_loadFromDB for DB data storage information .
        Automatically creates its own Chromophore instance.
        If encountered an error while loading data from database, it tries to default anything it can.

        @param name: ID  for the database search. If defaulted, an empty FRETEntity is created which can be filled with 
                any kind of chromophore or data by hand.
        @type name: string
        @param database: If defined it will load data from this file instead. 
        @type database: string
        @param structure: Atomic structure of the fret entity. Must implement "take"and "xyz" (see PDBModel).
        @type structure: PDBModel (by now...)
        @param chromo_autodef: If False chromophore is not defined.
        @type chromo_autodef: bool

        """

        self.name = name
        
        self.ehandler = EHandler
        
        self.structure= structure
        
        if self.name != "FRETEntity":
            
            if  self._loadFromDB(name,database) and chromo_autodef :
                
                try:
                    self.chromophore = Chromophore(self.name,self.source)
                except:
                    self.chromophore = None
            else:
                self.chromophore = None

    def __str__(self):
        """
        @return: String representation of a FRET entity.
        @rtype: string
        """
        desc += "\n["+str(self.name)+", lifetime: "+str(self.lifetime)+", abso. wavelenth: "+str(self.abswl)+", emis. wavelenth: "+str(self.emwl)+", quantum yield: "+str(self.quantumYield)+", PDB source file: "+self.source+"]"
        return desc


    def getChromoTransMoment(self):
        """
        Calculates the transition dipole moment of the chromophore using the parameters loaded from the DB.
        For more info about how dipole moment is defined into DB (and then how is it calculated here) please 
        see Chromophore::_loadFromDB .
        """
        
        if not self.structure:
            self.ehandler.warning( "Transition moment is not defined (Structure is not defined).")
            return [0,0,0]
        
        if not self.chromophore:
            self.ehandler.warning( "Transition moment is not defined (Chromophore is not defined).")
            return [0,0,0]
        
        k = 0
        self.vector = array([0,0,0])
        
        for i in self.chromophore.points:
            _points = self.structure.take(i).xyz
            p = (_points[1] - _points[0])*self.chromophore.modifiers[k]
            k = k+1
            self.vector = self.vector+p
            
        return self.vector

    def getAppPoint(self):
        """
        Function for getting the coordinates of the application point.
        """
        return self.structure.take([self.chromophore.apppoint]).xyz[0]

    def _loadFromDB(self,name,database = ""):
        """
        Loads FRET protein parameters from DB.
        
        A line of the database file ( 'fret_prots_db/single_parameters.db' ) contains a total of 7 entries separated by tabs plus an
        eigth comments entry.
        
        One example line is:	
        
        mCitrine2	3.6	516	529	90e3	0.62	1HUY	mCitrine V.2, changes in epsilon (latest)
        
        Where:
        - mCitrine2
        - 3.6 is the lifetime of the protein in absence of an acceptor.
        - 516 and 529 are the maximum absortion and emission wavelengths (nm).
        - 90e3 is the extinction coefficient at maximum absorvance (M^-1 cm^-1).
        - 0.62 is he quantum yield.
        - 1HUY is de PDB id of the structure model.
        
        And the rest of the line is for comments.
        
        Sections containing 'X' will be defaulted. An Error warning will raise if the missing parameter is mandatory.
        
        @param name: Is the ID of the protein. As the method of search returns the last hit, if one searches for 'mCitrine' 
                    in the DB and it contains 'mCitrine1' and 'mCitrine2', it will read the last.
        @type name: string
        @param database: If defined it will load data from this file instead. 
        @type database: string
        
        @return: @True if there were no problems during loading.
        @rtype: bool
        """
        
        if database == "" :
            database = self.DEFAULT_DB
        
        f = open (database,"r")
        
        lineas = f.readlines()
        f.close()
        
        line = -1
        if len(lineas)<=1:
            self.ehandler.warning( name+" is not an available name in database (empty or badly formatted database?)." )
            return False
        else:
            for i in range(len(lineas)):
                if name in lineas[i]:
                    line = i 
        if line == -1 :
            self.ehandler.warning(  name+" is not an available name in database (name not found)." )
            return False
        
        parameters = lineas[line].split('\t',7)
        self.lifetime = dbPre(parameters[1],'float',ehandler=self.ehandler)
        self.abswl = dbPre(parameters[2],'float',-1,ehandler=self.ehandler)
        self.emwl = dbPre(parameters[3],'float',-1,ehandler=self.ehandler)
        self.epsilon = dbPre(parameters[4],'float',ehandler=self.ehandler)
        self.quantumYield = dbPre(parameters[5],'float',ehandler=self.ehandler)
        self.source = parameters[6]
        
        if (self.source[0] =='X' and len(self.source)<4) or self.abswl == -1 or self.emwl == -1:
            self.ehandler.fatal( "Some mandatory parameters are not defined (X) in DB." )
            return False
        
        if len(parameters)==8:
            self.notes = parameters[7]
        else:
            self.notes = ""
        
        return True

##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit import PDBModel

class Test(BT.BiskitTest):
    """ Test cases for FRETEntity"""

    def prepare(self):
        EHandler.verbose = False
        EHandler.fails = False
        

    def cleanUp( self ):
        EHandler.verbose = True
        EHandler.fails = True


    def test_Creation(self):
        """Instantiation test """
        f = FRETEntity()
        self.assertEqual(f.name,"FRETEntity")

    def test_Loading(self):
        """Loading from database test """
        
        # Fully available data
        f = FRETEntity("LoadTest",T.testRoot() + "/fret/single.db",chromo_autodef = False)
        
        self.assertEqual (f.name, 'LoadTest' )
        self.assertEqual (f.lifetime, 1 )
        self.assertEqual (f.abswl, 2 )
        self.assertEqual (f.emwl, 3 )
        self.assertEqual (f.epsilon, 4 )
        self.assertEqual (f.quantumYield, 5 )
        self.assertEqual (f.source, 'FAKE' )
        
        # Not available data
        f = FRETEntity("foooo",T.testRoot() + "/fret/single.db",chromo_autodef = False)
        self.assertEqual ( "name not found" in f.ehandler.lastWarning,True) 
        
        # Not available db
        f = FRETEntity("LoadTest",T.testRoot() + "/fret/lol.db",chromo_autodef = False)
        self.assertEqual ( "empty or badly formatted database" in f.ehandler.lastWarning,True) 
        
        #Mandatory params
        
        f = FRETEntity("Mandatory1",T.testRoot() + "/fret/single.db",chromo_autodef = False)
        self.assertEqual ( "mandatory parameters are not defined" in f.ehandler.lastError,True) 
        
        f = FRETEntity("Mandatory2",T.testRoot() + "/fret/single.db",chromo_autodef = False)
        self.assertEqual ( "mandatory parameters are not defined" in f.ehandler.lastError,True) 

    def test_TransitionMoment(self):
        """ Calculating transition dipole moment test"""
        try:
            f = FRETEntity(name ='mCerulean',chromo_autodef = False, database = T.testRoot() + "/fret/single.db")
            f.chromophore = Chromophore("mCerulean","2Q57",T.testRoot() + "/fret/chromo.db")
        except:
            self.assertEqual (1,2)
        
        structure = PDBModel(T.testRoot() + "/fret/2Q57.pdb")
        
        f.structure = structure
        
        t = f.getChromoTransMoment()
        self.assertEqual (  t[0],2.0)
        self.assertEqual (  t[1],0.0)
        self.assertEqual (  t[2],7.0)

if __name__ == '__main__':

    BT.localTest()