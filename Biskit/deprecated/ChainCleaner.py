##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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
"""
Fill in missing atoms, and report low occupancies

Vintage code used by pdb2xplor script -- use L{Biskit.PDBCleaner} for
new implementations.
"""

from ChainSeparator import ChainSeparator
from Scientific.IO.PDB import *
import tools as T

class ChainCleaner:
    """
    Take with each call to next() one chain from ChainSeparator;
    Return this chain with completed residues.
    """

    def __init__(self, chainSeparator):
        """
        Clean up separate chains.
        
        @param chainSeparator: ChainSeparator
        @type  chainSeparator: module
        """
        self.reader = chainSeparator
        self.pdbname = self.reader.pdbname()   # take over pdb name
        self.log = self.reader.log             # take over log file

        self.aminoAcidDict = {'GLY':['N','CA','C','O'],
        'ALA':['N','CA','C','O','CB'],
        'VAL':['N','CA','C','O','CB','CG1','CG2'],
        'LEU':['N','CA','C','O','CB','CG','CD1','CD2'],
        'ILE':['N','CA','C','O','CB','CG1','CG2','CD1'],
        'MET':['N','CA','C','O','CB','CG','SD','CE'],
        'PRO':['N','CA','C','O','CB','CG','CD'],
        'PHE':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
        'TRP':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3',
               'CZ2','CZ3','CH2'],
        'SER':['N','CA','C','O','CB','OG'],
        'THR':['N','CA','C','O','CB','OG1','CG2'],
        'ASN':['N','CA','C','O','CB','CG','OD1','ND2'],
        'GLN':['N','CA','C','O','CB','CG','CD','OE1','NE2'],
        'TYR':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'],
        'CYS':['N','CA','C','O','CB','SG'],
        'LYS':['N','CA','C','O','CB','CG','CD','CE','NZ'],
        'ARG':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],
        'HIS':['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2'],
        'ASP':['N','CA','C','O','CB','CG','OD1','OD2'],
        'GLU':['N','CA','C','O','CB','CG','CD','OE1','OE2'],
        'ACE':['CA', 'C', 'O'],
        'NME':['N', 'CA'] }


    def _res2Terminal(self, aName_list):
        """
        Tweak list of allowed atom names to one expected for a
        C-terminal residue.

        @param aName_list: list of allowed atom names
                           e.g. ['N','CA','C','O','CB']
        @type  aName_list: list of strings
        
        @return: e.g. ['N','CA','C','OT1','CB','OT2']
        @rtype: list of str
        """
        result = []         # make local copy instead of taking reference
        result += aName_list
        try:
            result[result.index('O')] = 'OT1'
            result += ['OT2']
        except:
            pass  ## skip for CBX (Methyl amine)
        return result


    def _addMissing(self, residue, atom_name):
        """
        Add atom with given name to residue.
        
        @param residue: Scientific.IO.PDB.Residue object
        @type  residue: object
        @param atom_name: atom name
        @type  atom_name: str
        """
        if len(atom_name) < 2:
            # Atom.__init__ complaints about 1-char names
            atom_name = atom_name + ' '
        # new atom, element property needed for correct alignment in file
        residue.atom_list.append( Atom(atom_name, Vector(0,0,0),\
                                      element=atom_name[0] ) )


    def _completeResidues(self, chain):
        """
        Look for missing or unknown atom names, add missing atoms,
        report unknown atoms.
        
        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: object
        
        @return: Scientific.IO.PDB.PeptideChain object
        @rtype: chain object
        """
        chain.deleteHydrogens() ## delete all hydrogens
        i = 0
        self.log.add("Checking atoms of chain "+chain.segment_id)

        for res in chain:
            try:
                if i < len(chain)-1:            # normal residue
                    alowed = self.aminoAcidDict[res.name]
                else:                           # c-terminal residue
                    alowed = self._res2Terminal(self.aminoAcidDict[res.name])

                name_list = []

                for atom in res.atom_list:      # check for unknown atom names
                    # store for missing atom check
                    name_list = name_list + [atom.name]
                    if not (atom.name in alowed):
                        self.log.add('\tunknown atom: ' + atom.name + ' : '+\
                                     res.name+ str(res.number))

                for name in alowed:              # check for missing atoms
                    if not (name in name_list):
                        # add missing atom with 0 xyz
                        self._addMissing(res, name)  
                        self.log.add('\tadded missing atom -> '+ name+ ' : '+\
                                     res.name+ str(res.number))

            except:
               s = "\ncompleteResidues(): Error while checking atoms.\n"
               s = s + "residue " + str(i)+ " :"+ str(res) + "\n"
               s = s + T.lastError()
               T.errWriteln(
                   "Error while completing residues, check log for details.")
               self.log.add(s)

            i = i+1

        return chain


    def _checkOccupancy(self, chain):
        """
        Check and report atoms with ocupancies that is not 100% 
        Scientific.PDB.IO will only take one of the atoms even if there are
        alternate locations indicated in the PDB-file. The code below does only
        check for none 100% occupancies and report them to the log-file.

        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: chain object

        @return: Scientific.IO.PDB.PeptideChain object
        @rtype: chain object
        """
        self.log.add("Checking atom occupancies of chain "+chain.segment_id)
        for res in chain:
            for atom in res:
                if atom.properties.get('occupancy',1.0) != 1.0:
                    self.log.add('\tOccupancy: '+\
                                 str(atom.properties['occupancy']) \
                                 + ' : ' + atom.name + ' : '+ res.name+ ' ' +\
                                 str(res.number))
        return chain


    def _find_and_change(self, chain, oldAtomName, residueNum, newAtomName):
        """
        Change name of atoms in specified residue.

        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: chain object
        @param oldAtomName: atom name
        @type  oldAtomName: str
        @param residueNum: residue number
        @type  residueNum: int
        @param newAtomName: atom name
        @type  newAtomName: str
        
        @return: Scientific.IO.PDB.PeptideChain object
        @rtype: chain object
        """
        changeRes = chain.residues[residueNum]
        if changeRes.atoms.has_key(oldAtomName) == 1:
            changeRes.atoms[oldAtomName].name = newAtomName
        return chain


    def _correct_Cterm(self, chain):
        """
        Terminal amino acid can't have atom type O and OXT - have to be
        OT1 and OT2

        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: chain object

        @return: Scientific.IO.PDB.PeptideChain object
        @rtype: chain object
        """
        self._find_and_change(chain, 'O', -1, 'OT1')
        self._find_and_change(chain, 'OXT', -1, 'OT2')
        return chain


    def next(self):
        """
        Obtain next chain from ChainSeparator; Add missing atoms.
        
        @return: Scientific.IO.PDB.PeptideChain, completed chain OR
                 if no chain is left
        @rtype: chain object OR  None
        """
        chain = self.reader.next()
        if (chain is None):
            ## extract all waters into separate pdb
            self.reader.extractWaters( )

            return None
        self.log.add("\nCleaning up chain "+chain.segment_id+": ")

        # check for atoms that don't have 100% occupancy
        chain = self._checkOccupancy(chain) 

        # change OXT to OT2 and O to OT2 for C terminal
        chain = self._correct_Cterm(chain)  
        chain = self._completeResidues(chain)

        return chain


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """ Test ChainCleaner"""   

    def prepare(self):
        self.fname =   T.testRoot() + '/rec/1A2P_rec_original.pdb'
        self.outPath = T.tempDir()

    def cleanUp(self):
        T.tryRemove( self.cleaner.log.fname )
            

    def test_ChainCleaner( self):
        """ChainCleaner test"""
        self.cleaner = ChainCleaner( ChainSeparator(self.fname, self.outPath) )
        self.cleaner.next()
        
        if self.local:
            print 'Wrote log: %s'%(self.cleaner.log.fname)            

        ## return logfile contents
        f= open( self.cleaner.log.fname, 'r')
        self.result = f.readlines()
        f.close()

        self.assertEqual(self.result, self.EXPECTED)


    EXPECTED = ['\n', 'WARNING! The PDB-file contains coordinates for none water HETATMs.\n', 'If you want to keep the HETATM -  prepare the file for Xplor manualy \n', '\n', 'String found in line(s): \n', '\n', 'HETATM             \n', 'HETATM  881 ZN    ZN A 112      33.898  42.301  55.562  0.30 41.80          ZN  \n', 'HETATM 1891 ZN    ZN B 112       8.869  30.194  74.011  0.35 32.13          ZN  \n', 'HETATM 2908 ZN    ZN C 112      39.324  15.495  56.124  1.00 12.86          ZN  \n', '\n', '--------------------------------------------------------------------------------\n', '\n', 'Separate chains: \n', '------------------\n', '\n', "  Chain ID's of compared chains: ['A', 'B', 'C']\n", '  Cross-Identity between chains:\n', '[[ 1.  1.  1.]\n', ' [ 0.  1.  1.]\n', ' [ 0.  0.  1.]]\n', '  Identity threshold used: 0.9\n', 'NOTE: Identity matrix will be used for removing identical chains.\n', '2 chains have been removed.\n', '\n', '0 breaks found in chain (108 residues) A: []\n', 'changed segment ID of chain A to 1A2A\n', '\n', 'Cleaning up chain 1A2A: \n', 'Checking atom occupancies of chain 1A2A\n', '\tOccupancy: 0.3 : CD2 : TYR 17\n', '\tOccupancy: 0.8 : CE2 : TYR 17\n', '\tOccupancy: 0.5 : CZ : TYR 17\n', '\tOccupancy: 0.5 : OH : TYR 17\n', '\tOccupancy: 0.5 : OD1 : ASP 22\n', '\tOccupancy: 0.5 : OD2 : ASP 22\n', '\tOccupancy: 0.5 : CB : SER 28\n', '\tOccupancy: 0.5 : OG : SER 28\n', '\tOccupancy: 0.5 : CB : GLN 31\n', '\tOccupancy: 0.5 : CG : GLN 31\n', '\tOccupancy: 0.5 : CD : GLN 31\n', '\tOccupancy: 0.5 : OE1 : GLN 31\n', '\tOccupancy: 0.5 : NE2 : GLN 31\n', '\tOccupancy: 0.5 : CB : SER 38\n', '\tOccupancy: 0.5 : OG : SER 38\n', '\tOccupancy: 0.5 : N : ARG 59\n', '\tOccupancy: 0.5 : CA : ARG 59\n', '\tOccupancy: 0.5 : C : ARG 59\n', '\tOccupancy: 0.5 : O : ARG 59\n', '\tOccupancy: 0.0 : CB : ARG 59\n', '\tOccupancy: 0.0 : CG : ARG 59\n', '\tOccupancy: 0.0 : CD : ARG 59\n', '\tOccupancy: 0.0 : NE : ARG 59\n', '\tOccupancy: 0.0 : CZ : ARG 59\n', '\tOccupancy: 0.0 : NH1 : ARG 59\n', '\tOccupancy: 0.0 : NH2 : ARG 59\n', '\tOccupancy: 0.5 : OE1 : GLU 60\n', '\tOccupancy: 0.5 : OE2 : GLU 60\n', '\tOccupancy: 0.5 : CE : LYS 66\n', '\tOccupancy: 0.5 : NZ : LYS 66\n', '\tOccupancy: 0.5 : CB : SER 85\n', '\tOccupancy: 0.5 : OG : SER 85\n', '\tOccupancy: 0.5 : CB : ILE 96\n', '\tOccupancy: 0.5 : CG1 : ILE 96\n', '\tOccupancy: 0.5 : CG2 : ILE 96\n', '\tOccupancy: 0.5 : CD1 : ILE 96\n', 'Checking atoms of chain 1A2A\n']
    

if __name__ == '__main__':

    BT.localTest()

