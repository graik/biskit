## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

## class ChainSeperator:
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## $Revision$
## last $Author$
## last $Date$

"""
Seperate PDB into continuous peptide chains for XPlor. Remove duplicate
peptide chains. Required by pdb2xplor.py

This is vintage code. See L{Biskit.PDBCleaner} for a more recent
version (though yet lacking some functions).

@todo: Create an override for the chain comparison if one wants
          to keep identical chains (i.e homodimers)
"""

## from Blast2Seq import *   # compare 2 sequences
from molUtils import singleAA
import Biskit.tools as T
from LogFile import LogFile

from Scientific.IO.PDB import *
import numpy as N
import string
from difflib import SequenceMatcher
import re

class ChainSeparator:
    """
    Open PDB file; give back one chain whenever next() is
    called. This class is used by the pdb2xplor script.

    This class constitutes vintage code. See
    L{Biskit.PDBCleaner} and L{Biskit.Mod.TemplateCleaner} for a more
    recent implementation of PDB cleaning.

    @todo: The removal of duplicate chains should be transferred to
    the PDBCleaner so that this class can be retired
    """

    def __init__(self, fname, outPath='', chainIdOffset=0,
                 capBreaks=0, chainMask=0, log=None ):
        """
        @param fname: pdb filename
        @type  fname: str
        @param outPath: path for log file
        @type  outPath: str
        @param chainIdOffset: start chain numbering at this offset
        @type  chainIdOffset: int
        @param capBreaks: add ACE and NME to N- and C-term. of chain breaks [0]
        @type  capBreaks: 0|1
        @param chainMask: chain mask for overriding the default sequence identity [None]
        @type  chainMask: [1|0]
        @param log: LogFile object
        @type  log: object
        """
        self.pdb = Structure(fname);
        self.fname = fname
        self.outPath = T.absfile( outPath )
        self.chainIdOffset = chainIdOffset
        self.capBreaks = capBreaks
        self.log = LogFile( T.absfile(outPath)+'/' + self.pdbname()+'.log')
        if log:
            self.log = log

        self.chains = self.pdb.peptide_chains
        self.counter = -1
        self.threshold = 0.9 # sequence identity between multiple copies in PDB
        self._expressionCheck(
            "[^\n].*[Hh][Oo][Mm][Oo].?[Dd][Ii][Mm][eE][Rr].*\n", 'HOMODIMER')
        self._expressionCheck("[^\n].*[Tt][Rr][Ii][Mm][Ee][Rr].*\n", 'TRIMER')
        self._hetatomCheck()

        self.log.add("Separate chains: \n------------------")
        self._removeDuplicateChains(chainMask)   # keep only one copy of molecule
        self._separateChainBreaks()
        self._assign_seg_ids()          # new segment id for each chain


    def pdbname(self):
        """
        Extract pdb code from file name.
        
        @return: (assumed) pdb code
        @rtype: str
        """
        return T.stripFilename(self.pdb.filename)


    def _expressionCheck(self, findExpression, findClean):
        """
        Check and report if the regular expression 'findExpression'
        exists in the PDB-file. Use this to locate data in the REMARK
        section of a pdb file. Prints a warning to stdOut if the
        regular expression is found.

        @param findExpression: regular expression
        @type  findExpression: str
        @param findClean: clean name of regular expression
        @type  findClean: str
        """
        pdb = open(self.fname,'r')
        pdbFile = pdb.read()
        searchResult = re.findall(findExpression,pdbFile)

        warningMessage = """
WARNINGR! The text string'%s' was found in the PDB-file.
If this PDB-file contains a homodimer one of the chains will be
deleted by this script. To avoid this prepare the file for Xplor manualy \n""" %\
        ( findClean )
        warningMessage2 = """--------------------------------------------\n"""

        if len(searchResult) != 0:
            self.log.add(warningMessage)
            self.log.add("String found in line(s): \n")
            for i in range(0,len(searchResult)):
                self.log.add(searchResult[i])
            self.log.add(warningMessage2)
        pdb.close() 


    def _hetatomCheck(self):
        """
        Check and report if there are any none-water HETATMs in the PDB-file
        """
        pdb = open(self.fname,'r')
        pdbFile = pdb.read()
        findExpression = "HETATM.*\n"
        searchResult = re.findall(findExpression,pdbFile)
        i=0
        j = len(searchResult)
        while i<j:
            if searchResult[i][17:20] == "HOH" or \
               searchResult[i][0:6] != "HETATM" :
                del searchResult[i]
                i=i-1
                j=j-1
            i=i+1

        warningMessage = """
WARNING! The PDB-file contains coordinates for none water HETATMs.
If you want to keep the HETATM -  prepare the file for Xplor manualy \n"""
        warningMessage2 = "\n"+ 80*"-" + "\n"
        if len(searchResult) != 0:
            self.log.add(warningMessage)
            self.log.add("String found in line(s): \n")
            for i in range(0,len(searchResult)):
               self.log.add(searchResult[i][0:-1])
            self.log.add(warningMessage2)
        pdb.close()


    def _compareSequences( self, seq1, seq2 ):
        """
        @param seq1: sequence 1 to compare
        @type  seq1: str
        @param seq2: sequence 1 to compare
        @type  seq2: str
        @return: identity (0.0 - 1.0) between the two sequences
        @rtype : float
        """
        # compare the 2 sequences
##        blast = Blast2Seq( seq1, seq2 )
##        id = blast.run()
        matcher = SequenceMatcher( None, ''.join(seq1) , ''.join(seq2) )
        return matcher.ratio()


    def _removeDuplicateChains(self, chainMask=None):
        """
        Get rid of identical chains by comparing all chains with Blast2seq.

        @param chainMask: chain mask for overriding the
                          chain identity checking (default: None)
        @type  chainMask: [int]
        
        @return: number of chains removed
        @rtype: int
        """
        chainCount = len(self.chains)
        matrix = 1.0 * N.zeros((chainCount,chainCount))
        chain_ids = []

        ## create identity matrix for all chains against all chains
        for i in range(0, chainCount):
            chain_ids = chain_ids + [self.chains[i].chain_id] # collect for log file
            for j in range(i, len(self.chains)):

                # convert 3-letter-code res list into 1-letter-code String
                seq1 = singleAA( self.chains[i].sequence() )
                seq2 = singleAA( self.chains[j].sequence() )

##                 if len(seq1) > len(seq2):           # take shorter sequence
##                 # aln len at least half the len of the shortest sequence
##                     alnCutoff = len(seq2) * 0.5     
##                 else:
##                     alnCutoff = len(seq1) * 0.5
##                 if id['aln_len'] > alnCutoff:
##                     matrix[i,j] = id['aln_id']
##                 else:                           # aln length too short, ignore
##                     matrix[i,j] = 0

                matrix[i,j] = self._compareSequences( seq1, seq2 )

        ## report activity
        self.log.add("\n  Chain ID's of compared chains: "+str(chain_ids))
        self.log.add("  Cross-Identity between chains:\n"+str(matrix))
        self.log.add("  Identity threshold used: "+str(self.threshold))
        
        ## override the automatic chain deletion by supplying a
        ## chain mask to this function
        if chainMask:
            if len(chainMask) == chainCount:
                self.chains = N.compress(chainMask, self.chains)
                self.log.add("NOTE: chain mask %s used for removing chains.\n"%chainMask)
           
            else:
                self.log.add("########## ERROR ###############")
                self.log.add("# Chain mask is only %i chains long"%len(chainMask))
                self.log.add("# when a mask of length %i is needed"%chainCount)
                self.log.add("# No cleaning will be performed.\n")

        if not chainMask:
            ## look at diagonals in "identity matrix"
            ## (each chain against each)
            duplicate = len(self.chains)
            for offset in range(1,chainCount):
                diag = N.diagonal(matrix, offset ,0,1)
                # diagonal of 1's mark begin of duplicate
                avg = 1.0 * N.sum(diag)/len(diag)
                if (avg >= self.threshold):
                    duplicate = offset
                    break
            self.chains = self.chains[:duplicate]
            self.log.add("NOTE: Identity matrix will be used for removing identical chains.")

        ## report activit
        self.log.add(str(chainCount - len(self.chains))+\
                     " chains have been removed.\n")
        
        # how many chains have been removed?
        return (chainCount - len(self.chains))


    def _assign_seg_ids(self):
        """
        Assign new segment id to each chain.
        """
        counter = self.chainIdOffset
        for chain in self.chains:

            ## Assemble segid from pdb code + one letter out of A to Z
            chain.segment_id = self.pdbname()[:3] + string.uppercase[counter]
            counter = counter + 1
            try:                        # report changed segement ids
                chain_id = chain.chain_id
                self.log.add("changed segment ID of chain "+chain_id+\
                             " to "+chain.segment_id)
            except:
                T.errWriteln("_assign_seg_ids(): logerror")


    def _sequentialDist(self, chain, cutoff, atom):
        """
        Calculate sequential atom-atom distance, report residues with
        longer distance than cutoff (chain break positions).
        
        @param chain: Scientific.IO.PDB.PeptideChain object
        @type  chain: object
        @param cutoff: threshold for reporting gap (chain break)
        @type  cutoff: float
        @param atom: type of atoms to check (i.e. 'CA')
        @type  atom: str

        @return: list of chain break positions (residue index for each
                 first residue of two that are too distant)
        @rtype: list of int           
        """        
        distanceList = []
        v0 = Vector( 0,0,0 )
        jump = 1

        for res in range(0,len(chain)-2):

            try:
                v1 = Vector(chain[res][atom].position.array)

                ## ignore CA with 0,0,0 coordinate
                if v1 != v0:

                    jump = 1
                    v2 = Vector(chain[ res+jump ][atom].position.array)

                    ## look for next CA with non-zero coordinate
                    while v2 == v0 and jump + res < len( chain ):
                        jump += 1
                        v2 = Vector(chain[ res+jump ][atom].position.array)

                    if (v1 - v2).length() > cutoff * jump:
                        distanceList = distanceList + [res + jump - 1]

            except:
                self.log.add(
                    "_sequentialDist():\nError while checking CA-CA distance"+\
                    " between residues "+str(chain[res].name)+\
                    str(chain[res].number)+" and "+\
                    str(chain[res+jump].name)+\
                    str(chain[res+jump].number)+ " in chain "+chain.chain_id)
                self.log.add("Error: " + T.lastError() )

        return distanceList

##     def _sequentialDist(self, chain, cutoff, atom):
##         """
##         Calculate sequential atom-atom distance, report residues with
##         longer distance than cutoff (chain break positions).
##         chain  - PDB.PeptideChain
##         cutoff - float, threshold for reporting gap (chain break)
##         atom   - str, type of atoms to check (i.e. 'CA')

##         -> [int, int, ...], list of chain break positions (residue index
##                             for each first residue of two that are too distant)
##         """        
##         distanceList = []        
##         for residue in range(0,len(chain)-1):
##             # iterate through residue 1 to ter-1
##             try:
##                 vectorAtom1 = Vector(chain[residue][atom].position.array)
##                 vectorAtom2 = Vector(chain[residue+1][atom].position.array)

##                 if (vectorAtom1 - vectorAtom2).length() > cutoff:
##                     distanceList = distanceList + [residue]
##             except:
##                 self.log.add(
##                     "_sequentialDist():\nError while checking CA-CA distance"+ \
##                     " between residues "+str(chain[residue].name)+\
##                     str(chain[residue].number)+" and "+str(chain[residue+1].name)+\
##                     str(chain[residue+1].number)+ " in chain "+chain.chain_id)
##                 self.log.add("Error: " + T.lastError() )

##         return distanceList


    def _separateChainBreaks(self):
        """
        Separate chains with breaks into 2 chains.
        The new chain(s) is/are added to the internal PDB instance
        (self.chains).
        """
        fragments = []

        for chain in self.chains:
            # res number of residues before a break
            breaks = self._sequentialDist(chain, 4.5, 'CA')
            self.log.add(str(len(breaks)) + " breaks found in chain " +\
                         "(" + str(len(chain)) \
                         + " residues) " + chain.chain_id + ": "+str(breaks))

            previous = 0
            ncap_next = 0
            for breakRes in breaks:

                residues = chain.residues[previous:breakRes+1]
                previous = breakRes + 1

                chainNew = PeptideChain(residues, chain.chain_id,
                                        chain.segment_id)
                if ncap_next:
                    self.__nCap( chainNew )
                    ncap_next = 0

                if self.capBreaks:
                    ## add N-Methyl to c terminal
                    self.__cCap( chainNew )

                    ncap_next = 1

                fragments = fragments + [chainNew]

            chainNew = PeptideChain(chain.residues[previous:], chain.chain_id,
                                    chain.segment_id)
            if ncap_next:
                self.__nCap( chainNew )

            fragments = fragments + [chainNew]

        self.chains = fragments


    def __nCap( self, pep_chain ):
        """
        Add acetyl capping to N-terminal of peptide chain
        """
        n = (pep_chain[0].number or 1) - 1

        r = AminoAcidResidue('ACE', number=n, atoms=[Atom('CA', Vector(0,0,0),
                                                          element='C')])
        pep_chain.residues = [r] + pep_chain.residues

        self.log.add('Capping chain break with ACE %i' % n) 


    def __cCap( self, pep_chain ):
        """
        Add methyle amine capping to C-terminal of peptide chain
        """
        n = (pep_chain[-1].number or len(pep_chain)) + 1

        r = AminoAcidResidue('NME', number=n, atoms=[Atom('CA', Vector(0,0,0),
                                                          element='C')])
        pep_chain.residues = pep_chain.residues + [r]

        self.log.add('Capping chain break at with NME %i' % n) 


    def extractWaters(self):
        """
        Write waters into separate pdb file, called |pdbCode|_waters.pdb.
        """
        try:
            fTarget = self.outPath + '/' +\
                      self.pdbname()[:4] + '_waters.pdb'
            pdb = PDBFile( fTarget, mode='w' )

            waters = []
            for key in ['HOH', 'DOD']:

                if self.pdb.molecules.has_key( key ):

                    waters += self.pdb.molecules[ key ]

            pdb.nextChain(chain_id='', segment_id='1XWW')
            for w in waters:
                pdb.nextResidue('TIP3')
                ## XPLOR wants "ATOM" not "HETATM":
                pdb.het_flag = 0
                pdb.writeAtom('OH2', w.atoms['O'].position)

            ## keep TIP3 waters as well
            if len(waters) == 0:
                try:
                    TIP3_waters = self.pdb.molecules[ 'TIP3' ]
                except:
                    TIP3_waters = []
                    
                for w in TIP3_waters:
                    pdb.nextResidue('TIP3')
                    ## XPLOR wants "ATOM" not "HETATM":
                    pdb.het_flag = 0
                    pdb.writeAtom('OH2', w.atoms['OH2'].position)
                    pdb.writeAtom('H1', w.atoms['H1'].position)
                    pdb.writeAtom('H2', w.atoms['H2'].position)
            pdb.close()
            
        except:
            T.errWriteln("Error writing waters to %s: " % fTarget )
            T.errWriteln( T.lastError() )



    def next(self):
        """
        Return next 'clean', non-redundant, non-broken chain from PDB

        @return: Scientific.IO.PDB.PeptideChain, completed chain OR
                 if no chain is left
        @rtype: chain object OR None        
        """
        self.counter = self.counter + 1
        if (len(self.chains) > self.counter):
            return self.chains[self.counter]
        else:
            return None


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test ChainSeparator """

    def prepare(self):
        self.fname =   T.testRoot() + '/com/1BGS_original.pdb'
        self.outPath = T.tempDir()

    def cleanUp(self):
        T.tryRemove( self.sep.log.fname ) 

    def test_ChainSeparator( self ):
        """ChainSeparator test"""
        self.sep = ChainSeparator( self.fname, self.outPath, 1)  

        self.chain = self.sep.next()

        i=1
        all_chains = []
        while self.chain <> None:
            if self.local:
                print 'Chain %i:'%i, ''.join(singleAA(self.chain.sequence() ) )
            all_chains += self.chain.sequence()
            self.chain = self.sep.next()
            i += 1

        if self.local:
            print 'ChainSeparator log file written to: %s'%self.sep.log.fname

        r = ''.join( singleAA( all_chains ) )
        self.assertEqual(r, self.EXPECTED)


    EXPECTED='AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIRAQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIRAQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIRKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILSKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILSKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS'

        

if __name__ == '__main__':

    BT.localTest()





        
