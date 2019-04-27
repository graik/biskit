## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
##

"""
Clean PDB-files so that they can be used for MD. This module is a 
(still partial) re-implementation of the vintage pdb2xplor script.
"""

import biskit.molUtils as MU
import biskit.mathUtils as M
import biskit.tools as t
from biskit import PDBModel, StdLog

import biskit.core.oldnumeric as N0

import copy

class CleanerError( Exception ):
    pass

class CappingError( CleanerError ):
    pass

class PDBCleaner:
    """
    PDBCleaner performs the following tasks:
    
      * remove HETAtoms from PDB
      * replace non-standard AA by its closest standard AA
      * remove non-standard atoms from standard AA residues
      * delete atoms that follow missing atoms (in a chain)
      * remove multiple occupancy atoms (except the one with highest occupancy)
      * add ACE and NME capping residues to C- and N-terminals or chain breaks
        (see capTerminals(), this is NOT done automatically in process())

    Usage:
    =======

      >>> c = PDBCleaner( model )
      >>> c.process()
      >>> c.capTerminals( auto=True )

    This will modify the model in-place and report changes to STDOUT.
    Alternatively, you can specify a log file instance for the output.
    PDBCleaner.process accepts several options to modify the processing.
    
    Capping
    =======
    
    Capping will add N-methyl groups to free C-terminal carboxy ends
    or Acetyl groups to free N-terminal Amines and will thus 'simulate' the
    continuation of the protein chain -- a common practice in order to 
    prevent fake terminal charges. The automatic discovery of missing residues
    is guess work at best. The more conservative approach is to use,
    for example:
    
      >>> c.capTerminals( breaks=1, capC=[0], capN=[2] )
      
    In this case, only the chain break detection is used for automatic capping
    -- the last residue before a chain break is capped with NME and the first
    residue after the chain break is capped with ACE. Chain break detection
    relies on PDBModel.chainBreaks() (via PDBModel.chainIndex( breaks=1 )).
    The normal terminals to be capped are now specified explicitely. The first
    chain (not counting chain breaks) will receive a NME C-terminal cap and the
    third chain of the PDB will receive a N-terminal ACE cap. 
    
    Note: Dictionaries with standard residues and atom content are defined
          in Biskit.molUtils. This is a duplicate effort with the new strategy
          to parse Amber prep files for very similar information
          (AmberResidueType, AmberResidueLibrary) and should change once we 
          implement a real framework for better residue handling. 
    """
    
    #: these atoms always occur at the tip of of a chain or within a ring
    #: and, if missing, will not trigger the removal of other atoms
    TOLERATE_MISSING = ['O', 'CG2', 'CD1', 'CD2', 'OG1', 'OE1', 'NH1',
                        'OD1', 'OE1',
                        'H5T',"O5'", ]

    ## PDB with ACE capping residue
    F_ace_cap = t.dataRoot() + '/amber/leap/ace_cap.pdb'
    ## PDB with NME capping residue
    F_nme_cap = t.dataRoot() + '/amber/leap/nme_cap.pdb'

    def __init__( self, fpdb, log=None, verbose=True ):
        """
        :param fpdb: pdb file OR PDBModel instance
        :type  fpdb: str OR Biskit.PDBModel
        :param log: biskit.LogFile object (default: STDOUT)
        :type  log: biskit.LogFile
        :param verbose: log warnings and infos (default: True)
        :type  verbose: bool
        """
        self.model = PDBModel( fpdb )
        self.log = log or StdLog()
        self.verbose = verbose


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print(msg)

    def remove_multi_occupancies( self ):
        """
        Keep only atoms with alternate A field (well, or no alternate).
        """
        if self.verbose:
            self.logWrite( self.model.pdbCode +
                           ': Removing multiple occupancies of atoms ...')

        i = 0
        to_be_removed = []

        for a in self.model:

            if a['alternate']:
                try:
                    str_id = "%i %s %s %i" % (a['serial_number'], a['name'],
                                              a['residue_name'],
                                              a['residue_number'])

                    if a['alternate'].upper() in ['A', '1']:
                        a['alternate'] = ''

                    else:
                        if float( a['occupancy'] ) < 1.0:
                            to_be_removed += [ i ]
                            if self.verbose:
                                self.logWrite(
                                    'removing %s (%s %s)' %
                                    (str_id,a['alternate'], a['occupancy']))
                        else:
                            if self.verbose:
                                self.logWrite(
                                 ('keeping non-A duplicate %s because of 1.0 '+
                                  'occupancy') % str_id )

                except:
                    self.logWrite("Error removing duplicate: "+t.lastError() )
            i+=1

        try:
            self.model.remove( to_be_removed )
            if self.verbose:
                self.logWrite('Removed %i atoms' % len( to_be_removed ) )

        except:
            if self.verbose:
                self.logWrite('No atoms with multiple occupancies to remove' )


    def replace_non_standard_AA( self, amber=0, keep=[] ):
        """
        Replace amino acids with none standard names with standard
        amino acids according to :class:`MU.nonStandardAA`
        
        :param amber: don't rename HID, HIE, HIP, CYX, NME, ACE [0]
        :type  amber: 1||0
        :param keep: names of additional residues to keep
        :type keep:  [ str ]
        """
        standard = list(MU.atomDic.keys()) + keep

        if amber:
            standard.extend( ['HID', 'HIE', 'HIP', 'CYX', 'NME', 'ACE'] )

        replaced = 0

        if self.verbose:
            self.logWrite(self.model.pdbCode +
                          ': Looking for non-standard residue names...')

        resnames = self.model['residue_name']
        for i in self.model.atomRange():

            resname = resnames[i].upper()

            if resname not in standard:
                if resname in MU.nonStandardAA:
                    resnames[i] = MU.nonStandardAA[ resname ]

                    if self.verbose:
                        self.logWrite('renamed %s %i to %s' % \
                                     (resname, i, MU.nonStandardAA[ resname ]))
                else:
                    resnames[i] = 'ALA'

                    self.logWrite('Warning: unknown residue name %s %i: ' \
                                  % (resname, i ) )
                    if self.verbose:
                        self.logWrite('\t->renamed to ALA.')

                replaced += 1

        if self.verbose:
            self.logWrite('Found %i atoms with non-standard residue names.'% \
                          replaced )


    def __standard_res( self, resname, amber=0 ):
        """
        Check if resname is a standard residue (according to :class:`MU.atomDic`)
        if not return the closest standard residue (according to
        :class:`MU.nonStandardAA`).
        
        :param resname: 3-letter residue name
        :type  resname: str
        
        :return: name of closest standard residue or resname itself
        :rtype: str
        """
        if resname in MU.atomDic:
            return resname

        if resname in MU.nonStandardAA:
            return MU.nonStandardAA[ resname ]

        return resname


    def remove_non_standard_atoms( self ):
        """
        First missing standard atom triggers removal of standard atoms that
        follow in the standard order. All non-standard atoms are removed too.
        Data about standard atoms are taken from :class:`MU.atomDic` and symomym
        atom name is defined in :class:`MU.atomSynonyms`.
        
        :return: number of atoms removed
        :rtype: int
        """
        mask = []
        
        if self.verbose:
            self.logWrite("Checking content of standard amino-acids...")

        for res in self.model.resList():

            resname  = self.__standard_res( res[0]['residue_name'] ).upper()
            if resname == 'DC5':
                pass
            
            ## bugfix: ignore non-standard residues that have no matching 
            ## standard residue
            if resname in MU.atomDic:
                
                standard = copy.copy( MU.atomDic[ resname ] )
    
                ## replace known synonyms by standard atom name
                for a in res:
                    n = a['name']
                    if not n in standard and MU.atomSynonyms.get(n,0) in standard:
                        a['name'] = MU.atomSynonyms[n]
                        if self.verbose:
                            self.logWrite('%s: renaming %s to %s in %s %i' %\
                                          ( self.model.pdbCode, n, a['name'],
                                           a['residue_name'], a['residue_number']))
    
                anames   = [ a['name'] for a in res ]
                keep = 1
    
                ## kick out all standard atoms that follow a missing one
                rm = []
                for n in standard:
                    if (not n in anames) and not (n in self.TOLERATE_MISSING):
                        keep = 0
    
                    if not keep:
                        rm += [ n ]
    
                for n in rm:
                    standard.remove( n )
    
                ## keep only atoms that are standard (and not kicked out above)
                for a in res:
    
                    if a['name'] not in standard:
                        mask += [1]
                        if self.verbose:
                            self.logWrite('%s: removing atom %s in %s %i '%\
                                          ( self.model.pdbCode, a['name'],
                                           a['residue_name'], a['residue_number']))
                    else:
                        mask += [0]

        self.model.remove( mask )
        
        if self.verbose:
            self.logWrite('Removed ' + str(N0.sum(mask)) +
                          ' atoms because they were non-standard' +
                          ' or followed a missing atom.' )

        return N0.sum( mask )

    def capACE( self, model, chain, breaks=True, checkgap=True ):
        """
        Cap N-terminal of given chain.

        Note: In order to allow the capping of chain breaks,
        the chain index is, by default, based on model.chainIndex(breaks=True), 
        that means with chain break detection activated! This is not the 
        default behaviour of PDBModel.chainIndex or takeChains or chainLength. 
        Please use the wrapping method capTerminals() for more convenient 
        handling of the index.

        :param model: model
        :type  model: PDBMode
        :param chain: index of chain to be capped
        :type  chain: int
        :param breaks: consider chain breaks when identifying chain boundaries
        :type  breaks: bool
        
        :return: model with added NME capping
        :rtype : PDBModel
        """
        if self.verbose:
            self.logWrite('Capping N-terminal of chain %i with ACE' % chain )

        c_start = model.chainIndex( breaks=breaks )
        c_end = model.chainEndIndex( breaks=breaks)
        Nterm_is_break = False
        Cterm_is_break = False
        
        if breaks:
            Nterm_is_break = c_start[chain] not in model.chainIndex()
            Cterm_is_break = c_end[chain] not in model.chainEndIndex()
            
        m_ace = PDBModel( self.F_ace_cap )

        chains_before = model.takeChains( range(chain), breaks=breaks )
        m_chain       = model.takeChains( [chain], breaks=breaks )
        chains_after  = model.takeChains( range(chain+1, len(c_start)),
                                          breaks=breaks )

        m_term  = m_chain.resModels()[0]

        ## we need 3 atoms for superposition, CB might mess things up but
        ## could help if there is no HN
        ##        if 'HN' in m_term.atomNames():
        m_ace.remove( ['CB'] )  ## use backbone 'C' rather than CB for fitting 

        ## rename overhanging residue in cap PDB
        for a in m_ace:
            if a['residue_name'] != 'ACE':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]-1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## fit cap onto first residue of chain
        m_ace = m_ace.magicFit( m_term )

        cap = m_ace.resModels()[0]
        serial = m_term['serial_number'][0] - len(cap)
        cap['serial_number'] = list(range( serial, serial + len(cap)))

        ## concat cap on chain
        m_chain = cap.concat( m_chain, newChain=False )

        ## re-assemble whole model
        r = chains_before.concat( m_chain, newChain=not Nterm_is_break)
            
        r = r.concat( chains_after, newChain=not Cterm_is_break)
        
        if checkgap and len(c_start) != r.lenChains( breaks=breaks ):
            raise CappingError('Capping ACE would mask a chain break. '+\
                  'This typically indicates a tight gap with high risk of '+\
                  'clashes and other issues.')

        return r


    def capNME( self, model, chain, breaks=True, checkgap=True ):
        """
        Cap C-terminal of given chain. 

        Note: In order to allow the capping of chain breaks,
        the chain index is, by default, based on model.chainIndex(breaks=True), 
        that means with chain break detection activated! This is not the 
        default behaviour of PDBModel.chainIndex or takeChains or chainLength.
        Please use the wrapping method capTerminals() for more convenient 
        handling of the index.

        :param model: model
        :type  model: PDBMode
        :param chain: index of chain to be capped
        :type  chain: int
        :param breaks: consider chain breaks when identifying chain boundaries
        :type  breaks: bool
        
        :return: model with added NME capping residue
        :rtype : PDBModel
        """
        if self.verbose:
            self.logWrite('Capping C-terminal of chain %i with NME.' % chain )
        m_nme   = PDBModel( self.F_nme_cap )

        c_start = model.chainIndex( breaks=breaks )
        c_end = model.chainEndIndex( breaks=breaks)

        Nterm_is_break = False
        Cterm_is_break = False
        if breaks:
            Nterm_is_break = c_start[chain] not in model.chainIndex()
            Cterm_is_break = c_end[chain] not in model.chainEndIndex()
         
        chains_before = model.takeChains( range(chain), breaks=breaks )
        m_chain       = model.takeChains( [chain], breaks=breaks )
        chains_after  = model.takeChains( range(chain+1, len(c_start)),
                                          breaks=breaks )

        m_term  = m_chain.resModels()[-1]

        ## rename overhanging residue in cap PDB, renumber cap residue
        for a in m_nme:
            if a['residue_name'] != 'NME':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]+1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## chain should not have any terminal O after capping
        m_chain.remove( ['OXT'] )            

        ## fit cap onto last residue of chain
        m_nme = m_nme.magicFit( m_term )
        
        cap = m_nme.resModels()[-1]
        serial = m_term['serial_number'][-1]+1
        cap['serial_number'] = list(range( serial, serial + len(cap)))

        ## concat cap on chain
        m_chain = m_chain.concat( cap, newChain=False )

        ## should be obsolete now
        if getattr( m_chain, '_PDBModel__terAtoms', []) != []:
            m_chain._PDBModel__terAtoms = [ len( m_chain ) - 1 ]
        assert m_chain.lenChains() == 1

        ## re-assemble whole model
        r = chains_before.concat( m_chain, newChain=not Nterm_is_break)
        r = r.concat( chains_after, newChain=not Cterm_is_break)

        if checkgap and len(c_start) != r.lenChains( breaks=breaks ):
            raise CappingError('Capping NME would mask a chain break. '+\
                  'This typically indicates a tight gap with high risk of '+\
                  'clashes and other issues.')
        
        return r


    def convertChainIdsNter( self, model, chains ):
        """
        Convert normal chain ids to chain ids considering chain breaks.
        """
        if len(chains) == 0: 
            return chains
        i = N0.take( model.chainIndex(), chains ) 
        ## convert back to chain indices but this time including chain breaks
        return model.atom2chainIndices( i, breaks=1 )
        
    def convertChainIdsCter( self, model, chains ):
        """
        Convert normal chain ids to chain ids considering chain breaks.
        """
        if len(chains) == 0: 
            return chains
        ## fetch last atom of given chains
        index = N0.concatenate( (model.chainIndex(), [len(model)]) )
        i = N0.take( index, N0.array( chains ) + 1 ) - 1
        ## convert back to chain indices but this time including chain breaks
        return model.atom2chainIndices( i, breaks=1 )
    

    def unresolvedTerminals( self, model ):
        """
        Autodetect (aka "guess") which N- and C-terminals are most likely not
        the real end of each chain. This guess work is based on residue 
        numbering:
        
        * unresolved N-terminal: a protein residue with a residue number > 1

        * unresolved C-terminal: a protein residue that does not contain either
                               OXT or OT or OT1 or OT2 atoms
                               
        :param model: PDBModel
        
        :return: chains with unresolved N-term, with unresolved C-term
        :rtype : ([int], [int])
        """
        c_first = model.chainIndex()
        c_last  = model.chainEndIndex()
        
        capN = [ i for (i,pos) in enumerate(c_first)\
                 if model['residue_number'][pos] > 1 ]
        
        capN = [i for i in capN if model['residue_name'][c_first[i]] != 'ACE']
        
        capN = self.filterProteinChains( model, capN, c_first )
        
        capC = []
        for (i,pos) in enumerate(c_last):
            atoms = model.takeResidues(model.atom2resIndices([pos])).atomNames()
            
            if not( 'OXT' in atoms or 'OT' in atoms or 'OT1' in atoms or \
                    'OT2' in atoms ):
                capC += [ i ]

        capC = self.filterProteinChains( model, capC, c_last )
                  
        return capN, capC
    
    #@todo filter for protein positions in breaks=1

    def filterProteinChains( self, model, chains, chainindex ):
        maskProtein = model.maskProtein()
        chains = [ i for i in chains if maskProtein[ chainindex[i] ] ]
        return chains

    def capTerminals( self, auto=False, breaks=False, capN=[], capC=[],
                      checkgap=True):
        """
        Add NME and ACE capping residues to chain breaks or normal N- and 
        C-terminals. Note: these capping residues contain hydrogen atoms.
        
        Chain indices for capN and capC arguments can be interpreted either
        with or without chain break detection enabled. For example, let's
        assume we have a two-chain protein with some missing residues (chain
        break) in the first chain:
        
        A:   MGSKVSK---FLNAGSK
        B:   FGHLAKSDAK

        Then:
          capTerminals( breaks=False, capN=[1], capC=[1]) will add N-and 
          C-terminal caps to chain B.
        However:
          capTerminals( breaks=True, capN=[1], capC=[1]) will add N- and 
          C-terminal caps to the second fragment of chain A.
          
        
        Note: this operation *replaces* the internal model.
        
        :param auto: put ACE and NME capping residue on chain breaks
                     and on suspected false N- and C-termini (default: False)
        :type  auto: bool
        :param breaks: switch on chain break detection before interpreting
                       capN and capC
        :type  breaks: False
        :param capN: indices of chains that should get ACE cap (default: [])
        :type  capN: [int]
        :param capC: indices of chains that should get NME cap (default: [])
        :type  capC: [int]
        """
        m = self.model
        c_len = m.lenChains()
        i_breaks = m.chainBreaks()
            
        if auto:
            if not breaks:
                capN = self.convertChainIdsNter( m, capN )
                capC = self.convertChainIdsCter( m, capC )

            breaks=True
            capN, capC = self.unresolvedTerminals( m )
        
            end_broken = m.atom2chainIndices( m.chainBreaks(), breaks=1 )
            
            capC = M.union( capC, end_broken )
            capN = M.union( capN, N0.array( end_broken ) + 1 )
            
        capN = self.filterProteinChains(m, capN, m.chainIndex(breaks=breaks))
        capC = self.filterProteinChains(m, capC, m.chainEndIndex(breaks=breaks))

        for i in capN:
            m = self.capACE( m, i, breaks=breaks, checkgap=checkgap )
            assert m.lenChains() == c_len, '%i != %i' % \
                   (m.lenChains(), c_len)
            assert len(m.chainBreaks(force=True)) == len(i_breaks)
            assert m['serial_number'].dtype == N0.Int32, 'serial_number not int'

        for i in capC:
            m = self.capNME( m, i, breaks=breaks, checkgap=checkgap )
            assert m.lenChains() == c_len
            assert len(m.chainBreaks(force=True)) == len(i_breaks)
        
        self.model = m
        return self.model

    
    
    def process( self, keep_hetatoms=0, amber=0, keep_xaa=[] ):
        """
        Remove Hetatoms, waters. Replace non-standard names.
        Remove non-standard atoms.
        
        :param keep_hetatoms: option
        :type  keep_hetatoms: 0||1
        :param amber: don't rename amber residue names (HIE, HID, CYX,..)
        :type  amber: 0||1
        :param keep_xaa: names of non-standard residues to be kept
        :type  keep_xaa: [ str ]
        
        :return: PDBModel (reference to internal)
        :rtype: PDBModel
        
        :raise CleanerError: if something doesn't go as expected ...
        """
        try:
            if not keep_hetatoms:
                self.model.remove( self.model.maskHetatm() )

            self.model.remove( self.model.maskH2O() )

            self.model.remove( self.model.maskH() )

            self.remove_multi_occupancies()

            self.replace_non_standard_AA( amber=amber, keep=keep_xaa )

            self.remove_non_standard_atoms()


        except KeyboardInterrupt as why:
            raise KeyboardInterrupt( why )
        except Exception as why:
            self.logWrite('Error: '+t.lastErrorTrace())
            raise CleanerError( 'Error cleaning model: %r' % why )

        return self.model



#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    def prepare(self):
        from biskit import LogFile
        import tempfile


    def test_PDBCleaner( self ):
        """PDBCleaner general test"""
        
        ## Loading PDB...
        self.c = PDBCleaner( t.testRoot() + '/rec/1A2P_rec_original.pdb',
                             log=self.log,
                             verbose=self.local)
        
        self.m = self.c.process()

        self.assertAlmostEqual( self.m.mass(), 34029.0115499993, 7 )
        
    def test_DNACleaning( self ):
        """PDBCleaner DNA test"""
        ## Loading PDB...
        self.c = PDBCleaner( t.testRoot() + 'amber/entropy/0_com.pdb',
                             log=self.log, verbose=self.local )
        
        self.dna = self.c.process(amber=True)

        self.assertAlmostEqual( self.dna.mass(), 26953.26, 1 )
        
        
    def test_Capping( self ):
        """PDBCleaner.capTerminals test"""
        ## Loading PDB...
        self.model = PDBModel(t.testRoot() + '/rec/1A2P_rec_original.pdb')

        self.c = PDBCleaner( self.model, log=self.log, verbose=self.local )       
        self.m2 = self.c.capTerminals( breaks=True )
        self.assertTrue( self.m2.atomNames() == self.model.atomNames() )
        
        self.m3 = self.model.clone()
        self.m3.removeRes( [10,11,12,13,14,15] )
        self.m4 = self.m3.clone()
        
        self.c = PDBCleaner( self.m3, log=self.log, verbose=self.local )
        self.m3 = self.c.capTerminals( breaks=True, capC=[0], capN=[0,1])
        self.assertEqual( self.m3.takeChains([0]).sequence()[:18], 
                          'XVINTFDGVADXXKLPDN' )
        
        if self.local:
            self.log.add( '\nTesting automatic chain capping...\n' )
        
        self.c = PDBCleaner( self.m4, log=self.log, verbose=self.local )
        self.m4 = self.c.capTerminals( auto=True )
        self.assertEqual( self.m4.takeChains([0]).sequence()[:18], 
                          'XVINTFDGVADXXKLPDN' )
        
        
    def test_capping_internal(self):
        self.m3 = PDBModel(t.testRoot('pdbclean/foldx_citche.pdb'))
        self.m3 = self.m3.takeChains([1])  # pick second chain; has chain break
        self.c3 = PDBCleaner(self.m3, verbose=self.local, log=self.log) 
        
        m = self.c3.capACE(self.m3, 1, checkgap=False)
        
        self.assertEqual(m.lenChains(breaks=True) -1, 
                         self.m3.lenChains(breaks=True))
        
class FailingTest(BT.BiskitTest):
    """Test class """
    
    TAGS = [BT.NORMAL, BT.FAILS]
        
    def test_capping_extra( self ):
        """PDBCleaner.capTerminals extra challenge"""
        self.m2 = PDBModel( t.testRoot() + '/pdbclean/foldx_citche.pdb' )
        self.c = PDBCleaner( self.m2, verbose=self.local, log=self.log)
        self.assertRaises(CappingError, self.c.capTerminals, auto=True)
        if self.local:
            self.log.add('OK: CappingError has been raised indicating clash.' )
        
        self.assertEqual( len(self.m2.takeChains([1]).chainBreaks()), 1 )

if __name__ == '__main__':

    BT.localTest()
