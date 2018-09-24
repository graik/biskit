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

"""
Calculates the secondary structure using DSSP.
"""

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import tempfile
import biskit.core.oldnumeric as N0
import numpy as N

import biskit.tools as T
import biskit.molUtils as MU
from biskit.errors import BiskitError
from biskit import EHandler

from .executor import Executor, TemplateError

class ParseDSSPTable(object):
    """Parse DSSP result table assuming fixed-width columns"""
    
    ## (<name>, <start>, <stop>, <conversion_function>)
    fields = [('i',0,5,int), 
              ('residue_number',5,10,int),('insertion_code', 10,11,None),
              ('chain',11,12,None), 
              ('residue',12,14,None),
              ('ss',16,17,None), ('acc',35,39,int),
              ('phi',103,109,float), ('psi',109,115,float),
              ('x',115,122,float), ('y',122,129,float), ('z',129,136,float)]


    def parse_line(self, l):
        d = {}
        try:
            for f in self.fields:
                name, start, end, conversion = f
                if conversion is None:
                    d[name] = l[start:end].strip()
                else:
                    d[name] = conversion( l[start:end].strip() )
                    
        except:
            return {}
            
        return d

class Dssp_Error( BiskitError ):
    pass

class Dssp( Executor ):
    """
    Run Dssp
    ========
    
    Tested with mkdssp version 2.2.1
    
    The DSSP program will define the secondary structure of a given
    structure. The secondary structure elements defined are::

      H = alpha helix
      B = residue in isolated beta-bridge
      E = extended strand, participates in beta ladder
      G = 3-helix (3/10 helix)
      I = 5 helix (pi helix)
      T = hydrogen bonded turn
      S = bend
      . = loop or irregular

    Example usage
    -------------
        >>> d = Dssp( model )
        >>> result = d.run()
        >>> result['dssp']
           ['.','.','E','E','E','.','H','H',...]
        
    Returns a copy of the input model with four new residue profiles:
        - dssp -- the "summary" secondary structure code for each residue
        - dssp_acc -- accessible surface area calculated by dssp
        - dssp_phi -- PHI angle (as defined by IUPAC), if undefined: 360.0
        - dssp_psi -- PSI angle (as defined by IUPAC), if undefined: 360.0
    
    Residues without DSSP data are marked as '.' in the 'dssp' profile.
    
    Note: additional residue profiles are "residue", "chain", "insertion_code"
    which are used by DSSP to match result residue entries to the source model.

    References
    ----------
       - U{http://swift.cmbi.ru.nl/gv/dssp/}
       - Kabsch W, Sander C (1983) Dictionary of protein secondary
         structure: pattern recognition of hydrogen-bonded and geometrical
         features. Biopolymers Dec;22(12):2577-637. 
    """
    
    def __init__( self, model, **kw ):
        """
        @param model: model analyze
        @type  model: PDBModel

        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        self.model = model

        ## temporary pdb-file
        self.f_pdb = tempfile.mktemp( '_dssp.pdb')
#        self.f_out = tempfile.mktemp( '_dssp.out')
        
        ## will hold content of result file
        self.raw_result = []

        Executor.__init__( self, 'dsspcmbi',
                           args='-i %s'%self.f_pdb,
                           catch_err=1, **kw )


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.model = self.model.compress( self.model.maskHeavy() )
        if self.model.lenAtoms() == N0.sum(self.model.maskCA):
            raise Dssp_Error('The structure you want to calculate the secondary structure for seems to be a carbon alpha trace. Terminating')
        self.model.writePdb( self.f_pdb )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )


    def parse_file( self ):
        """
        Parse the secondary structure from tha DSSP output file.

        @return: a list with the standard DSSP secondary structure
                 description (one letter code) with the exception that
                 a blank " " has been replaced by a dot "."
        @rtype: [str]
        """
        ## check that the outfut is there and seems valid
        try:
            lines = self.output.split('\n')
        except:
            raise Dssp_Error('Dssp result missing.')
        if len(lines) == 0:
            raise Dssp_Error('Dssp result empty')
        if len(lines) < 9:
            raise Dssp_Error('Dssp result contains no secondary structure data')
        
        return lines

    def parse_lines(self, lines):

        ## Collect secondary structure data. Note that:
        ##
        ## 1. If Dssp detects a chain break or a residue with an
        ##    incomplete backbone it inserts a line with an
        ##    exclamation mark "!" in the residue column.
        ##
        ##    Example: 496  524 C V    <     
        ##             497  525 C P            
        ##             498        !             
        ##             499  539 C F             
        ##             500  540 C N     >  -    
        ##        
        ## 2. A termini (a switch of chain identifier) is marked
        ##    by a asterisk "*" in the column next to the "!"
        ##
        ##    Example: 721  789 C Q  T <       
        ##             722  790 C A    <       
        ##             723        !*          
        ##             724   33 Z E     >     
        ##             725   34 Z K  T  4  +    
        ##             726   35 Z L  T >> S+   
        ##
        ## 3. If DSSP detects an incomplete residue in the terminus
        ##    we only get a single line with a "!" and "*" not a
        ##    for the incomplete residue and one for the terminus
        ##    as woulld be expected!
        ##
        ##    Example: 84  101 A A           
        ##             85  102 A E         
        ##             86        !*     <-- residue 103 skipped!     
        ##             87   16 B S           
        ##             88   17 B G  E

        ## skip header (may be interesting to get total area and bonds though)
        while lines and not lines[0][:12]=='  #  RESIDUE':
            lines.pop(0)
        lines.pop(0)

        table = ParseDSSPTable()
        r = [ table.parse_line(l) for l in lines ]
        r = [ res for res in r if res ]  ## filter out TER and invalid records
        
        m = self.model
        len_res = m.lenResidues()
        
        m.residues.set('dssp', ['.']*len_res ,default='.')
        m.residues.set('dssp_acc', N.zeros(len_res), default=0.0)
        m.residues.set('dssp_phi', N.ones(len_res)*360.0, default=360.0 )
        m.residues.set('dssp_psi', N.ones(len_res)*360.0, default=360.0 )
        m.residues['residue'] = MU.singleAA( m.atom2resProfile('residue_name'),
                                             unknown='X')
        m.residues['insertion_code'] = m.atom2resProfile('insertion_code')
        m.residues['chain'] = m.atom2resProfile('chain_id')

        i = 0
        for res in r:
            while i < (len_res-1) and \
                  not ( m.residues['residue'][i] == res['residue'] and\
                        m.residues['chain'][i] == res['chain'] and\
                        m.residues['insertion_code'][i] == res['insertion_code'] ):
                i += 1
            m.residues['dssp'][i] = res['ss'] if res['ss'] else '.'
            m.residues['dssp_acc'][i] = res['acc']
            m.residues['dssp_phi'][i] = res['phi']
            m.residues['dssp_psi'][i] = res['psi']
            
            i += 1


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.raw_result = self.parse_file()
        self.parse_lines( self.raw_result )
        self.result = self.model


#############
##  TESTING        
#############

import biskit.test as BT
class Test(BT.BiskitTest):
    """DSSP test"""

    TAGS = [BT.EXE]

    def prepare(self):
        self.f = T.testRoot()+"/com/1BGS.pdb"
        
    def generic_test(self, pdb, expected=None, proteinonly=False):
        """generic DSSP test"""

        from biskit import PDBModel

        if self.local: print('Loading PDB...')
        self.m = PDBModel(pdb)
        
        if proteinonly:
            self.m = self.m.compress( self.m.maskProtein() )

        if self.local:  print('Starting DSSP')
        self.dssp = Dssp( self.m, verbose=self.local, debug=self.DEBUG )

        if self.local: print('Running DSSP')

        self.result = self.dssp.run()  ## returns modified PDBModel
        self.result = self.result.compress( self.result.maskProtein() )
        self.result = ''.join(self.result['dssp'])

        if self.local:
            print("Sequence :", self.m.sequence())
            print("Secondary:", self.result)            
        
        if expected:
            self.assertEqual( self.result, expected)        


    EXPECTED =  '.....SHHHHHHHHHHHSS..TTEE.HHHHHHHT..GGGT.HHHHSTT.EEEEEEE..TT..S...TT..EEEEE.S..SSS..S.EEEEETT..EEEESSSSSS.EE...EEEEETTT..SHHHHHHHHHHHHT..TT..SSHHHHHHHHHHT..SSEEEEEE.HHHHHHHTTTTHHHHHHHHHHHHHHT..EEEEE.'
    def test_DSSP( self ):
        """DSSP test"""
        self.generic_test(self.f, expected=self.EXPECTED, proteinonly=False)

    def test_DSSP_2W3A( self ):
        """DSSP test 2W3A"""
        self.generic_test('2W3A', proteinonly=False)

    def test_DSSP_3TGI( self ):
        """DSSP test 3TGI"""
        ## handling of residues with insertion code (184A, 186A, 221A)
        self.generic_test('3TGI', proteinonly=False)

    def test_DSSP_1R4Q( self ):
        """DSSP test"""
        self.generic_test('1R4Q', proteinonly=False)    

    def test_DSSP_1AM7(self):
        self.generic_test('1AM7')

    def test_DSSP_5FB8(self):
        self.generic_test('5fb8')
        
    def test_DSSP_3j4q(self):
        self.generic_test(T.testRoot('dssp') + '/3j4q.pdb')

class TestLong(BT.BiskitTest):
    """DSSP testing full list of trouble-makers"""

    TAGS = [BT.EXE, BT.LONG, BT.EXTRA]

    def test_DSSP_alltroublemakers(self):
        from biskit import PDBModel
        lst = '''1c48.pdb  3rlg.pdb  4gv5.pdb  4tsp.pdb
1p9g.pdb  3v03.pdb  4i0n.pdb  4tsq.pdb
1rd9.pdb  3v09.pdb  4idj.pdb  4wdc.pdb
1rdp.pdb  3vwi.pdb  4iw1.pdb  4whn.pdb
2anv.pdb  3w9p.pdb  4iya.pdb  4wvp.pdb
2anx.pdb  3weu.pdb  4j0l.pdb  4wx3.pdb
2f9r.pdb  3wev.pdb  4j2v.pdb  4wx5.pdb
2j8f.pdb  3wmu.pdb  4jk4.pdb  4zbq.pdb
2wc8.pdb  3wmv.pdb  4nzl.pdb  4zbr.pdb
2wcb.pdb  4bct.pdb  4ot2.pdb  5abw.pdb
2wce.pdb  4e99.pdb  4pho.pdb  5bpg.pdb
2xrs.pdb  4emx.pdb  4phq.pdb  5cr6.pdb
2xsc.pdb  4f5s.pdb  4po0.pdb  5dby.pdb
3a57.pdb  4f5t.pdb  4q7g.pdb  5elc.pdb
3hdf.pdb  4f5u.pdb  4rw3.pdb  5ele.pdb
3kve.pdb  4f5v.pdb  4rw5.pdb  5elf.pdb
3lu6.pdb  4fba.pdb  4tsl.pdb  5id7.pdb
3lu8.pdb  4fzm.pdb  4tsn.pdb  5id9.pdb
3mxg.pdb  4fzn.pdb  4tso.pdb'''.split()

        for s in lst:
            dssp = Dssp(PDBModel(s[:4]))
            r = dssp.run()

            if self.local:
                print(s, ':')
                print(r.sequence())
                print(''.join( r.compress(r.maskProtein())['dssp'] ))
                print()

if __name__ == '__main__':

    BT.localTest()
