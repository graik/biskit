## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

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
## last $Author$
## last $Date$
## $Revision$

"""
Calculates the secondary structure using DSSP.
"""

import tempfile
import numpy as N
from Biskit import Executor, TemplateError
import Biskit.tools as T
import Biskit.molUtils as MU
from Errors import BiskitError

class Dssp_Error( BiskitError ):
    pass


class Dssp( Executor ):
    """
    Run Dssp
    ========
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
#        self.model = model.clone( deepcopy=1 )

        ## temporary pdb-file
        self.f_pdb = tempfile.mktemp( '_dssp.pdb')
        self.f_out = tempfile.mktemp( '_dssp.out')

        Executor.__init__( self, 'dsspcmbi',
                           args='-na %s'%self.f_pdb,
                           catch_err=1, **kw )


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.model = self.model.compress( self.model.maskHeavy() )
        if self.model.lenAtoms() == N.sum(self.model.maskCA):
            raise Dssp_Error, 'The structure you want to calculate the secondary structure for seems to be a carbon alpha trace. Terminating'
        self.model.writePdb( self.f_pdb )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )


    def parse_result( self ):
        """
        Parse the secondary structure from tha DSSP output file.

        @return: a list with the standard DSSP secondary structure
                 description (one letter code) with the exception that
                 a blank " " has been replaced by a dot "."
        @rtype: [str]
        """
        ## check that the outfut file is there and seems valid
        try:
            out_file = open( self.f_out )
            lines = out_file.readlines()
            out_file.close()
        except:
            raise Dssp_Error,\
                  'Dssp result file %s does not exist.'%self.f_out
        if len(lines) == 0:
            raise Dssp_Error,\
                  'Dssp result file %s empty'%self.f_out
        if len(lines) < 9:
            raise Dssp_Error,\
                  'Dssp result file %s contains no secondary structure data'%self.f_out

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

        ## don't parse the header
        for i in range( len(lines) ):
            if lines[i][:12]=='  #  RESIDUE':
                start = i+1

        ## collect DSSP data
        ss, seq, term = [], [], []
        for i in range( start, len(lines) ):
            ss   += [ lines[i][16:17] ]
            term += [ lines[i][14:15] ]
            seq  += [ lines[i][13:14] ]

        def __completeBB( res ):
            """
            Check that residue have all backbone atoms
            CA, N, C and O or OXT
            """
            atoms = [ a['name'] for a in res ]
            count = atoms.count('CA') + atoms.count('N') + \
                  atoms.count('C') + atoms.count('O')
            if count == 4:
                return 1

        secStruc = []

        resDic = self.model.resList()
        i = 0
        j = 0
        while i<len(ss) or j<len(resDic):

            complete = __completeBB( resDic[j] )
##            res_name = MU.singleAA( [resDic[j][0]['residue_name']] )[0]

            ## assign irregular if not complete residue, DSSP
            ## skipps these residues
            if not complete:
                secStruc += ['.']
                j += 1

            ## termini, only in DSSP output
            elif ( seq[i] == '!' ) and ( term[i] == '*' ):
                i += 1

            ## chain break, only in DSSP output
            elif seq[i] == '!':
                i += 1

            ## normal data
            elif seq[i] != '!':
                ## replace ' ' with '.'
                if ss[i] == ' ':
                    secStruc += ['.']
                else:
                    secStruc += [ss[i]]
                i += 1
                j += 1

        ## check that the entire sequence has a secondary structure assigned
        assert len(secStruc) == self.model.lenResidues()

        return ''.join(secStruc)        


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_result( )



#############
##  TESTING        
#############

import Biskit.test as BT
class Test(BT.BiskitTest):
    """DSSP test"""

    TAGS = [BT.EXE]

    def prepare(self):
        self.f = T.testRoot()+"/com/1BGS.pdb"


    def test_DSSP( self ):
        """DSSP test"""

        from Biskit import PDBModel

        if self.local: print 'Loading PDB...'
        self.m = PDBModel(self.f)
        self.m = self.m.compress( self.m.maskProtein() )

        if self.local:  print 'Starting DSSP'
        self.dssp = Dssp( self.m )

        if self.local: print 'Running DSSP'

        self.result = self.dssp.run()

        if self.local:
            print "Sequence :", self.m.sequence()
            print "Secondary:", self.result            

        self.assertEquals( self.result, self.EXPECTED)


    EXPECTED =  '.....SHHHHHHHHHHHSS..TTEE.HHHHHHHT..GGGT.HHHHSTT.EEEEEEE..TT..S...TT..EEEEE.S..SSS..S.EEEEETT..EEEESSSSSS.EE...EEEEETTT..SHHHHHHHHHHHHT..TT..SSHHHHHHHHHHT..SSEEEEEE.HHHHHHHTTTTHHHHHHHHHHHHHHT..EEEEE.'

if __name__ == '__main__':

    BT.localTest()
