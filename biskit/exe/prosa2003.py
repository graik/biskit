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
Analyze a structure using  Prosa2003.
"""
## allow relative imports when calling module as main script for testing https://www.python.org/dev/peps/pep-0366/
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import os.path
import time
import subprocess
import tempfile

import biskit.core.oldnumeric as N0
import biskit.tools as T
from biskit.errors import BiskitError

from .executor import Executor, TemplateError

class Prosa2003_Error( BiskitError ):
    pass


class Prosa2003( Executor ):
    """
    Analyze energies using Prosa2003 for a given PDBModels.
    =======================================================

    Reference
    ---------
       - U{http://lore.came.sbg.ac.at/}
       - U{http://www.proceryon.com/solutions/prosa.html}
       - Sippl,M.J. Recognition of Errors in Three-Dimensional
         Structures of Proteins.Proteins, 17: 355-362 (1993).

    Example usage
    -------------
        >>> x = Prosa2003( ref_model, [ m1, m2 ], verbose=1 )
        >>> result = x.run()

    Settings
    --------
       - I{lower_k = |integer| } and I{upper_k = |integer| }
          - Pair interaction energies are calculated for residue pairs
            whose distance k along the sequence is lower_k < k < upper_k.
            Default values are lower_k = 1, upper_k = 600. For example if
            you want to see only the short range energy contributions
            (e.g. sequence separation k < 9) set lower_k = 1, upper_k= 9.
            These parameters affect pair interactions only.

       - I{pot_lb = |float| } and I{pot_ub = |float| }
         - Pair energies are calculated in the distance range
           [pot_lb, pot_ub] A. Outside this range energies are zero.
           You can set the desired distance range of pair potential
           calculations. For example if you are not interested in
           energy contributions of close contacts, set pot lb = 4. The
           energy of pair interactions in the intervall [0; pot_lb] is
           then zero. Default: pot_lb = 0, pot_ub = 15. In addition the
           upper bound depends on the pair potentials used. The current
           maximum range is 15A.
           
    Structure of the input file 
    ---------------------------
    ::
       pair potential --- read in other than default potential
       surface potential
       pdb_dir --- directory with pdb file
       read pdb |filename| |object name|
       lower_k |int|
       upper_k |int|
       pot_lb |float|
       pot_ub |float|
       analyse energy |object name| 
       print energy |object name| |filename(.ana)|
       exit --- exit without launching graphical mode
    """ 


    inp = \
"""pair potential $PROSA_BASE/%(pairPot)s pcb
surface potential $PROSA_BASE/%(surfPot)s scb
pdb_dir = %(temp_dir)s
read pdb %(prosaPdbFile)s %(objectName)s 
lower_k = %(lower_k)i 
upper_k = %(upper_k)i 
pot_lb = %(pot_lb)3.1f 
pot_ub = %(pot_ub)3.1f 
analyse energy %(objectName)s 
print energy %(objectName)s %(prosaOutput)s
exit\n
"""

    def __init__(self, models, **kw ):
        """
        @param models: if more than one model is given they are
                       concatenated and the energy is calculated
                       for the two together.
        @type  models: PDBModels
                 

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

        self.models = models

        ## Potentials to use, pII3.0 is the default setting
        self.pairPot = 'prosa2003.pair-cb' # default: pII3.0.pair-cb 
        self.surfPot = 'prosa2003.surf-cb' # default: pII3.0.surf-cb

        ## temp files for prosa pdb file and prosa output file
        self.prosaPdbFile = tempfile.mktemp('_prosa2003.pdb')
        self.prosaOutput = tempfile.mktemp('_prosa2003.out')
        self.temp_dir = T.tempDir()

        prosaInput = tempfile.mktemp('_prosa2003.inp')

        ## set default values
        self.objectName = 'obj1'
        self.lower_k = 1
        self.upper_k = 600
        self.pot_lb = 0.
        self.pot_ub = 15.

        Executor.__init__( self, 'prosa2003', template=self.inp,
                           f_in=prosaInput, **kw )

        ## check the path to the potential files
        self.checkPotentials()


    def execute( self, inp=None ):
        """
        Run Prosa2003.
        
        @note: Was forced to overwrite execute function of Executor to get
        prosa2003 running (the code below equals an os.system call).

        @return: duration of calculation in seconds
        @rtype: float

        @raise Prosa2003_Error: if execution failed
        """
        start_time = time.time()

        cmd = self.command()

        shellexe = None
        if self.exe.shell and self.exe.shellexe:
            shellexe = self.exe.shellexe

        stdin = stdout = stderr = None

        if self.exe.pipes:
            stdin = subprocess.PIPE
            stdout= subprocess.PIPE
            stderr= subprocess.PIPE
        else:
            inp = None
            if self.f_in:
                stdin = open( self.f_in )
            if self.f_out:
                stdout= open( self.f_out, 'w' )
            if self.f_err and self.catch_err:    
                stderr= open( self.f_err, 'w' )                

        if self.verbose:
            self.log.add('executing: %s' % cmd)
            self.log.add('in folder: %s' % self.cwd ) 
            self.log.add('input:  %r' % stdin )
            self.log.add('output: %r' % stdout )
            self.log.add('errors: %r' % stderr )
            self.log.add('wrapped: %r'% self.exe.shell )
            self.log.add('shell: %r'  % shellexe )
            self.log.add('environment: %r' % self.environment() )
            if self.exe.pipes and inp:
                self.log.add('%i byte of input pipe' % len(str(inp)))

        try:
            retcode = subprocess.call( '%s %s > %s'\
                                       %(self.exe.bin, self.f_in, self.f_out),
                                       shell=True)
            if retcode < 0:
                raise Prosa2003_Error("Child was terminated by signal" \
                      + str(retcode))

        except Prosa2003_Error as why:
            raise Prosa2003_Error("Execution of Prosa2003 failed: " + str(why))

        return time.time() - start_time


    def checkPotentials( self ):
        """
        Check that the environment variable $PROSA_BASE is correct
        and that the potential files are there.
        """
        sysPotDir = os.environ['PROSA_BASE']
        altPotDir = '/Bis/shared/centos-3/prosa2003/ProSaData/'

        if not os.path.exists( sysPotDir + self.pairPot + '.frq' ):
            if os.path.exists( altPotDir + self.pairPot + '.frq' ):
                os.environ['PROSA_BASE'] = altPotDir

        if not os.path.exists( sysPotDir + self.surfPot + '.sff' ):
            if os.path.exists( altPotDir + self.surfPot + '.sff' ):
                os.environ['PROSA_BASE'] = altPotDir


    def prepare( self ):
        """
        Make and write a PROSA compatible pdb file.
          - consecutive residue numbering
          - same chainId troughout the model

        @note: Overrides Executor method.
        """
        model=self.models[0]

        if len(self.models)> 1:
            for i in range( 1, len( self.models )):
                model.concat(self.models[i])

        model.renumberResidues()
        #for a in model.getAtoms():
            #a['chain_id'] = 'P'
        model['chain_id'] = ['P'] * len( model )
        model.writePdb( self.prosaPdbFile, ter=0 )


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.prosaPdbFile )
            T.tryRemove( self.f_in)
            T.tryRemove( self.prosaOutput + '.ana' )


    def parse_result( self ):
        """
        Parse the Prosa2003 output file.
        
        @return: dictionary with the calculated potential profiles
                 and the parameters used
        @rtype: dict
        """
        prosa_pair = [] 
        prosa_surf = []
        prosa_tot = []  

        prosaout = self.prosaOutput + '.ana'

        lines = []
        try:
            lines = open( prosaout ).readlines()
            if not lines:
                raise IOError('File %s is empty'%( prosaout ))
        except IOError as why:
            raise IOError("Couldn't read Prosa result: " + str( why ) \
                  + '\n Check the Prosa license!')

        ## comment lines starts with '#'
        for i in range( len(lines) ):
            if lines[i][0] != '#':
                nr, pair, surf, tot = str.split( lines[i])
                prosa_pair += [ float( pair ) ]
                prosa_surf += [ float( surf ) ]
                prosa_tot += [ float( tot ) ]

        ## create dictionary with residue profiles and calc. info
        result = {'prosa_pair':N0.array(prosa_pair),
                  'prosa_surf':N0.array(prosa_surf),
                  'prosa_tot':N0.array(prosa_tot),
                  'ProsaInfo':{ 'lower_k':self.lower_k,
                                'upper_k':self.upper_k,
                                'pot_lb':self.pot_lb,
                                'pot_ub':self.pot_ub } }

        return result


    def prosaEnergy( self ):
        """
        Sum of energy profiles.
        
        @return: sum of the three energy profiles
                 [ E_prosa_pair, E_prosa_surf, E_prosa_tot]
        @rtype: [float]
        """
        ## calc. energies
        r = []
        for k in ['prosa_pair', 'prosa_surf', 'prosa_tot']:
            r += [ self.result[k] ]

        return N0.sum( r, 1 )


    def isFailed( self ):
        """
        @note: Overrides Executor method
        """
        return not self.error is None 


    def finish( self ):
        """
        @note: Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_result( )


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.EXE, BT.OLD ]

    def test_Prosa2003(self):
        """Prosa2003 test"""
        from biskit import PDBModel
        
        ## Loading PDB...
        self.ml = PDBModel( T.testRoot()+'/lig/1A19.pdb' )
        self.ml = self.ml.compress( self.ml.maskProtein() )

        self.mr = PDBModel( T.testRoot()+'/rec/1A2P.pdb' )
        self.mr = self.mr.compress( self.mr.maskProtein() )

        ## Starting Prosa2003
        self.prosa = Prosa2003( [self.ml, self.mr], debug=0, verbose=0 )

        ## Running
        self.ene = self.prosa.run()

        self.result = self.prosa.prosaEnergy()

        if self.local:
            print("Result: ", self.result)

        self.assertTrue( N0.sum(self.result - [ -94.568,  -64.903, -159.463 ] ) \
                      < 0.0000001 )

        
if __name__ == '__main__':

    BT.localTest()
