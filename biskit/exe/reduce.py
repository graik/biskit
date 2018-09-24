##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg
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

"""
Wrapper for structure protonation program reduce
"""
## allow relative imports when calling module as main script for testing https://www.python.org/dev/peps/pep-0366/
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import tempfile
import numpy as N

import biskit
from biskit import PDBModel, PDBCleaner
import biskit.tools as T

from .executor import Executor

class ReduceError( Exception ):
    pass

class Reduce( Executor ):
    """
    A reduce wrapper. Creates a new PDBModel with hydrogens added and some
    polar groups flipped so as to optimize the H-bond network. The result 
    structure (in Reduce.result) is a new PDBModel with, potentially,
    some more atom changes than only added H. All existing H atoms are stripped
    before reduce is called.
    
    Profiles associated to residues are carried over into the result model.
    The method Reduce.merge() can be used to rescue atom profiles  
    from the input model, but use with care because the atom content has
    changed.
    
    The informative stderr output of the program is captured in Reduce.error.
    The run is reported as failed (and an exception is raised) if the keyword
    "ERROR" appears in this program output.
    
    Note: reduce appears to loosely follow the Amber convention of naming
    H atoms (e.g. HD12). You can switch to the xplor convention (2HD1) with
    PDBModel.amber2xplor().
    
    Capping of terminals
    =====================
    
    See L{Biskit.PDBCleaner} for a discussion of protein chain capping.
    
    Usage
    =====

    >>> red = Reduce( model )
    >>> result = red.run()

    @note: Command configuration: biskit/Biskit/data/defaults/exe_reduce.dat
    """

    def __init__( self, model, tempdir=None, args='', 
                  autocap=False, capN=[], capC=[],
                  **kw ):
        """
        @param model: input structure
        @type  model: PDBModel
        @param tempdir: create dedicated temporary folder (default: None)
                        see Executor
        @param tempdir: str | 0|1
        @param args: additional command line arguments for reduce (default:'')
                     example: '-OLDpdb'
        @type  args: str
        
        @param autocap: add capping NME and ACE residues to any (auto-detected)
                        false N- or C-terminal (default: False)
        @type  autocap: bool
        
        @param capN: cap N-terminal of these chains (indices) with ACE ([])
        @type  capN: [ int ]
        
        @param capC: cap C-terminal of these chains (indices) with NME ([])
        @type  capN: [ int ]

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
        tempdir = self.newtempfolder( tempdir )
        
        self.f_pdbin = tempfile.mktemp( '_in.pdb', 'reduce_', dir=tempdir )
        f_out= tempfile.mktemp( '_out.pdb', 'reduce_', dir=tempdir)
        self.f_db = T.dataRoot() + '/reduce/reduce_wwPDB_het_dict.txt' 
        
        self.autocap = autocap
        self.capN = capN
        self.capC = capC
        
        lenchains = model.lenChains()

        Executor.__init__( self, 'reduce', 
                           args= '%s -BUILD -Nterm%i -DB %s %s' %\
                           (args, lenchains, self.f_db, self.f_pdbin),
                           f_out=f_out, catch_err=True,
                           tempdir=tempdir,
                           **kw )

        self.model = model

    def capTerminals( self ):
        c = PDBCleaner( self.model, verbose=self.verbose )
        self.model = c.capTerminals( auto=self.autocap, 
                                     capN=self.capN, capC=self.capC )
    
    def prepare( self ):
        """
        Overrides Executor method.
        """
        if self.autocap or len(self.capN) or len(self.capC):
            self.capTerminals()
            
        self.model = self.model.compress( self.model.maskHeavy() )
        self.model.writePdb( self.f_pdbin )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        Does nothing. No temporary files are created.
        """        
        Executor.cleanup( self )
        if not self.debug:
            T.tryRemove( self.f_pdbin )
            T.tryRemove( self.f_out)


    def isFailed( self ):
        """
        Overrides Executor method
        """
        for l in self.error:
            if 'ERROR' in l:
                return True
        return False

    def fail( self ):
        """
        Overrides Executor method. Called when execution fails.
        """
        s = 'reduce failed. Please check the program output in the '+\
          'field `output` of this Reduce instance (e.g. `print x.output`)!'
        self.log.add( s )

        raise ReduceError(s)

    def merge( self, inmodel=None ):
        """
        Rescue profiles and info records from input model into result model.
        """
        r = self.result.clone()
        m = inmodel or self.model
        i1, i2 = r.compareAtoms( m )

        mask1 = N.zeros( len( r ), N.int )
        N.put( mask1, i1, 1 )

        r.atoms.update( m.atoms, mask=mask1 )
        r.fileName = m.fileName
        r.source = m.source
        
        return r


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = PDBModel( self.f_out )

        ## renumber atoms 
        self.result['serial_number'] = N.arange( len( self.result ) )
        ## rescue non-atom informations
        self.result.pdbCode = self.model.pdbCode
        self.result.info.update( self.model.info )
        self.result.residues.update( self.model.residues )
        self.result.info['reduce'] = 'hydrogens added/replaced by '\
            + 'Reduce v' + biskit.__version__


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    def test_reduce( self ):
        """Reduce test"""
        if self.local: self.log.add('Loading PDB...')

        self.m1 = PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        
        rec = T.load( T.testRoot('com/rec.model'))
        lig = T.load( T.testRoot('com/lig.model'))
        
        self.m2 = rec.concat(lig)

        if self.local: self.log.add('Starting Reduce')
        self.x = Reduce( self.m1, debug=self.DEBUG,
                         verbose=self.local )
        if self.local:
            self.log.add( 'Running')
        self.r = self.x.run()

        if self.local:
            self.log.add("Result: ")
            self.log.add(self.r.report(prnt=False))

        if self.local:
            self.log.add('Reduce protein complex')
        self.x = Reduce( self.m2, debug=self.DEBUG, verbose=self.local,
                         autocap=True )
        self.r2 = self.x.run()
        

if __name__ == '__main__':

    BT.localTest(debug=False)
