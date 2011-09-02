##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg
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

## last $Author: graik $
## last $Date: 2009-05-09 14:17:28 +0200 (Sat, 09 May 2009) $
## $Revision: $
"""
Wrapper for structure protonation program reduce

Application: 
"""

import tempfile, re
import numpy as N

from Biskit import Executor, TemplateError, PDBModel
## import Biskit.settings as S
import Biskit.tools as T

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
    
    Usage
    =====

    >>> red = Reduce( model )
    >>> result = red.run()

    @note: Command configuration: biskit/Biskit/data/defaults/exe_reduce.dat
    """

    def __init__( self, model, **kw ):
        """
        @param model: structure to be aligned to reference
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
        self.f_pdbin = tempfile.mktemp( '_reduce_in.pdb' )
        f_out= tempfile.mktemp( '_reduce_out.pdb')
        self.f_db = T.dataRoot() + '/reduce/reduce_wwPDB_het_dict.txt' 

        kw['catch_err'] = True
        
        Executor.__init__( self, 'reduce', 
                           args= '-BUILD -DB %s %s'%(self.f_db,self.f_pdbin),
                           f_out=f_out,
                           **kw )

        self.model = model

    def version(self):
        return 'Reduce $Revision: $'

    def prepare( self ):
        """
        Overrides Executor method.
        """
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

        raise ReduceError, s

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
            + self.version()


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    def test_reduce( self ):
        """TMAlign test"""
        if self.local: print 'Loading PDB...'

        self.m1 = PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )


        if self.local: print 'Starting Reduce'
        self.x = Reduce( self.m1, debug=self.DEBUG,
                         verbose=self.local )

        if self.local:
            print 'Running'

        self.r = self.x.run()

        if self.local:
            print "Result: "
            print self.r.report()



if __name__ == '__main__':

    BT.localTest(debug=False)