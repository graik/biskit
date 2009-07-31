##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2009-07-02 17:36:36 +0200 (Thu, 02 Jul 2009) $
## $Revision: 810 $
"""Wrap TMAlign structure alignment program"""

import tempfile
from Biskit import Executor, TemplateError, BiskitError, PDBModel
import Biskit.tools as T
import Biskit.molUtils as MU
import re
import numpy as N

## Example output of TMAlign (for parsing)
"""

 **************************************************************************
 *                               TM-align                                 *
 * A protein structural alignment algorithm based on TM-score             *
 * Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9  *
 * Comments on the program, please email to: yzhang@ku.edu                *
 **************************************************************************

Chain 1:1MQ1_A.pdb  Size=  91
Chain 2:1MWN_A.pdb  Size=  91 (TM-score is normalized by   91)

Aligned length=  87, RMSD=  2.52, TM-score=0.68780, ID=0.828

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     64.4096909452   0.5070101014  -0.0228654336   0.8616367733
 2     -1.0448182313   0.4864477233  -0.8176450973  -0.3079368562
 3      0.1858238387   0.7115541930   0.5752683434  -0.4034314856

(":" denotes the residue pairs of distance < 5.0 Angstrom)
SELEKAMVALIDVFHQYSG-REGDKHKLKKSELKELINNELSHFLEEIKEQ-EVVDKVMETLDNDGDGECDFQEFMAFVAMVTTACHEFFEHE--
::::::::::::::::::: ::: :::::::::::::::::::::: .::: ::::::::::::..:::::::::::::::::: ::: :::.
SELEKAMVALIDVFHQYSGREGD-KHKLKKSELKELINNELSHFLE-EIKEQEVVDKVMETLDEDGDGECDFQEFMAFVSMVTT-ACH-EFFEHE
"""

class TMAlign_Error( BiskitError ):
    pass

class TMAlign( Executor ):
    """
    Wrapper for TM-Align Structure alignment program.

    Usage::
    
       x = TMAlign( model1, model2 )   ## model1 and 2 are PDBModel instances
       result = x.run()
    
    Returns dictionary with the following keys:
       'matrix'  ... N.array of float, rotation translation matrix
       'rmsd'    ... float, rmsd calculated by TMAlign
       'length'  ... int, alignment length
       'score'   ... float, TMAlign score (normalized by length of target)
       'id'      ... float, percent identity of structure alignment
       'aln_model'  ... str, alignment string for model1 '-' indicates gap
       'aln_target' ... str, alignment string for target '-' indicates gap
       'aln_match'  ... N.array of boolean, CA within 5A of target in alignment
       
    Requires:
       * installation of TMAlign
    
    Configuration:
       * Biskit/data/defaults/exe_tmalign.dat (modify copy in ~/.biskit/)
    """
    
    re_score = re.compile( '.*length=\s*(?P<length>\d+)' +
                           '.*RMSD=\s*(?P<rmsd>\d+\.\d+)' +
                           '.*TM-score=\s*(?P<score>\d+\.\d+)' +
                           '.*ID=\s*(?P<id>\d\.\d+)'
                           )
    
    re_matrix = re.compile( '\s+([-]?\d+\.\d+)' )

    def __init__( self, model, target, **kw ):
        """
        @param model: model to be aligned
        @type  model: PDBModel
        @param target: target model of the alignment
        @type  target: PDBModel

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
        self.target = target
        
        assert isinstance( self.model, PDBModel), 'expecting PDBModel' 
        assert isinstance( self.target, PDBModel), 'expecting PDBModel'

        ## temporary pdb-file
        self.f_model  = tempfile.mktemp( '_tmalign_model.pdb')
        self.f_target = tempfile.mktemp( '_tmalign_target.pdb')
##        self.f_super  = tempfile.mktemp( '_tmalign_super.pdb' )
        
        Executor.__init__( self, 'tmalign',
                           args='%(f_model)s %(f_target)s' % self.__dict__,
                           catch_err=1, **kw )
        
        self.result = {}

    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.model.writePdb( self.f_model )
        self.target.writePdb(self.f_target)
    
    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_model )
            T.tryRemove( self.f_target )
 
    def isfailed( self ):
        if not self.output:
            return True
        if self.output[:10] == "\n ********":
            return False
        return True
    
    def __parse_score( self, lines ):
        l = lines.pop()
        r = self.re_score.match( l ).groupdict()
        r['length'] = int( r['length'])
        r['score']  = float( r['score'] )
        r['id']     = float( r['id'] )
        r['rmsd']   = float( r['rmsd'] )
        
        return r
    
    def __parse_matrix( self, lines ):
        lines.pop()
        lines.pop()
        
        r = []
        for i in range( 3 ):
            l = lines.pop()
            xyzt = self.re_matrix.findall( l )
            xyzt = [ float( x ) for x in xyzt ]
            r.append( xyzt )
        
        return {'matrix' : N.array( r ) }

    def __parse_aln( self, lines ):
        lines.pop()
        
        l1 = lines.pop()
        l2 = lines.pop()
        l3 = lines.pop()

        match = N.array( [ x == ':' for x in l2 ] ) ## all matching positions
        
        return {'aln_model':l1, 'aln_target':l3, 'aln_match': match }
            
 
    def parse_output( self, lines ):
        
        while lines:
            l = lines[-1]

            if l.startswith( 'Aligned' ):
                self.result.update( self.__parse_score( lines ) )

            elif l.startswith( ' -------- rotation' ):
                self.result.update( self.__parse_matrix( lines ) )

            elif l.startswith( '(":"' ):
                self.result.update( self.__parse_aln( lines ) )

            else:
                lines.pop()
            

    def finish( self ):
        """
        Called if program finished successfully (override postprocess otherwise)
        """
        if not self.exe.pipes:
            self.output = ''.join( open(self.f_out).readlines() )
        
        assert self.output, 'no output'
        
        lines = self.output.split( '\n' )
        lines.reverse()
        
        self.parse_output( lines )
        

if __name__ == '__main__':
    ## needs proper test case including failure result
    
    m1 = PDBModel( T.testRoot() + '/Mod/project/templates/modeller/1MQ1_A.pdb' )
    m2 = PDBModel( T.testRoot() + '/Mod/project/templates/modeller/1MWN_A.pdb' )
    
    x = TMAlign( m1 , m2 )
    r = x.run()
    


    
