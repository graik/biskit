##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg
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
Wrapper for structure alignment program TM-Align

Application: http://zhang.bioinformatics.ku.edu/TM-align/
Reference: Y. Zhang, J. Skolnick, TM-align: A protein structure alignment
algorithm based on TM-score , Nucleic Acids Research, 2005 33: 2302-2309
"""

import tempfile, re
import numpy as N

from Biskit import Executor, TemplateError
## import Biskit.settings as S
import Biskit.tools as T

class TMAlignError( Exception ):
    pass

class TMAlign( Executor ):
    """
    A TM-Align wrapper for sequence-independent structure alignments of 
    proteins.

    Usage
    =====

    >>> tm = TMAlign( model, refmodel )
    >>> result = tm.run()

    TMAlign takes two PDBModel instances as input and returns a dictionary 
    'result' with the following keys:

    'rt'    ... rotation /translation matrix to superposition model on refmodel
                (numpy.array of float)
    'score' ... TM-Align score (float)
    'rmsd'  ... RMSD calculated by TM-Align (float)
    'len'   ... length (in residues) of the structure alignment (float)
    'id'    ... sequence identity based on the structure alignment (float)

    @note: the alignment itself is currently not parsed in. Adapt parse_tmalign
           if you want to recover this information, too.
    @note: Command configuration: biskit/Biskit/data/defaults/exe_tmalign.dat
    """

    re_rt1 = re.compile(' 1\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)')
    re_rt2 = re.compile( ' 2\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)')
    re_rt3 = re.compile( ' 3\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)')

    re_info = re.compile( 'Aligned length=\s*([0-9]+), RMSD=\s*([0-9\.]+), TM-score=\s*([0-9\.]+), ID=\s*([0-9\.]+)')

    def __init__( self, model, refmodel, **kw ):
        """
        @param model: structure to be aligned to reference
        @type  model: PDBModel
        @param refmodel: reference structure
        @type  refmodel: PDBModel

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
        self.f_pdbin = tempfile.mktemp( '_tmalign_in.pdb' )
        self.f_pdbref= tempfile.mktemp( '_tmalign_ref.pdb' )

        Executor.__init__( self, 'tmalign', 
                           args= '%s %s' % (self.f_pdbin, self.f_pdbref),
                           **kw )

        self.refmodel = refmodel
        self.model = model


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.model.writePdb( self.f_pdbin )
        self.refmodel.writePdb( self.f_pdbref )

    def cleanup( self ):
        """
        Tidy up the mess you created.
        Does nothing. No temporary files are created.
        """        
        Executor.cleanup( self )
        if not self.debug:
            T.tryRemove( self.f_pdbin )
            T.tryRemove( self.f_pdbref)

    def __translate_rt( self, rt ):
        """
        TM-Align reports translation vector in first row -> push it to last row
        """
        rt = rt[:, (1,2,3,0)]
        return rt

    def parse_tmalign( self, output ):
        """
        Parse TM-Align output
        @param output: STDOUT result of TM-Align run
        @type  output: [str]

        @return: rotation/translation matrix
        @rtype: N.array

        @raise TMAlignError: if no result
        """
        if  output.count('\n') < 7:
            raise TMAlignError, 'no TM-Align result'

        r = {}

        try:
            rt1 = self.re_rt1.findall( output )[0]
            rt2 = self.re_rt2.findall( output )[0]
            rt3 = self.re_rt3.findall( output )[0]            
        except IndexError, why:
            raise TMAlignError(
                'Could not find rotation matrix in TMAlign output')

        try:
            info = self.re_info.findall( output )[0]
            info = [ float( i ) for i in info ]
        except IndexError, why:
            raise TMAlignError(
                'Could not find score values in TMAlign output')

        r.update( dict( zip( ['len', 'rmsd', 'score', 'id'], info ) ) )        
        r['rt'] = self.__translate_rt( N.array( (rt1, rt2, rt3), float ) )

        return r

    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.error != '' or \
               'lately' in self.output

    def fail( self ):
        """
        Overrides Executor method. Called when execution fails.
        """
        s = 'TMAlign failed. Please check the program output in the '+\
          'field `output` of this TMAlign instance (e.g. `print x.output`)!'
        self.log.add( s )

        raise TMAlignError, s

    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_tmalign( self.output )


    def applyTransformation( self, model=None ):
        """
        Apply transformation from structure alignment to given model.
        @param model: PDBModel, external model (default: input model)
        @return: transformed PDBModel
        @rtype: PDBModel
        """
        model = model or self.model
        assert self.result, 'call TMAlign.run() first!'

        r = model.transform( self.result['rt'] )
        r.info['tm_score'] = self.result['score']
        r.info['tm_id'] = self.result['id']
        r.info['tm_len'] = self.result['len']
        r.info['tm_rmsd']= self.result['rmsd']

        return r



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    def test_tmalign( self ):
        """TMAlign test"""
        from Biskit import PDBModel

        if self.local: print 'Loading PDB...'

        self.m1 = PDBModel( T.testRoot( 'tmalign/1huy_citrine.model' ) )
        self.m2 = PDBModel( T.testRoot('tmalign/1zgp_dsred_dimer.model' ) )


        if self.local: print 'Starting TMAlign'
        self.x = TMAlign( self.m1, self.m2, debug=self.DEBUG,
                          verbose=self.local )

        if self.local:
            print 'Running'

        self.r = self.x.run()

        if self.local:
            print "Result: "
            for key, value in self.r.items():
                print '\t', key, ':\t', value

        self.assertEqual( self.r['rmsd'], 1.76 )

    def test_tmalignTransform( self ):
        """TMAlign.applyTransformation test"""
        m = T.load( T.testRoot( 'tmalign/1huy_citrine.model' ) )
        ref = T.load( T.testRoot( 'tmalign/1zgp_dsred_dimer.model' ) )
        ref = ref.takeChains( [0] )

        tm = TMAlign( m, ref )
        tm.run()

        self.maligned = tm.applyTransformation()

        diff = self.maligned.centerOfMass() - ref.centerOfMass()

        if self.VERBOSITY > 2 or self.local:
            print 'center of mass deviation: \n%r' % diff
            self.maligned.concat( ref ).plot()

        self.assert_( N.all( N.absolute(diff) < 1 ),
                      'superposition failed: \n%r' % diff)



if __name__ == '__main__':

    BT.localTest()
