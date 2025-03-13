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
Wrapper for structure alignment program TM-Align

Application: http://zhang.bioinformatics.ku.edu/TM-align/
Reference: Y. Zhang, J. Skolnick, TM-align: A protein structure alignment
algorithm based on TM-score , Nucleic Acids Research, 2005 33: 2302-2309
"""

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import tempfile, re
import numpy as N

from .executor import Executor, TemplateError
## import Biskit.settings as S
import biskit.tools as T

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

    - `rt`  ... rotation /translation matrix to superposition model on refmodel\
                (numpy.array of float), see `PDBModel.transform`
    - `score` ... TM-Align score (float)
    - `rmsd`  ... RMSD calculated by TM-Align (float)
    - `len`   ... length (in residues) of the structure alignment (float)
    - `id`    ... sequence identity based on the structure alignment (float)

    Optionally, TMAlign can also apply the transformation to the input
    (or any other) structure:

    >>> model_transformed = tm.applyTransformation(model)

    The TM result values (score, rmsd, etc) are added to `model_transformed.info`. 

    .. note:: 
    
        The alignment itself is currently not parsed in. Adapt parse_tmalign
        if you want to recover this information, too.
        
        Command configuration: `biskit/Biskit/data/defaults/exe_tmalign.dat`
    """
    _f = r'\s+([-0-9]+\.[0-9]+)'
    re_rt1 = re.compile(r'\s*1' + _f + _f + _f + _f)
    re_rt2 = re.compile(r'\s*2' + _f + _f + _f + _f)
    re_rt3 = re.compile(r'\s*3' + _f + _f + _f + _f)

    re_info = re.compile( r'Aligned length=\s*(?P<len>[0-9]+), RMSD=\s*(?P<rmsd>[0-9\.]+), Seq_ID=.+=\s*(?P<id>[0-9\.]+)')
    re_score= re.compile( r'TM-score=\s*(?P<score>[0-9\.]+)\s+\(if normalized by length of Chain_2\)' ) 

    def __init__( self, model, refmodel, **kw ):
        """
        :param model: structure to be aligned to reference
        :type  model: PDBModel
        :param refmodel: reference structure
        :type  refmodel: PDBModel

        :param kw: additional key=value parameters are passed on to \
            `Executor.__init__`. For example:
            ::
                debug    - 0|1, keep all temporary files (default: 0)
                verbose  - 0|1, print progress messages to log (log != STDOUT)
                node     - str, host for calculation (None->local) NOT TESTED
                                (default: None)
                nice     - int, nice level (default: 0)
                log      - biskit.LogFile, program log (None->STOUT) (default: None)
        """
        self.f_pdbin = tempfile.mktemp( '_tmalign_in.pdb' )
        self.f_pdbref= tempfile.mktemp( '_tmalign_ref.pdb' )
        self.f_matrix= tempfile.mktemp( '_tmalign_matrix.out' )

        Executor.__init__( self, 'tmalign', 
                           args= '%s %s -m %s' % (self.f_pdbin, self.f_pdbref, self.f_matrix),
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
            T.tryRemove( self.f_matrix )

    def __translate_rt( self, rt ):
        """
        TM-Align reports translation vector in first row -> push it to last row
        """
        rt = rt[:, (1,2,3,0)]
        return rt
    
    def parse_matrix(self):
        try:
            lines = open(self.f_matrix, 'r').readlines()
            s = ''.join(lines)
            
            rt1 = self.re_rt1.findall( s )[0]
            rt2 = self.re_rt2.findall( s )[0]
            rt3 = self.re_rt3.findall( s )[0]            

            r = self.__translate_rt( N.array( (rt1, rt2, rt3), float ) )            
            return r
        
        except IOError as why:
            raise TMAlignError('Cannot find matrix output file.')

        except IndexError as why:
            raise TMAlignError(
                'Could not find rotation matrix in TMAlign output')


    def parse_tmalign( self, output ):
        """
        Parse TM-Align output
        
        :param output: STDOUT result of TM-Align run
        :type  output: [str]

        :return: rotation/translation matrix
        :rtype: N.array

        :raise TMAlignError: if no result
        """
        if  output.count('\n') < 7:
            raise TMAlignError('no TM-Align result')

        r = {}
        try:
            r = self.re_info.search( output ).groupdict()
            r.update( self.re_score.search(output).groupdict() )
            r = { k : float(v) for k, v in r.items() }
        except AttributeError as why:
            raise TMAlignError(
                'Could not find score and rmsd values in TMAlign output')

        r['rt'] = self.parse_matrix()

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

        raise TMAlignError(s)

    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_tmalign( self.output )


    def applyTransformation( self, model=None ):
        """
        Apply transformation from structure alignment to given model.
        
        :param model: external model (default: input model)
        :type model: PDBModel
        :return: transformed PDBModel
        :rtype: PDBModel
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
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    def test_tmalign( self ):
        """TMAlign test"""
        from biskit import PDBModel

        if self.local: print('Loading PDB...')

        self.m1 = PDBModel( T.testRoot( 'tmalign/1huy_citrine.model' ) )
        self.m2 = PDBModel( T.testRoot('tmalign/1zgp_dsred_dimer.model' ) )


        if self.local: print('Starting TMAlign')
        self.x = TMAlign( self.m1, self.m2, debug=self.DEBUG,
                          verbose=self.local )

        if self.local:
            print('Running')

        self.r = self.x.run()

        if self.local:
            print("Result: ")
            for key, value in self.r.items():
                print('\t', key, ':\t', value)

        self.assertEqual( self.r['rmsd'], 1.76 )

    def test_tmalignTransform( self ):
        """TMAlign.applyTransformation test"""
        m = T.load( T.testRoot( 'tmalign/1huy_citrine.model' ) )
        ref = T.load( T.testRoot( 'tmalign/1zgp_dsred_dimer.model' ) )
        ref = ref.takeChains( [0] )

        tm = TMAlign( m, ref, debug=self.DEBUG, verbose=self.local )
        tm.run()

        self.maligned = tm.applyTransformation()

        diff = self.maligned.centerOfMass() - ref.centerOfMass()

        if self.VERBOSITY > 2 or self.local:
            print('center of mass deviation: \n%r' % diff)
            self.maligned.concat( ref ).plot()

        self.assertTrue( N.all( N.absolute(diff) < 1 ),
                      'superposition failed: \n%r' % diff)



if __name__ == '__main__':

    BT.localTest(debug=False)
