##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
Get binding energy from fold_X
"""

from Biskit import Executor, TemplateError
from Biskit import BiskitError

import Biskit.tools as T

import re, string, tempfile

import Biskit.molUtils as molUtils


class Fold_XError( BiskitError ):
    pass


class Fold_X( Executor ):
    """
    Run fold_X with given PDBModel
    ==============================
      Returns a dictionary with energy terms.

      Example usage
      -------------
         >>> x = Fold_X( model, verbose=1 )
         >>> result = x.run()

      Energy terms (kcal/mol)
      -----------------------
      This is the order of the energy terms in the output file. There
      are two additional values that are unlabeled in the output
      (version 2.5.2). The first is the entropy cost of making a complex.
      This is zero if not 'AnalyseComplex' is run. The final value are
      the number of residues.

      We'll try to keep the energy labels of Fold-X version 1 if possible.
      ::
        ene    - total energy
        bb_hb  - Backbone Hbond
        sc_hb  - Sidechain Hbond
        vw     - Van der Waals
        el     - Electrostatics
        sol_p  - Solvation Polar  
        sol_h  - Solvation Hydrophobic
        vwcl   - Van de Waals clashes
        sc     - entropy side chain
        mc     - entropy main chain
        sloop  - sloop_entropy
        mloop  - mloop_entropy
        cis    - cis_bond
        tcl    - torsional clash   
        bbcl   - backbone clash
        dip    - helix dipole
        wtbr   - water bridge
        disu   - disulfide
        el_kon - electrostatic kon
        p_cov  - partial covalent bonds
             ( - complex entropy )
             ( - number of residues )        

      Reference
      ---------
        U{http://foldx.embl.de} for more info

    @note: September 11 - 2006: 
      Tested with Fold-X 2.5.2
    """

    def __init__(self, model, **kw ):
        """
        @param model: reference PDBModel
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
        self.temp_pdb     = tempfile.mktemp('_foldx_.pdb')
        self.temp_command = tempfile.mktemp('_foldx_.command')
        self.temp_option  = tempfile.mktemp('_foldx_.option')
        self.temp_result  = tempfile.mktemp('_foldx_.result')
        self.temp_runlog  = tempfile.mktemp('_foldx_.log')
        self.temp_errlog  = tempfile.mktemp('_foldx_.err')

        Executor.__init__( self, 'fold_X', args='-manual %s %s %s'\
                           %(self.temp_pdb, self.temp_option,
                             self.temp_command), **kw )

        self.model = model.clone()

        ## fold-X-allowed atoms for each res in standard order
        self.aminoAcidDict = molUtils.aaAtoms
        for k in self.aminoAcidDict:
            if 'HN' not in self.aminoAcidDict[ k ]:
                self.aminoAcidDict[ k ] += ['HN']


    def __prepareModel( self, model, f_pdb_out ):
        """
        Prepare model for fold_X and write to file::
          - remove all hydrogens but H
          - rename H to HN
          - rename H1 to HN (remove H2 and H3)
          - ( terminal oxygens should be OXT )
          - sort atoms according to self.aminoAcidDict

        @param model: reference PDBModel
        @type  model: PDBModel
        @param f_pdb_out: remaned pdb file name
        @type  f_pdb_out: str

        @raise Fold_XError: if preparation fail
        """

        def __cmpAtoms( a1, a2 ):
            """
            compare atoms.

            @param a1: atom name
            @type  a1: str        
            @param a2: atom name
            @type  a2: str
            """
            res = a1['residue_name']
            target = self.aminoAcidDict[ res ]
            try:
                return cmp(target.index( a1['name'] ),
                           target.index( a2['name'] ))
            except ValueError, why:
                s = "Unknown atom for %s %i: %s or %s" % \
                  (res, a1['residue_number'], a1['name'], a2['name'] )
                raise Fold_XError( s )

        ## make a copy
    #       model = model.take( range(model.lenAtoms()) )

        ## mask for all heavy atoms, H and H3
        heavy_mask = model.maskHeavy()
    #       HN_mask = model.mask( lambda a: a['name'] == 'H' )
    #       H3_mask = model.mask( lambda a: a['name'] == 'H3' )

            ## rename H and H3 -> HN
    #       atm_dic = model.getAtoms()
    #       for i in N.nonzero( HN_mask + H3_mask ):
    #           atm_dic[i]['name'] = 'HN'

        ## remove all none backbone hydrogens
        keep_mask = heavy_mask #+ HN_mask + H3_mask
        model = model.compress( keep_mask )

        ## consecutive residue numbering
        model.renumberResidues()

        ## sort atoms
        model = model.sort( model.argsort( __cmpAtoms ) )

        model.writePdb( f_pdb_out )


    def prepare( self ):
        """
        Write a Fold-X compatible pdb file to disc.

        @note: Overrides Executor method.
        """
        self.__prepareModel( self.model, self.temp_pdb )

        f_com = open( self.temp_command, 'w')
        f_com.writelines(['<TITLE>FOLDX_commandfile;\n',
                          '<Stability>%s;\n'%self.temp_result])
        f_com.close()

        f_opt = open( self.temp_option, 'w')
        f_opt.writelines(['<TITLE> FOLDX_optionfile;\n',
                          '<logfile_name> %s;'%self.temp_runlog,
                          '<errorfile_name> %s;'%self.temp_errlog])
        f_opt.close()


    def cleanup( self ):
        """
        Remove files created for and by the calculation. 
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.temp_pdb )
            T.tryRemove( self.temp_command )
            T.tryRemove( self.temp_option )
            T.tryRemove( self.temp_result )
            T.tryRemove( self.temp_runlog )
            T.tryRemove( self.temp_errlog )
            ## Fold-X writes a file called "runlog.txt"
            ## to local directory. Try to remove it.
            T.tryRemove( 'runlog.txt' )
            ## and even though the error log is supposed
            ## to be written to self.temp_errlog, I get a
            ## 'errorfile.txt' in the local directory. Remove.
            T.tryRemove( 'errorfile.txt' )


    def parse_foldx( self):
        """
        Extract energies from output.

        @return: dictionary with energy terms as keys
        @rtype: dict
        """
        energy = {}

        ## extract energy terms from output
        l= open(self.temp_result).readlines()

        if len(l) == 0:
            raise Fold_XError, 'ERROR: fold_X Segmentation fault'

        if re.match('^%s'%self.temp_pdb, l[-1]):

            ## get the energies
            E = string.split( l[-1])[1:]

            ## energy labels can't be parsed from file, so
            ## we'll have to use this list (and lets hope
            ## they don't change the order of the items).
            keys = [ 'ene',  'bb_hb', 'sc_hb',  'vw',  
                     'el',   'sol_p', 'sol_h',  'vwcl',
                     'sc',   'mc',    'sloop',  'mloop',
                     'cis',  'tcl',   'bbcl',   'dip',
                     'wtbr', 'disu',  'el_kon', 'p_cov' ]

            ## add energies to dictionary
            for i in range( len( keys ) ):
                energy[ keys[i] ] = float( E[i] )
        else:
            print 'No fold_X energy to extract'
            energy = 0

        return energy


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
        self.result = self.parse_foldx()

#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Fold_X test"""

    TAGS = [ BT.EXE ]

    def test_Fold_X( self):
        """Fold_X test"""
        from Biskit import PDBModel

        ## Loading PDB...
        self.m = PDBModel( T.testRoot() + '/rec/1A2P.pdb' )
        self.m = self.m.compress( self.m.maskProtein() )

        ## Starting fold_X
        self.x = Fold_X( self.m, debug=self.DEBUG,
                         verbose=(self.VERBOSITY>2) )

        ## Running
        self.r = self.x.run()

        if self.local:
            print "Result: ", self.r

        for k, v in self.r.items():
            self.assertAlmostEqual( v, self.EXPECTED[k], 6 )


    EXPECTED = {'el': -13.766400000000001, 'wtbr': -4.8059700000000003, 'ene': -18.475000000000001, 'mc': 160.28800000000001, 'sloop': 0.0, 'dip': 0.00177626, 'sol_p': 167.71100000000001, 'disu': 0.0, 'tcl': 0.72696700000000003, 'cis': 0.0, 'p_cov': 0.0, 'sol_h': -134.613, 'bb_hb': -87.362499999999997, 'sc_hb': -48.350000000000001, 'vw': -116.67100000000001, 'sc': 58.089300000000001, 'el_kon': 0.0, 'mloop': 0.0, 'vwcl': 0.27728599999999998, 'bbcl': 0.37019200000000002}



if __name__ == '__main__':

    BT.localTest()


