##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$
"""Get binding energy from fold_X"""

from Biskit import Executor, TemplateError
import Biskit.settings as S
import Biskit.tools as T

import re, string, tempfile

import Biskit.molUtils as molUtils


class Fold_XError( Exception ):
    pass


class Fold_X( Executor ):
    """
    Run fold_X with given PDBModel, return dictionary with energy terms

    Example usage:
        x = Fold_X( model, verbose=1 )
        result = x.run()

    Energy terms (kcal/mol):
     DDG  - deltaG 
     vwcl - Van der Waals clash
     mc	  - Main chain entropy
     sc	  - Side chain entropy
     vw	  - Van der Waals interaction
     hb	  - Hydrogen bond energy
     hyd  - Hydrophobic desolvation
     ene  - Total energy
     wtbr - Water Bridge Energy
     pol  - Polar desolvation
     el	  - Electrostatic energy
    
    See http://foldx.embl.de for more info.
   
    August 23 - 2005:
    Documentation is sparse here as a new wersion of FoldX is in the
    works and the old version is unavaliable for download. The new
    version seems to have quite som changes to the interface.
    """

    def __init__(self, model, **kw ):
        """
        model - PDBModel, reference

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        """
        temp_pdb = tempfile.mktemp('_foldx_.pdb')
        
        Executor.__init__( self, 'fold_X', args='-b -u %s'%temp_pdb,
                           f_in=temp_pdb, **kw )
        
        self.model = model.clone( deepcopy=1 )
            
        ## fold-X-allowed atoms for each res in standard order
        self.aminoAcidDict = molUtils.aaAtoms
        for k in self.aminoAcidDict:
            self.aminoAcidDict[ k ] += ['HN']


    def __prepareModel( self, model, f_pdb_out ):
        """
        Prepare model for fold_X:
        - remove all hydrogens but H
        - rename H to HN
        - rename H1 to HN (remove H2 and H3)
        (- terminal oxygens should be OXT )
        - sort atoms according to self.aminoAcidDict
        """
                                        
        def __cmpAtoms( a1, a2 ):
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
 #       model = model.take( range(model.lenAtoms()), deepcopy=1  )

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
        """Overrides Executor method."""

        self.__prepareModel( self.model, self.f_in )


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_in )
            T.tryRemove( 'foldx_results' )

    def parse_foldx( self, output ):
        """
        Extract energies from output.
        -> dictionary with energy terma as keys
        """
        energy = {}
        
        ## extract energy terms from output
        l= open(self.f_out).readlines()

        if len(l) == 0:
            raise Fold_XError, 'ERROR: fold_X Segmentation fault (core dumped)'
        
        if re.match('^.+nbhb.+vw', l[-2]):
            E = string.split( l[-2])[2:]
            E = [ string.split( e, ':') for e in E ]
            ## add energies to dictionary
            for i in range( len(E) ):
                ## change nan to 0.0
                if E[i][1] == 'nan':
                    E[i][1] = 0.0
                energy[str(E[i][0])] = float( E[i][1] )
        else:
            print 'No fold_X energy to extract'
            energy = 0
            
        return energy


    def isFailed( self ):
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_foldx( self.output )
        

##########
## test ##
    
if __name__ == '__main__':

    from Biskit import PDBModel
    import Biskit.tools as T
    import glob

    print "Loading PDB..."

    f = glob.glob( T.testRoot()+'/rec_pc2_00/pdb/*_1_*pdb.gz' )[1]
    m = PDBModel(f)
    m = m.compress( m.maskProtein() )

    print "Starting fold_X"

    x = Fold_X( m, debug=0, verbose=1 )

    print "Running"
    r = x.run()

    print "Result: ", r




