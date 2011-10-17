##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2011 Raik Gruenberg
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
Map (Amber) atomic partial charges on a structure.
"""
import numpy as N

from Biskit import PDBModel, StdLog
from Biskit import AmberResidueLibrary
from Biskit import Reduce


class AtomCharger( object ):
    """
    Assign atomic partial charges to each atom of a structure. Charges are
    taken from a dictionary of L{Biskit.AmberResidueType} instances. This
    dictionary can be created with the AmberResidue class from standard Amber
    topology files (bundled with Biskit in Biskit/data/amber/residues). The
    content of several Amber topology files is bundled up into a
    L{Biskit.AmberResidueLibrary} instance and passed on to AtomCharger.
    
    Residues in the input PDB are matched against Amber residues with exactly 
    the same atom content (same set of atom names). That means, the name of 
    the residue can differ (as is the case for e.g. CALA, CTHR or NALA, etc.).

    Residues for which there was no atom match in the AmberResidueLibrary, get 
    charge 0 assigned to every atom. A warning is reported to the log instance
    if verbose==True.
    
    Note, the default Amber residue topologies are all-atom including
    hydrogens. Most likely, you will need to add hydrogens to your PDBModel
    before passing it to AtomCharger. Use L{Biskit.Reduce} for this purpose.
    Example:
    
    >>> m = PDBModel( '3tgi' )
    >>> r = Reduce( m, autocap=True )
    >>> m = r.run()
    
    >>> ac = AtomCharger()
    >>> ac.charge( m )
    
    ...This will first autodetect and "cap" any protein chain breaks or
    missing N- and C-terminals with ACE and NME residues, then run the reduce
    program to optimize hydrogen bonding networks and add hydrogen atoms, and
    only then call AtomCharger to match each residue to the Amber residue
    description with exactly the same atom content. Partial charges are then
    assigned to the atom profile m['charges'] (one value for each atom). That
    means the charge (float value) of the third atom can be accessed as:
    
    >>> m['partial_charge'][2]
    
    And the numpy array of all charges can be summed into the total charge of 
    the molecule:
    
    >>> q = numpy.sum( m['partial_charge'] )
    
    See also the testing code for a more elaborate example of using
    AtomCharger.
    """
    
    def __init__(self, reslibrary=None, log=None, verbose=False ):
        """
        @param reslibrary: Collection of AmberResidueType instances to be 
                           matched against the residues of input models.
        @type  reslibrary: AmberResidueLibrary
        
        @param log: log (default STDOUT)
        @type  log: Biskit.LogFile
        
        @param verbose: give status updates and warnings to log
        @type  verbose: bool
        """
        self.reslib = reslibrary or AmberResidueLibrary(log=log, 
                                                        verbose=verbose)
        self.verbose = verbose
        self.log = log or StdLog()
    

    def charge(self, model):
        """
        Assign atom charges to given model
        @param model: input structure
        @type  model: PDBModel
        """
        q = N.array( [], float )
        failed = []
        
        for res in model.resModels():
            refres = self.reslib.byAtoms( res )

            if refres is None:
                failed += [ res ]
                if self.verbose:
                    self.log.add('Warning: could not find residue type for '+\
                                 '%s %i' % (res['residue_name'][0],
                                            res['residue_number'][0]) )
                refres = res
                refres['partial_charge'] = N.zeros( len(refres) )
            
            refres = refres.clone()
            refres['residue_name'] = len(refres) * [ res['residue_name'][0] ]
            iref, i = refres.compareAtoms( res )
            
            qres = N.take( refres['partial_charge'], iref )
            
##            assert len(qres) == len(res), 'missmatch of atom number'
            
            q = N.concatenate( (q, qres ) )
        
        assert len(q) == len(model), 'AtomCharger: missing charge records'
        model['partial_charge'] = q

#############
##  TESTING        
#############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """Test class"""

    def test_reduce( self ):
        """AtomCharger test"""
        if self.local: self.log.add('\nLoading PDB...')

        self.m1 = PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        self.m2 = T.load( T.testRoot( 'com/ref.complex' ) )
        self.m2 = self.m2.model()

        if self.local: self.log.add('\nRunning Reduce...')
        self.x = Reduce( self.m1, debug=self.DEBUG, verbose=self.local,
                         log=self.log,
                         autocap=True )
        self.m1 = self.x.run()

        if self.local:
            self.log.add('\nReduce protein complex')
        self.x = Reduce( self.m2, debug=self.DEBUG, verbose=self.local,
                         log=self.log,
                         autocap=True )
        self.m2 = self.x.run()
        
        if self.local:
            self.log.add('\nSetup Residue Library\n')
        
        ac = AtomCharger()
        
        if self.local:
            self.log.add('match residues to Amber topology')
        
        ac.charge( self.m1 )
        ac.charge( self.m2 )

        self.assertAlmostEqual( N.sum(self.m1['partial_charge']), -6, 2 )
        self.assertAlmostEqual( N.sum(self.m2['partial_charge']), -4, 2 )
        self.assert_( N.all(self.m1['partial_charge'] != 0), 'unmatched atoms 1' )
        self.assert_( N.all(self.m2['partial_charge'] != 0), 'unmatched atoms 2' )
        

if __name__ == '__main__':

    BT.localTest(debug=False)


