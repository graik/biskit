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
Map (Amber) atomic partial charges on a structure.
"""
import numpy as N

from biskit import PDBModel, StdLog, AmberResidueLibrary
import biskit.mathUtils as MA


class AtomCharger( object ):
    """
    Assign atomic partial charges to each atom of a structure. Charges are
    taken from a collection of L{Biskit.AmberResidueType} instances.    
    L{Biskit.AmberResidueLibrary} bundles the content of several amber
    topology files and identifies matching residue types by atom content and
    (if that fails) by name. The original topology files are bundled with 
    Biskit in Biskit/data/amber/residues. 

    Residues in the input PDB are first matched against Amber residues with
    exactly the same atom content (same set of atom names). That means, the
    name of the residue can differ (as is the case for e.g. CALA, CTHR or
    NALA, etc.).

    Residues for which there was no exact atom match in the
    AmberResidueLibrary, will be matched against a residue type of the same
    name. Missing atoms will be ignored which will lead to an incomplete
    charge. However, if there is any unknown atoms (atoms that cannot be
    matched to any atom in the reference residue type) charge 0 is assigned to
    every atom of the *whole* residue. 

    Any missmatches are reported as a warning to the log instance if
    verbose==True.

    Use:
    ====

    >>> m = PDBModel( 'my.pdb' )
    >>> ac = AtomCharger( verbose=True )
    >>> ac.charge( m )

    Partial charges are then assigned to the atom profile m['partial_charge']
    (one value for each atom). That means the charge (float value) of the
    third atom can be accessed as:

    >>> m['partial_charge'][2]

    And the numpy array of all charges can be summed into the total charge of 
    the molecule:

    >>> q = numpy.sum( m['partial_charge'] )

    We chose 'partial_charge' because the standard PDB data model already 
    contains a record 'charge', which is set to a string and may contain ugly
    things like '+2' or '-1' but cannot hold float values.

    Note:
    =====

    The default Amber residue topologies are all-atom including hydrogens.
    Most likely, you will need to add hydrogens to your PDBModel before
    passing it to AtomCharger. Use L{Biskit.Reduce} for this purpose. Example:

    >>> m = PDBModel( '3tgi' )
    >>> r = Reduce( m )
    >>> m = r.run()

    >>> ac = AtomCharger()
    >>> ac.charge( m )

    ...This will run the reduce program to optimize hydrogen bonding networks
    and add hydrogen atoms, and only then call AtomCharger to match each
    residue to the Amber residue description with exactly the same atom
    content.

    Customization:
    ==============

    If you want to transfer charges from a different set of topologies, 
    create your own L{AmberResidueLibrary} instance and pass it to the 
    constructor. For example:

    >>> topofiles = ['uni_amino03', 'uni_aminoct03', 'uni_aminont03']
    >>> reslib = AmberResidueLibrary( topofiles, verbose=True )
    >>> ac = AtomCharger( reslib )

    This will match residues against the united atom topology files that are
    also found in Biskit/data/amber/residues. You can specify full file
    names to point to your own prep files.

    See also the testing code for a more elaborate example of using
    AtomCharger.

    @see L{AmberResidueType}
    @see L{AmberResidueLibrary}
    @see L{Reduce}
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


    def lookupResidue(self, model):
        """
        Identify best-matching ResidueType. First try by atom content then by
        name.
        @param model: single residue PDBModel
        @type  model: PDBModel
        @return: residue type with name adjusted to the model's residue name
        @rtype : AmberResidueType OR None
        """
        resname = model['residue_name'][0].upper()

        refres = self.reslib.byAtoms( model )

        if refres is None:
##            if self.verbose:
##                self.log.add('Warning: no exact atom match for residue '+\
##                             '%3s %3i'% (resname, model['residue_number'][0])+\
##                             '\nNow trying by residue name...')
            try:
                refres = self.reslib.byName( resname )
            except KeyError as what:
                return None

        refres = refres.clone()
        refres['residue_name'] = len(refres) * [resname]

        return refres

    def __resinfo( self, res, i ):
        r = '#%-3i (%3s %3i %1s)' % (i, res['residue_name'][0], 
                                     res['residue_number'][0], res['chain_id'][0])
        return r

    def charge(self, model):
        """
        Assign atom charges to given model
        @param model: input structure
        @type  model: PDBModel
        """
        q = N.array( [], float )

        for j,res in enumerate(model.resModels()):
            refres = self.lookupResidue( res )

            if refres is None:
                if self.verbose:
                    self.log.add('Warning: no residue type for '+\
                                 self.__resinfo(res,j) +\
                                 '\nWhole residue will be set to charge 0.')
                refres = res
                refres['partial_charge'] = N.zeros( len(refres) )

            iref, i = refres.compareAtoms( res )

            qres = N.take( refres['partial_charge'], iref )

            if len(iref) < len(refres):
                if self.verbose:
                    self.log.add('Warning: %i atoms missing from residue %s:'%\
                                 (len(refres)-len(iref), self.__resinfo(res,j)) )

                    missing = MA.difference( list(range(len(refres))), iref )
                    self.log.add('\t'+str(N.take( refres['name'], missing )) )
                    self.log.add(
                        '\tIgnored. Residue will have incomplete charge.')

            if len(i) < len( res ):
                if self.verbose:
                    self.log.add('Warning: %i unknown  atoms in residue %s:' %\
                                 (len(res)-len(i), self.__resinfo(res,j)) )

                    missing = MA.difference( list(range(len(res))), i )
                    self.log.add('\t'+str(N.take( res['name'], missing )) )
                    self.log.add('\tGiving up. Whole residue will be set to charge 0.')

                qres = N.zeros( len(res) )

##            assert len(qres) == len(res), 'missmatch of atom number'

            q = N.concatenate( (q, qres ) )

        assert len(q) == len(model), 'AtomCharger: missing charge records'
        model['partial_charge'] = q

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    def test_atomcharger( self ):
        """AtomCharger test"""
        import biskit.tools as T

        if self.local: self.log.add('\nLoading PDB...')

        self.m1 = PDBModel( T.testRoot('delphi/1A19_reduced.model'))
        self.m2 = PDBModel( T.testRoot('delphi/complex_reduced.model'))

        if self.local:
            self.log.add('\nSetup Residue Library\n')

        ac = AtomCharger(log=self.log, verbose=self.local)

        if self.local:
            self.log.add('match residues to Amber topology')

        ac.charge( self.m1 )
        ac.charge( self.m2 )

        self.assertAlmostEqual( N.sum(self.m1['partial_charge']), -6, 2 )
        self.assertAlmostEqual( N.sum(self.m2['partial_charge']), -4, 2 )
        self.assertTrue(N.all(self.m1['partial_charge'] != 0),'unmatched atoms 1')
        self.assertTrue(N.all(self.m2['partial_charge'] != 0),'unmatched atoms 2')

        if self.local:
            self.log.add('\nNow test handling of atom miss-matches:\n')

        self.m3 = PDBModel(self.m1.clone())
        self.m3.remove( [0,3,100,101,102,200] )
        ac.charge( self.m3 )

        self.assertAlmostEqual( N.sum(self.m3['partial_charge']),-8.21, 2)


if __name__ == '__main__':

    BT.localTest(debug=False)


