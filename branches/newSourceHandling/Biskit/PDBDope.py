## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## $Revision$
## last $Author$
## last $Date$
"""
Calculate and add various properties to PDBModel
"""

import numpy.oldnumeric as N
import Biskit.tools as T

from Biskit.WhatIf import WhatIf 
from Biskit.Hmmer import Hmmer
from Biskit.DSSP import Dssp
from Biskit.Fold_X import Fold_X
from Biskit.SurfaceRacer import SurfaceRacer


class PDBDope:
    """
    Decorate a PDBModel with calculated properties (profiles)
    """

    def __init__( self, model ):
        """
        @param model: model to dope
        @type  model: PDBModel
        """
        self.m = model

    def version( self ):
        """
        @return: version of class
        @rtype: str
        """
        return 'PDBDope $Revision$'

    def model( self ):
        """
        @return: model
        @rtype: PDBModel
        """
        return self.m


    def addASA( self ):
        """
        Add profiles of Accessible Surface Area: 'relASA', 'ASA_total',
        'ASA_sc', 'ASA_bb'. See L{Biskit.WhatIf}

        @note: Using WhatIf to calculate relative accessabilities is not
               nessesary any more. SurfaceRacer now also adds a profile
               'relAS' (conataining relative solvent accessible surf) and
               'relMS' (relative molecular surface).
        
        @raise ProfileError: if WhatIf-returned atom/residue lists don't match.
                             Usually that means, WhatIf didn't recognize some
                             residue name
        """
        w = WhatIf( self.m )

        atomRelAcc, resASA, resMask = w.run()

##         normalAtoms = N.logical_not( N.logical_or(self.m.maskHetatm(),
##                                                   self.m.maskSolvent() ) )

        normalAtoms = self.m.maskProtein( standard=1 )

        normalRes = self.m.atom2resMask( normalAtoms )

        self.m.atoms.set( 'relASA', atomRelAcc, ## normalAtoms, 0,
                               comment='relative accessible surface area in %',
                               version= T.dateString() + ' ' + self.version() )

        self.m.residues.set( 'ASA_total', resASA[:,0], normalRes, 0,
                              comment='accessible surface area in A^2',
                              version= T.dateString() + ' ' + self.version() )

        self.m.residues.set( 'ASA_sc', resASA[:,1], normalRes, 0,
                           comment='side chain accessible surface area in A^2',
                           version= T.dateString() + ' ' + self.version() )

        self.m.residues.set( 'ASA_bb', resASA[:,2], normalRes, 0,
                           comment='back bone accessible surface area in A^2',
                           version= T.dateString() + ' ' + self.version() )


    def addSurfaceMask( self, pname='relAS' ):
        """
        Adds a surface mask profie that contains atoms with > 40% exposure
        compared to a random coil state.

        @param pname: name of relative profile to use
                      (Whatif-relASA OR SurfaceRacer - relAS)
                      (default: relAS)
        @type  pname: str
        """
        r = self.m.profile2mask( pname, cutoff_min=40 )
        self.m.residues.set( 'surfMask',  self.m.atom2resMask(r),
                              comment='residues with any atom > 40% exposed',
                              version= T.dateString() + ' ' + self.version() )


    def addSecondaryStructure( self ):
        """
        Adds a residue profile with the secondary structure as
        calculated by the DSSP program.

        Profile code::
          B = residue in isolated beta-bridge
          E = extended strand, participates in beta ladder
          G = 3-helix (3/10 helix)
          I = 5 helix (pi helix)
          T = hydrogen bonded turn
          S = bend
          . = loop or irregular

        @raise ExeConfigError: if external application is missing
        """
        dssp = Dssp( self.m )
        ss = dssp.run()
        
        self.m.residues.set( 'secondary',  ss,
                              comment='secondary structure from DSSP',
                              version= T.dateString() + ' ' + self.version() )
        

    def addConservation( self, pfamEntries=None, verbose=0, log=None):
        """
        Adds a conservation score profile from pFam HMMs. See L{Biskit.Hmmer}
	The theoretically most useful one is 'cons_ent' which gives the relative
	entropy of the residue distribution with respect to the background 
	distribution of amino acids (Kullback-Leibler distance) in swissprot.
	See PMID 16916457.
        
        @param pfamEntries: External hmmSearch result, list of
                            (non-overlapping) profile hits.
                            (default: None, do the search) Example::
                              [{'ribonuclease': [[1, 108]]},..]
                              [{profileName : [ [startPos, endPos],
                                                [start2, end2]]}]
                             - startPos, endPos as reported by hmmPfam
                               for PDB sequence generated from this model
        @type  pfamEntries: [{dict}]

        @param verbose: verbosity level (default: 0)
        @type  verbose: 1|0
        @param log: Log file for messages [STDOUT]
        @type  log: Biskit.LogFile
        
        @raise ExeConfigError: if external application is missing
        """
        ## mask out solvent and other troublemakers
	mask = self.m.maskProtein()
	resmask = self.m.atom2resMask( mask )
	
	m = self.m
	if not N.alltrue( mask ):
	    m = self.m.compress( mask )

        h = Hmmer( verbose=verbose, log=log )
        h.checkHmmdbIndex()

        p, hmmHits = h.scoreAbsSum( m, hmmNames=pfamEntries )

        self.m.residues.set( 'cons_abs', p, hmmHits=hmmHits, mask=resmask,
              comment="absolute sum of all 20 hmm scores per position",
              version= T.dateString() + ' ' + self.version() )

        p, hmmHits = h.scoreMaxAll( m, hmmNames=hmmHits )

        self.m.residues.set( 'cons_max', p, hmmHits=hmmHits, mask=resmask,
              comment="max of 20 hmm scores (-average / SD) per position",
              version= T.dateString() + ' ' + self.version() )

        p,  hmmHits = h.scoreEntropy( m, hmmNames=hmmHits )

        self.m.residues.set( 'cons_ent', p, hmmHits=hmmHits, mask=resmask,
              comment="relative entropy (Kullback-Leibler distance) between "\
	      +"observed and background amino acid distribution "\
	      +"(high -> high conservation/discrimination)",
              version= T.dateString() + ' ' + self.version() )


    def addDensity( self, radius=6, minasa=None, profName='density' ):
        """
        Count the number of heavy atoms within the given radius.
        Values are only collected for atoms with |minasa| accessible surface
        area.
        
        @param minasa: relative exposed surface - 0 to 100%
        @type  minasa: float
        @param radius: in Angstrom
        @type  radius: float
        """
        mHeavy = self.m.maskHeavy()

        xyz = N.compress( mHeavy, self.m.getXyz(), 0 )

        if minasa and self.m.profile( 'relAS', 0 ) == 0:
            self.addASA()

        if minasa:
            mSurf = self.m.profile2mask( 'relAS', minasa )
        else:
            mSurf = N.ones( self.m.lenAtoms() )

        ## loop over all surface atoms
        surf_pos = N.nonzero( mSurf )
        contacts = []

        for i in surf_pos:
            dist = N.sum(( xyz - self.m.xyz[i])**2, 1)
            contacts += [ N.sum( N.less(dist, radius**2 )) -1]

        self.m.atoms.set( profName, contacts, mSurf, default=-1,
                               comment='atom density radius %3.1fA' % radius,
                               version= T.dateString() + ' ' + self.version() )

    def addFoldX( self ):
        """
        Adds dict with fold-X energies to PDBModel's info dict.
        See L{Biskit.Fold_X}
        @raise ExeConfigError: if external application is missing
        """
        x = Fold_X( self.m )
        self.m.info['foldX'] = x.run()


    def addSurfaceRacer( self, probe=1.4, vdw_set=1, probe_suffix=0 ):
        """
        Always adds three different profiles as calculated by fastSurf::
           curvature - average curvature (or curvature_1.4 if probe_suffix=1)
           MS - molecular surface area   (or MS_1.4 if probe_suffix=1)
           AS - accessible surface area  (or AS_1.4 if probe_suffix=1)
           
        If the probe radii is 1.4 Angstrom and the Richards vdw radii
        set is used the following two profiles are also added::
           relAS - Relative solvent accessible surface
           relMS - Relative molecular surface
           
        See {Bikit.SurfaceRacer}
        
        @param probe: probe radius
        @type  probe: float
        @param vdw_set: defines what wdv-set to use (1-Richards, 2-Chothia)
        @type  vdw_set: 1|2
        @param probe_suffix: append probe radius to profile names
        @type  probe_suffix: 1|0

        @raise ExeConfigError: if external application is missing
        """
        name_MS   = 'MS' + probe_suffix * ('_%3.1f' % probe)
        name_AS   = 'AS' + probe_suffix * ('_%3.1f' % probe)
        name_curv = 'curvature' + probe_suffix * ('_%3.1f' % probe)

        ## hydrogens are not allowed during FastSurf calculation
        mask = self.m.maskHeavy()
        
        fs = SurfaceRacer( self.m, probe, vdw_set=vdw_set )
        fs_dic = fs.run()

        fs_info= fs_dic['surfaceRacerInfo']

        self.m.atoms.set( name_MS, fs_dic['MS'], mask, 0,
                               comment='Molecular Surface area in A',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        self.m.atoms.set( name_AS, fs_dic['AS'], mask, 0,
                               comment='Accessible Surface area in A',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        self.m.atoms.set( name_curv, fs_dic['curvature'], mask, 0,
                               comment='Average curvature',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        if round(probe, 1) == 1.4 and vdw_set == 1:
            self.m.atoms.set( 'relAS', fs_dic['relAS'], mask, 0,
                                   comment='Relative solvent accessible surf.',
                                   version= T.dateString()+' ' +self.version(),
                                   **fs_info )

            self.m.atoms.set( 'relMS', fs_dic['relMS'], mask, 0,
                                   comment='Relative molecular surf.',
                                   version= T.dateString()+' '+self.version(),
                                   **fs_info )


    def addIntervor( self, cr=[0], cl=None, mode=2, breaks=0, **kw ):
        """
        Triangulate a protein-protein interface with intervor.
        
        @param model: Structure of receptor, ligand and water
        @type  model: Biskit.PDBModel
        @param cr: receptor chains (default: [0] = first chain)
        @type  cr: [ int ]
        @param cl: ligand chains (default: None = all remaining protein chains)
        @type  cl: [ int ]
        @param breaks: consider chain breaks (backbone gaps) (default: 0)
        @type  breaks: bool or 1|0
        @param mode: what to calculate (default 2, = all with shelling order)
        @type  mode: int

        @return: Intervor instance
        @rtype: Biskit.Dock.Intervor
        """
	from Biskit.Dock.Intervor import Intervor

        x = Intervor( self.m, cr=cr, cl=cl, mode=mode, breaks=breaks, **kw)
        x.run()
        
        return x

#############
##  TESTING        
#############
import Biskit.test as BT 

class Test(BT.BiskitTest):
    """Test class """

    TAGS = [ BT.EXE ]
    
    M     = None #: cache processed PDB between tests

    def prepare(self):
        from Biskit import PDBModel
        self.f = T.testRoot() + '/com/1BGS.pdb'

        if not Test.M:
            if self.local: print "Loading PDB...",
            Test.M = PDBModel( self.f )
            Test.M = self.M.compress( Test.M.maskProtein() )

            if self.local: print "Done"
        
        self.d = PDBDope( self.M )

    def test_addFoldX( self ):
        """PDBDope.addFoldX test"""
        self.d.addFoldX()

    def test_addSurfaceRacer(self):
        """PDBDope.addSurfaceRacer/addSurfaceMask test"""
        if self.local: print "Adding SurfaceRacer curvature...",
        self.d.addSurfaceRacer( probe=1.4 )
        if self.local: print 'Done.'

        if self.local: print "Adding surface mask...",
        self.d.addSurfaceMask()
        if self.local: print 'Done.'

    def test_addSecondaryStructure(self):
        """PDBDope.addSecondaryStructure test"""
        self.d.addSecondaryStructure()

    def test_addDensity(self):
        """PDBDope.addDensity test"""
        self.d.addDensity()

##     def test_addIntervor( self ):
##         """PDBDope.addIntervor test // needs waters in structure"""
##         x = self.d.addIntervor()
##         if self.local:
##             x.visualize()

    def test_model(self):
        """PDBDope test final model"""
        from Biskit import PDBModel

        if self.local:
            print '\nData added to info record of model (key -- value):'
            for k in self.d.m.info.keys():
                print '%s -- %s'%(k, self.d.m.info[k])
                
            print '\nAdded atom profiles:'
            print self.M.atoms
            
            print '\nAdded residue  profiles:'
            print self.M.residues

            ## check that nothing has changed
            print '\nChecking that models are unchanged by doping ...'

        m_ref = PDBModel( self.f )
        m_ref = m_ref.compress( m_ref.maskProtein() )
        for k in m_ref.atoms.keys():
            #ref = [ m_ref.atoms[i][k] for i in m_ref.atomRange() ]
            #mod = [ self.M.atoms[i][k] for i in self.M.atomRange() ]
            self.assert_( N.all( m_ref[k] == self.M[k]) )
                
        ## display in Pymol
        if self.local:
            print "Starting PyMol..."
            from Biskit.Pymoler import Pymoler

            pm = Pymoler()
            pm.addPdb( self.M, 'm' )
            pm.colorAtoms( 'm', N.clip(self.M.profile('relAS'), 0.0, 100.0) )
            pm.show()

class LongTest( BT.BiskitTest ):

    TAGS = [ BT.EXE, BT.LONG ]


    def prepare(self):
        from Biskit import PDBModel
        self.f = T.testRoot() + '/com/1BGS.pdb'

        self.M = PDBModel( self.f )
        self.M = self.M.compress( self.M.maskProtein() )
        
        self.d = PDBDope( self.M )

    def test_conservation(self):
        """PDBDope.addConservation (Hmmer) test"""
        self.local = 0

        if self.local: print "Adding conservation data...",
        self.d.addConservation()
        if self.local: print 'Done.'

        ## display in Pymol
        if self.local:
            print "Starting PyMol..."
            from Biskit.Pymoler import Pymoler

            pm = Pymoler()
            pm.addPdb( self.M, 'm' )
            pm.colorAtoms( 'm', N.clip(self.M.profile('cons_ent'), 0.0, 100.0) )
            pm.show()


class OldTest( BT.BiskitTest ):

    TAGS = [ BT.EXE, BT.OLD ]

    def prepare(self):
        from Biskit import PDBModel
        self.f = T.testRoot() + '/com/1BGS.pdb'

        self.M = PDBModel( self.f )
        self.M = self.M.compress( self.M.maskProtein() )

        self.d = PDBDope( self.M )

    def _test_addAsa(self):
        """PDBDope.addAsa (Whatif, obsolete) test"""
        self.d.addASA()


if __name__ == '__main__':

    BT.localTest()  ## pushes test PDBDope instance as 'd' into namespace

##     from Biskit import PDBModel
##     f = T.testRoot() + '/com/1BGS.pdb'

##     M = PDBModel( f )
##     M = M.compress( M.maskProtein() )

    
##     d = PDBDope( M )

##     d.addSecondaryStructure()

##     d.addFoldX()
    
##     d.addDensity()
    
##     d.addSurfaceRacer( probe=1.4 )

##     d.addSurfaceMask()
    
