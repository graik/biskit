## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
Calculate and add various properties to PDBModel
"""

import biskit.core.oldnumeric as N0
import numpy as N

import biskit
import biskit.tools as T

## from biskit.exe.hmmer import Hmmer
from biskit.exe.dssp import Dssp
from biskit.exe.surfaceRacer import SurfaceRacer
from biskit.exe.delphi import Delphi, DelphiError


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
        return 'PDBDope ' + biskit.__version__

    def model( self ):
        """
        @return: model
        @rtype: PDBModel
        """
        return self.m


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
        rmodel = dssp.run()

        self.m.residues.set( 'secondary',  rmodel['dssp'], 
                             comment='secondary structure from DSSP',
                             version= T.dateString() + ' ' + self.version(),
                             default='.')
        self.m.residues.set('dssp_acc', rmodel['dssp_acc'],
                            comment='accessible surface area from DSSP')
        self.m.residues.set('dssp_phi', rmodel['dssp_phi'],
                            comment='PHI angle from DSSP',
                            default=360.0)
        self.m.residues.set('dssp_psi', rmodel['dssp_psi'],
                            comment='PSI angle from DSSP',
                            default=360.0)
                            


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
        if not N0.alltrue( mask ):
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

        xyz = N0.compress( mHeavy, self.m.getXyz(), 0 )

        if minasa and self.m.profile( 'relAS', 0 ) == 0:
            self.addASA()

        if minasa:
            mSurf = self.m.profile2mask( 'relAS', minasa )
        else:
            mSurf = N0.ones( self.m.lenAtoms() )

        ## loop over all surface atoms
        surf_pos = N0.nonzero( mSurf )
        contacts = []

        for i in surf_pos:
            dist = N0.sum(( xyz - self.m.xyz[i])**2, 1)
            contacts += [ N0.sum( N0.less(dist, radius**2 )) -1]

        self.m.atoms.set( profName, contacts, mSurf, default=-1,
                          comment='atom density radius %3.1fA' % radius,
                          version= T.dateString() + ' ' + self.version() )


    def addSurfaceRacer( self, probe=1.4, vdw_set=1, probe_suffix=0, mask=None ):
        """
        Always adds three different profiles as calculated by fastSurf::
           curvature - average curvature (or curvature_1.4 if probe_suffix=1)
           MS - molecular surface area   (or MS_1.4 if probe_suffix=1)
           AS - accessible surface area  (or AS_1.4 if probe_suffix=1)

        If the probe radii is 1.4 Angstrom and the Richards vdw radii
        set is used the following two profiles are also added::
           relAS - Relative solvent accessible surface
           relMS - Relative molecular surface

        See {Biskit.SurfaceRacer}

        @param probe: probe radius
        @type  probe: float
        @param vdw_set: defines what wdv-set to use (1-Richards, 2-Chothia)
        @type  vdw_set: 1|2
        @param probe_suffix: append probe radius to profile names
        @type  probe_suffix: 1|0
        @param mask: optional atom mask to apply before calling surface racer
                     (default: heavy atoms AND NOT solvent)
        @type mask: [ bool ]

        @raise ExeConfigError: if external application is missing
        """
        name_MS   = 'MS' + probe_suffix * ('_%3.1f' % probe)
        name_AS   = 'AS' + probe_suffix * ('_%3.1f' % probe)
        name_curv = 'curvature' + probe_suffix * ('_%3.1f' % probe)

        ## hydrogens + waters are not allowed during FastSurf calculation
        mask = mask if mask is not None else \
            self.m.maskHeavy() * N0.logical_not( self.m.maskSolvent() )
        
        fs = SurfaceRacer( self.m, probe, vdw_set=vdw_set, mask=mask )
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

        if round(probe, 1) == 1.4 and vdw_set == 1 and 'relAS' in fs_dic:
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
        from biskit.dock.Intervor import Intervor

        x = Intervor( self.m, cr=cr, cl=cl, mode=mode, breaks=breaks, **kw)
        x.run()

        return x
    
    def addDelphi( self, **kw ):
        """
        Calculate electrostatic potentials and potential maps with Delphi.
        See L{Biskit.Delphi} for details. The same options apply. By default,
        the DelPhi wrapper uses the program reduce to assign hydrogens and then
        maps Amber partial atomic charges into the structure. The result of the
        calculation is stored in the info dict of the model. Example:
        
        >>> dope = PDBDope( model )
        >>> dope.addDelphi()
        >>> print model.info['delphi']
        
            {'scharge' :  1.4266   # surface charge
            'egrid' :  9105.51    # total grid energy
            'ecoul' :  -9849.664  # couloumb energy
            'erxn'  :  -664.7469   # corrected reaction field energy
            'erxnt' :  -21048.13  # total reaction field energy
            'eself' :  -20383.39 }  # self reaction field energy
            
        The same dictionary is also returned by this method.
        
        @param f_map   : output file name for potential map [None= discard]
        @type  f_map   : str
        @param addcharge: build atomic partial charges with AtomCharger
                          [default: True]
        @type  addcharge: bool
        
        @param protonate: (re-)build hydrogen atoms with reduce program (True)
                          see L{Biskit.Reduce}
        @type  protonate: bool
        @param autocap: add capping NME and ACE residues to any (auto-detected)
                        false N- or C-terminal and chain breaks (default: False)
                        see L{Biskit.Reduce} and L{Biskit.PDBCleaner}
        @type  autocap: bool

        @param indi: interior dilectric (4.0)
        @param exdi: exterior dielectric (80.0)
        @param salt: salt conc. in M (0.15)
        @param ionrad: ion radius (2)
        @param prbrad: probe radius (1.4) 
        @param bndcon: boundary condition (4, delphi default is 2)
        @param scale:  grid spacing (2.3)
        @param perfil: grid fill factor in % (for automatic grid, 60) 

        @raise ExeConfigError: if external application (delphi, reduce) is 
                               missing
        
        @return: dict with delphi results
        @rtype: {str: float}
        """
        d = Delphi( self.m, **kw )
        r = d.run()
        
        self.m.info['delphi'] = r
        return r

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    TAGS = [ BT.EXE ]

    M     = None #: cache processed PDB between tests

    def prepare(self):
        from biskit import PDBModel
        self.f = T.testRoot() + '/com/1BGS.pdb'

        if not Test.M:
            if self.local: print("Loading PDB...", end=' ')
            Test.M = PDBModel( self.f )
            Test.M = self.M.compress( Test.M.maskProtein() )

            if self.local: print("Done")

        self.d = PDBDope( self.M )

    def test_addSurfaceRacer(self):
        """PDBDope.addSurfaceRacer/addSurfaceMask test"""
        if self.local: print("Adding SurfaceRacer curvature...", end=' ')
        self.d.addSurfaceRacer( probe=1.4 )
        if self.local: print('Done.')

        if self.local: print("Adding surface mask...", end=' ')
        self.d.addSurfaceMask()
        if self.local: print('Done.')
        
    def test_surfaceRacerBug(self):
        """PDBDope.addSurfaceRacer mask handling bugfix"""
        import os
        from biskit import PDBModel
        m = PDBModel(os.path.join(T.testRoot(),'lig/1A19.pdb'))
        d = PDBDope(m)
        d.addSurfaceRacer()

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
        from biskit import PDBModel

        if self.local:
            print('\nData added to info record of model (key -- value):')
            for k in list(self.d.m.info.keys()):
                print('%s -- %s'%(k, self.d.m.info[k]))

            print('\nAdded atom profiles:')
            print(self.M.atoms)

            print('\nAdded residue  profiles:')
            print(self.M.residues)

            ## check that nothing has changed
            print('\nChecking that models are unchanged by doping ...')

        m_ref = PDBModel( self.f )
        m_ref = m_ref.compress( m_ref.maskProtein() )
        for k in list(m_ref.atoms.keys()):
            #ref = [ m_ref.atoms[i][k] for i in m_ref.atomRange() ]
            #mod = [ self.M.atoms[i][k] for i in self.M.atomRange() ]
            self.assertTrue( N.all( m_ref[k] == self.M[k]) )

        ## display in Pymol
        if self.local:
            print("Starting PyMol...")
            from biskit.exe.pymoler import Pymoler

            pm = Pymoler()
            pm.addPdb( self.M, 'm' )
            pm.colorAtoms( 'm', N0.clip(self.M.profile('relAS'), 0.0, 100.0) )
            pm.show()

class LongTest( BT.BiskitTest ):

    TAGS = [ BT.EXE, BT.LONG ]


    def prepare(self):
        from biskit import PDBModel
        self.f = T.testRoot() + '/com/1BGS.pdb'

        self.M = PDBModel( self.f )
        self.M = self.M.compress( self.M.maskProtein() )

        self.d = PDBDope( self.M )

    ## test de-activated as we don't have a running Hmmer module any longer
    def __test_conservation(self):
        """PDBDope.addConservation (Hmmer) test"""
        if self.local: print("Adding conservation data...", end=' ')
        self.d.addConservation()
        if self.local: print('Done.')

        ## display in Pymol
        if self.local:
            print("Starting PyMol...")
            from biskit.exe.pymoler import Pymoler

            pm = Pymoler()
            pm.addPdb( self.M, 'm' )
            pm.colorAtoms( 'm', N0.clip(self.M.profile('cons_ent'), 0.0, 100.0) )
            pm.show()

    def test_delphi(self):
        """PDBDope.addDelphi test"""
        if self.local:
            self.log.add( 'Calculating Delphi electrostatic potential' )
            self.log.add( '' )

        self.d.addDelphi( scale=1.2 )
        
        self.assertAlmostEqual( self.M['delphi']['scharge'], 0.95, 1 )
                

    ## EXPECT_1RQ4 = '.EEEEE.SSHHHHHHHHHHHHHHHEEEEEEEE.SS.EEEEE..SS...EEEEEEE.SSTTT.....EEEEEESSS..EEEEEETTTTEEEE.GGGTT...TT.EEEE.SS.SSHHHHHHHHTS.STT.EE.HHHHHHHHHHHHT..SS...HHHHHHHHHHHHHTHHHHH.HHHHHHHHGGGT.TT...EE..HHHHHHHTTHHHHHHHGGG.SS.SEEEETTEEESSHHHHHTT..EE.......SS.SSSS..EEEETTEEEEHHHHHHH.....EEEEE.SSHHHHHHHHHHHHHHHEEEEEEEEETTEEEEEE...S...EEEEEEE.SSSSS.....EEEEEETTT..EEEEEETTTTEEEE.GGGTT...TT.EEEE.SS..SHHHHHHHS...TTT.EE.TTHHHHHHHHHHT..SS...HHHHHHHHHHHHHHHHHHH.HHHHHHHGGGSS.TT...EE..HHHHHHHHTHHHHHHHGGG.SS.SEEEETTEEE.SHHHHHTT..EE........S.STTSS.EEESSSEEEEHHHHHHH.......EEEEE.EEEEEE.TTS.EEEEESS.EEEE..HHHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEE...EEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT..EEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHHT.EEEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEE.EEEEEE.SSS.EEEEETTEEEEE..HHHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..HHHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEE.EEEEEE.TTS.EEEEETTB..EE..TTHHHHHHHHHHHT..EEEE.SS.STT.B..EEEE.....EEE.EEEEEE.TTS.EEEEETTEEEEE..THHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEE.EEEEEE.SSS.EEEEETTEEEEE..TTHHHHHHHHHTTT..EEEE.S..STT.B..EEEE.'
    EXPECT_1RQ4 = '.EEEEE.SSHHHHHHHHHHHHHHHEEEEEEEE.SS.EEEEE..SS...EEEEEEE.SSTTT.....EEEEEESSS..EEEEEETTTTEEEE.GGGTT...TT.EEEE.SS.SSHHHHHHHHTS.STT.EE.HHHHHHHHHHHHT..SS...HHHHHHHHHHHIIIIIHHH.HHHHHHHHGGGT.TT...EE..HHHHHHHTTHHHHHHHGGG.SS.SEEEETTEEESSHHHHHTT..EE.......SS.SSSS..EEEETTEEEEHHHHHHH.....EEEEE.SSHHHHHHHHHHHHHHHEEEEEEEEETTEEEEEE...S...EEEEEEE.SSSSS.....EEEEEETTT..EEEEEETTTTEEEE.GGGTT...TT.EEEE.SS..SHHHHHHHS...TTT.EE.TTHHHHHHHHHHT..SS...HHHHHHHHHHHIIIIIHHH.HHHHHHHGGGSS.TT...EE..HHHHHHHHTHHHHHHHGGG.SS.SEEEETTEEE.SHHHHHTT..EE........S.STTSS.EEESSSEEEEHHHHHHH.......EEEEE.EEEEEE.TTS.EEEEESS.EEEE..HHHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEE...EEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT..EEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHHT.EEEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEE.EEEEEE.SSS.EEEEETTEEEEE..HHHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..HHHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEEEEEEEEE.TTS.EEEEETTEEEEE..TTHHHHHHHHHHTT.EEEEE.S..STT.B..EEEE...EEEEE.EEEEEE.TTS.EEEEETTB..EE..TTHHHHHHHHHHHT..EEEE.SS.STT.B..EEEE.....EEE.EEEEEE.TTS.EEEEETTEEEEE..THHHHHHHHHHHHT..EEEE.S..STT.B..EEEE...EEEEE.EEEEEE.SSS.EEEEETTEEEEE..TTHHHHHHHHHTTT..EEEE.S..STT.B..EEEE.'
    def test_addSecondaryNonFiltered(self):
        from biskit import PDBModel
        m = PDBModel('1R4Q')
        d = PDBDope(m)
        d.addSecondaryStructure()
        
        r = m.compress(m.maskProtein())['secondary']
        self.assertEqual(''.join(r), self.EXPECT_1RQ4 )
    

if __name__ == '__main__':

    BT.localTest()  ## pushes test PDBDope instance as 'd' into namespace

