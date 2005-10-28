##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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

import Numeric as N
import Biskit.tools as T
import os.path

from Biskit.WhatIf import WhatIf 
from Biskit.Hmmer import Hmmer

from Biskit.Fold_X import Fold_X
from Biskit.SurfaceRacer import SurfaceRacer


class PDBDope:
    """
    Decorate a PDBModel with calculated properties (profiles)
    """

    def __init__( self, model ):
        self.m = model

    def version( self ):
        return 'PDBDope $Revision$'

    def model( self ):
        return self.m

    def addASA( self ):
        """
        Add profiles of Accessible Surface Area: 'relASA', 'ASA_total',
        'ASA_sc', 'ASA_bb'
        !! ProfileError, if whatif-returned atom/residue lists don't match.
           Usually that means, whatif didn't recognize some residue name
        """
        w = WhatIf( self.m )
        
        atomRelAcc, resASA, resMask = w.run()

##         normalAtoms = N.logical_not( N.logical_or(self.m.maskHetatm(),
##                                                   self.m.maskSolvent() ) )

        normalAtoms = self.m.maskProtein( standard=1 )

        normalRes = self.m.atom2resMask( normalAtoms )

        self.m.setAtomProfile( 'relASA', atomRelAcc, ## normalAtoms, 0,
                               comment='relative accessible surface area in %',
                               version= T.dateString() + ' ' + self.version() )

        self.m.setResProfile( 'ASA_total', resASA[:,0], normalRes, 0,
                              comment='accessible surface area in A^2',
                              version= T.dateString() + ' ' + self.version() )

        self.m.setResProfile( 'ASA_sc', resASA[:,1], normalRes, 0,
                           comment='side chain accessible surface area in A^2',
                           version= T.dateString() + ' ' + self.version() )

        self.m.setResProfile( 'ASA_bb', resASA[:,2], normalRes, 0,
                           comment='back bone accessible surface area in A^2',
                           version= T.dateString() + ' ' + self.version() )
        

    def addSurfaceMask( self ):
        r = self.m.profile2mask( 'relASA', cutoff_min=40 )
        self.m.setResProfile( 'surfMask',  self.m.atom2resMask(r),
                              comment='residues with any atom > 40% exposed',
                              version= T.dateString() + ' ' + self.version() )
        

    def addConservation( self, pfamEntries=None ):
        """
        pfamEntries - external hmmSearch result, list of (non-overlapping)
                      profile hits
                      e.g. [{'ribonuclease': [[1, 108]]},..]
                      [{profileName : [ [startPos, endPos],[start2, end2]]}]
                    - startPos, endPos as reported by hmmPfam for PDB sequence
                      generated from this model
                    - default: None, do the search
        """
        ## mask with normal AA also used for HMM search
        mask = self.m.maskCA()
        mask = self.m.atom2resMask( mask )

        h = Hmmer()
        h.checkHmmdbIndex()

        p, hmmHits = h.scoreAbsSum( self.m, hmmNames=pfamEntries )

        self.m.setResProfile( 'cons_abs', p, mask, 0, hmmHits=hmmHits,
              comment="absolute sum of all 20 hmm scores per position",
              version= T.dateString() + ' ' + self.version() )

        p, hmmHits = h.scoreMaxAll( self.m, hmmNames=hmmHits )

        self.m.setResProfile( 'cons_max', p, mask, 0, hmmHits=hmmHits,
              comment="max of 20 hmm scores (-average / SD) per position",
              version= T.dateString() + ' ' + self.version() )

        p,  hmmHits = h.scoreEntropy( self.m, hmmNames=hmmHits )

        self.m.setResProfile( 'cons_ent', p, mask, 0, hmmHits=hmmHits,
              comment="entropy of emmission probabilities per position "+
                              "(high -> high conservation/discrimination)",
              version= T.dateString() + ' ' + self.version() )


    def addDensity( self, radius=6, minasa=None, profName='density' ):
        """
        Count the number of heavy atoms within the given radius.
        Values are only collected for atoms with |minasa| accessible surface
        area.
        minasa - float, 0 to 100%
        radius - float, in A
        """
        mHeavy = self.m.maskHeavy()

        xyz = N.compress( mHeavy, self.m.getXyz(), 0 )
        
        if minasa and self.m.profile( 'relASA', 0 ) == 0:
            self.addASA()

        if minasa:
            mSurf = self.m.profile2mask( 'relASA', minasa )
        else:
            mSurf = N.ones( self.m.lenAtoms() )

        ## loop over all surface atoms
        surf_pos = N.nonzero( mSurf )
        contacts = []

        for i in surf_pos:
            dist = N.sum(( xyz - self.m.xyz[i])**2, 1)
            contacts += [ N.sum( N.less(dist, radius**2 )) -1]

        self.m.setAtomProfile( profName, contacts, mSurf, default=-1,
                               comment='atom density radius %3.1fA' % radius,
                               version= T.dateString() + ' ' + self.version() )

    def addFoldX( self ):
        """
        adds dict with fold-X energies to PDBModel's info dict.
        """
        x = Fold_X( self.m )
        self.m.info['foldX'] = x.run()


    def addSurfaceRacer( self, probe=1.4, vdw_set=1, probe_suffix=0 ):
        """
        probe        - float, probe radius
        probe_suffix - 1|0, append probe radius to profile names
        Adds three different profiles as calculated by fastSurf
           curvature - average curvature (or curvature_1.4 if probe_suffix=1)
           MS - molecular surface area   (or MS_1.4 if probe_suffix=1)
           AS - accessible surface area  (or AS_1.4 if probe_suffix=1)
        """
        name_MS   = 'MS' + probe_suffix * ('_%3.1f' % probe)
        name_AS   = 'AS' + probe_suffix * ('_%3.1f' % probe)
        name_curv = 'curvature' + probe_suffix * ('_%3.1f' % probe)

        ## hydrogens are not allowed during FastSurf calculation
        m = self.m.clone()
        mask = m.maskHeavy()
        m = m.compress( mask )

        
        fs = SurfaceRacer( m, probe, vdw_set=vdw_set )
        fs_dic = fs.run()

        fs_info= fs_dic['surfaceRacerInfo']

        
        self.m.setAtomProfile( name_MS, fs_dic['MS'], mask, 0,
                               comment='Molecular Surface area in A',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        self.m.setAtomProfile( name_AS, fs_dic['AS'], mask, 0,
                               comment='Accessible Surface area in A',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        self.m.setAtomProfile( name_curv, fs_dic['curvature'], mask, 0,
                               comment='Average curvature',
                               version= T.dateString() + ' ' + self.version(),
                               **fs_info )

        if probe == 1.4:
            self.m.setAtomProfile( 'relAS', fs_dic['relAS'], mask, 0,
                                   comment='Relative solvent accessible surf.',
                                   version= T.dateString()+' ' +self.version(),
                                   **fs_info )
        
            self.m.setAtomProfile( 'relMS', fs_dic['relMS'], mask, 0,
                                   comment='Relative molecular surf.',
                                   version= T.dateString()+' '+self.version(),
                                   **fs_info )
        
if __name__ == '__main__':

    from Biskit import PDBModel
    import glob

    print "Loading PDB..."
    f = glob.glob( T.testRoot()+'/lig_pcr_00/pcr_00/*_1_*pdb' )[1]
    f = T.testRoot()+"/com/1BGS.pdb"
    m = PDBModel(f)
    m = m.compress( m.maskProtein() )

    print "Initiating PDBDope...",
    d = PDBDope( m )
    print 'Done.\n'

    print "Adding FoldX energy...",
    d.addFoldX()
    print 'Done.'
    
    print "Adding WhatIf ASA...",
    d.addASA()
    print 'Done.'

    print "Adding surface mask...",
    d.addSurfaceMask()
    print 'Done.'
    
    print "Adding conservation data...",
    d.addConservation()
    print 'Done.'
    
    print "Adding surface density...",
    d.addDensity()
    print 'Done.'
    
    print "Adding SurfaceRacer curvature...",
    d.addSurfaceRacer( probe=1.4 )
    print 'Done.'
    
    print d.m.info

    ## check that nothing has changed
    print '\nChecking that models are unchanged by doping ...'
    m_ref = PDBModel(f)
    m_ref = m_ref.compress( m_ref.maskProtein() )
    for k in m_ref.atoms[0].keys():
        ref = [ m_ref.atoms[i][k] for i in range( m_ref.lenAtoms() ) ]
        mod = [ m.atoms[i][k] for i in range( m.lenAtoms() ) ]
        if not ref == mod:
            print 'CHANGED!! ', k
        if ref == mod:
            print 'Unchanged ', k


    ## display in Pymol
    print "Starting PyMol..."
    from Biskit.Pymoler import Pymoler

    pm = Pymoler()
    pm.addPdb( m, 'm' )
    pm.colorAtoms( 'm', N.clip(m.profile('relAS'), 0.0, 100.0) )
    pm.show()
