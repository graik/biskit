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
Read binary Amber trajectory files (NetCDF format) using parmed library.

.. seealso:: `biskit.md.Trajectory`, `biskit.core.TrajParserFactory`,
"""
import logging
import numpy as N

from parmed.amber import NetCDFTraj

import biskit as B
import biskit.tools as T
import biskit.core.trajparser as P
from biskit.core.pdbparseFile import PDBParseFile
from biskit.core.pdbparseModel import PDBParseModel

class TrajParsePDBs(P.TrajParser):
    """
    Generate ``Biskit.md.Trajectory` objects from list of PDB files.
    The parser attempts to correct for minor atom content missmatches. 
    """
    
    ## short free text description of the supported format
    description = 'PDB files with trajectory frames'

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance:

        >>> if TrajParser.supports(['frame_01.pdb', 'frame_02.pdb']):
        >>>     ...

        The test is looking at the first 10 members of the given list and
        verifies that each of them can be parsed either as a PDB file or
        is already a PDBModel instance.
        
        Returns:
            bool: True if the given source is supported by this parser
                  implementation
        """
        ## source is iterable and not empty
        if not (hasattr(source, '__iter__') and len(source) > 0):
            return False
        
        ## test first 10 items in source only
        for f in source[:10]:
            if not (PDBParseFile.supports(f) or PDBParseModel.supports(f)):
                return False
        
        return True
        

    def parse2new( self, source, ref=None, traj=None ):
        """
        Create / Replace Trajectory from the source list of PDBModels or PDBs.
        
        Args:
            source (str): list of file names or PDBModel instances
            ref (str or PDBModel): reference structure instance or file
            traj (Biskit.md.Trajectory): existing instance to be updated

        Returns:
           Biskit.Trajectory: new Trajectory instance
        """
        r = traj
        if traj is None:
            import biskit.md
            r = biskit.md.Trajectory()
        
        r.setRef(B.PDBModel(ref or source[0]))
        n_frames = len(source)
        
        if self.rmwat:
            r.ref = r.ref.compress( N.logical_not(r.ref.maskSolvent()) )

        r.resIndex = r.ref.resMap()
        refNames = r.ref.atomNames()  ## cache for atom checking

        if self.verbose: T.errWrite('reading %i pdbs...' % n_frames )

        r.frames = N.zeros((n_frames,r.ref.lenAtoms(),3)) ## target coordinate array
        r.frameNames = [ '#%i07' % i for i in range(n_frames) ]
        
        atomCast = None
        reportIntervall = 1 if n_frames < 100 else round(n_frames/100)
        
        for i, f in enumerate(source):

            m = B.PDBModel(f)

            ## compare atom order & content of first frame to reference pdb
            if self.analyzeEach or i==0:
                atomCast, castRef = m.compareAtoms( r.ref )

                if castRef != list(range(r.ref.lenAtoms())):
                    ## we can remove/reorder atoms from each frame but not from ref
                    raise P.TrajParserError("Reference PDB doesn't match %s."
                                    %m.fileName)

                if N.all(atomCast == list(range( len( m ))) ):
                    atomCast = None   ## no casting necessary
                else:
                    if self.verbose: T.errWrite(' casting ')

            ## assert that frame fits reference
            if atomCast:
                m = m.take( atomCast )

            ## additional check on each 100st frame
            if i % reportIntervall == 0 and m.atomNames() != refNames:
                raise P.TrajParserError("%s doesn't match reference pdb."%m.fileName )

            r.frames[i] = m.xyz

            if type(f) is str:  ## save original file name 
                r.frameNames[i] = T.stripFilename(f)

            if i % reportIntervall == 0 and self.verbose:
                T.errWrite('#')

        if self.verbose: T.errWrite( 'done\n' )        
        return r
    

##############
# Empty test #
##############

import biskit.test as BT

class Test(BT.BiskitTest):
    
    def test_TrajParsePDBs(self):
        """TrajParsePDBs test"""
        import os
        f = T.testRoot('amber/md_pdbs/')
        allfiles = os.listdir( f )
        pdbs = []
        for fn in allfiles:
            try:
                if (fn[-4:].upper() == '.PDB'):
                    pdbs += [f + fn]
            except:
                pass

        ref = pdbs[0]

        self.assertTrue( TrajParsePDBs.supports(pdbs))
        p = TrajParsePDBs(verbose=self.local, rmwat=True, analyzeEach=False)
        t = p.parse2new(pdbs, ref=ref)
        
        self.assertEqual(t.lenAtoms(),876)
        self.assertEqual(len(t),10)
        self.assertEqual(t.frameNames, [T.stripFilename(f) for f in pdbs])

if __name__ == '__main__':
   
    BT.localTest()
