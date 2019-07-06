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
Read classic (text) Amber trajectory files using parmed library.

.. seealso:: `biskit.md.Trajectory`, `biskit.core.TrajParserFactory`,
"""
import logging
import numpy as N

from parmed.amber import AmberMdcrd

import biskit as B
import biskit.tools as T
import biskit.core.trajparser as P

class TrajParseAmberCrd(P.TrajParser):
    """
    Generate ``Biskit.md.Trajectory` from Amber ASCII coordinates file.
    """
    
    ## short free text description of the supported format
    description = 'Amber ASCII-format trajectory file'
    
    def __init__(self, verbose=False, rmwat=False, analyzeEach=False, 
                 hasbox=True):
        super().__init__(verbose=verbose, rmwat=rmwat, analyzeEach=analyzeEach)
        self.hasbox = hasbox

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance:

        >>> if TrajParser.supports('myfile.crd'):
        >>>     ...

        Simplistic implementation that only looks for file ending "ncdf" or 
        "ncrd".

        Returns:
            bool: True if the given source is supported by this parser
                  implementation
        """
        return (type(source) is str or isinstance(source, B.LocalPath)) and \
            (source[-4:].upper() == '.CRD')
        

    def parse2new( self, source, ref, traj=None ):
        """
        Replace content of a new Trajectory from the source.
        
        Args:
            source (str): file name or other input object
            ref (str or PDBModel): reference structure instance or file
            traj (Biskit.md.Trajectory): existing instance to be updated

        Returns:
           Biskit.Trajectory: new Trajectory instance
        """
        r = traj
        if traj is None:
            import biskit.md
            r = biskit.md.Trajectory()
        
        ref = B.PDBModel(ref)
        
        src = AmberMdcrd(source, natom=ref.lenAtoms(),hasbox=self.hasbox)
        
        r.frames = src.coordinates
        
        r.setRef(ref)
        r.resIndex = r.ref.resMap()
        assert N.shape(r.frames) == (src.frame, src.natom, 3)
        assert r.lenAtoms() == src.natom
        
        return r
    

##############
# Empty test #
##############

import biskit.test as BT

class Test(BT.BiskitTest):
    
    def test_TrajParserAmberCrd(self):
        """TrajParserNetCDF test"""
        fcrd = T.testRoot('lig_pcr_00/raw/traj.crd')
        fpdb = T.testRoot('lig_pcr_00/raw/traj_ref.pdb')
        self.assertTrue( TrajParseAmberCrd.supports(fcrd))

        p = TrajParseAmberCrd(hasbox=False)
        t = p.parse2new( fcrd, fpdb)

        self.assertEqual(t.lenAtoms(),876)
        self.assertEqual(len(t),100)

if __name__ == '__main__':
   
    BT.localTest()
