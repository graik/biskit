##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner

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

## last $Revision: 680 $
## last $Date: 2009-04-25 21:06:41 +0200 (Sat, 25 Apr 2009) $
## last $Author: graik $
"""Convert xplor trajectory into bunch of PDBs"""

import Biskit.tools as t
from Biskit import Xplorer, XplorerError
import os.path as osp
import os

class Xplor2pdb( Xplorer ):
    """
    Simple extraction of PDB files from xplor trajectory files.

    The only complication here is that we have to work around a 80 char
    limit of xplor strings. If input file names are longer than 80 characters,
    we try to shorten them by making them relative to the current working
    directory.

    The PDBs can be converted into a trajectory object using pdb2traj.py.
    """

    def __init__( self, f_psf, f_traj, f_pdb='', skip=0,
                  ascii=True, **kw ):
        """
        @param f_psf: str, psf file for trajectory
        @param f_traj: str, trajectory file name
        @param f_pdb: str, BASE of the output pdbs, we will add ' _xx.pdb'
                      where xx is a running number starting with 0.
        @param skip : int, frame interval [0]
        """
        f_template = t.dataRoot()+'/xplor/crd2pdb.inp'
        Xplorer.__init__( self, template=f_template, **kw )

        ## separate folder from file to avoid xplor's 133 char length limit
        self.psf_folder,   self.f_psf  = osp.split( t.absfile( f_psf ) )
        self.traj_folder,  self.f_traj = osp.split( t.absfile( f_traj ) )
        self.out_folder, self.f_pdb    = osp.split( t.absfile( f_pdb ) )

        self.psf_folder  = self.relativize( self.psf_folder ) + '/'
        self.traj_folder = self.relativize( self.traj_folder ) + '/'
        self.out_folder  = self.relativize( self.out_folder ) + '/'

        self.skip  = skip

        self.ascii = 'true' if ascii else 'false' 


    def relativize( self, path, buffer=15 ):
        """
        Try to shorten path by making it relative to working directory.
        """
        if len( path ) + buffer < 80:
            return path

        cwd = self.cwd or t.absfile( os.curdir )
        path = t.relpath( cwd, path )

        if len( path ) + buffer < 133:
            return path

        raise XplorerError, 'input file name exceeds xplor input length: %s' \
              % path

##     def cleanup( self ):
##         Xplorer.cleanup( self )

##     def finish( self ):
##         pass


if __name__ == '__main__':

    x = Xplor2pdb( t.testRoot()+'/xplor/dynamics.psf',
                   t.testRoot()+'/xplor/dynamics.dcd',
                   f_pdb = t.testRoot() + '/xplor/snapshot', 
                   skip=50,
                   verbose=True,
                   debug=False)
    x.run()
    

    
    
