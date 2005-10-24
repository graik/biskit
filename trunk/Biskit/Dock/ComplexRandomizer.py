## generate random orientations of receptor and ligand
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
##
## $Revision$
## last $Author$
## last $Date$

from Complex import Complex
import Biskit.mathUtils as ma
import Biskit.molUtils as mol
import Biskit.tools as t
import RandomArray as ra
import Numeric as N
from Biskit import Xplorer, PCRModel

import tempfile

class ComplexRandomizer:
    """
    Create Biskit.Dock.Complex(es) with random orientation
    """

    def __init__( self, mrec, mlig, rec_out=None, lig_out=None, debug=0 ):
        """
        mrec - PCRModel
        mlig - PCRModel
        """
        ## rec and lig with centered coordinates
        self.rec = self.__center_model( mrec )
        self.lig = self.__center_model( mlig )

        ## this way they will be unique in ComplexList
        if rec_out:
            self.rec.saveAs( rec_out )
        if lig_out:
            self.lig.saveAs( lig_out )

        ## get max. center-to-atom distances
        self.d_max_rec = self.__max_distance( self.rec )
        self.d_max_lig = self.__max_distance( self.lig )

        self.xp_log = tempfile.mktemp('rb_min_xplor_log')

        ## keep temporary xplor files
        self.debug = debug

    def __center_model( self, model ):
        """
        translate PDBModel so that it's center is in 0,0,0
        model - PDBModel
        -> PDBModel (clone of model)
        """
        r = model.clone()
        r.keep( N.nonzero( N.logical_not( r.maskH2O() ) ) )
        center, mass = r.centerOfMass(), r.mass()
        r.setXyz( r.getXyz() - center )

        return r


    def __max_distance( self, model ):
        """
        largest center <-> atom distance
        model - PDBModel, with centered coordinates
        -> float
        """
        center, mass = model.centerOfMass(), model.mass()
        dist = N.sqrt( N.sum( ( model.getXyz()-center )**2 , 1 ) )
        
        return max( dist )


    def __random_translation( self ):
        """
        Random translation on a sphere around 0,0,0 with fixed radius
        The radius is the sum of the (max) radius of receptor and ligand
        -> array 3 x 1 of float
        """
        radius = (self.d_max_rec + self.d_max_lig) / 2.0
        xyz = ra.random( 3 ) - 0.5

        scale = radius*1.0 / N.sqrt( N.sum( xyz**2 ) )

        return scale * xyz


    def __random_matrix( self ):
        """
        -> 4 x 4 array of float, random rotation and translation matrix
        """
        r = ma.randomRotation()
##         r = N.array([[1,0,0],[0,1,0],[0,0,1]],'f')
        t = self.__random_translation()

        ## create 3 x 4 matrix: 0:3, 0:3 contains rot; 3,0:3 contains trans
        result = N.concatenate( (r, N.transpose( [ t.tolist() ] )), 1)

        ## make it square
        result = N.concatenate( (result, N.array([[0,0,0,1]],'f')), 0 )

        return result


    def random_complex_remote( self ):
        """
        -> Complex, rec & lig spaced r_rec + r_lig apart in random orientation
        """
        return Complex( self.rec, self.lig,
                        ligMatrix= self.__random_matrix() )


    def __minimize_complex( self, com ):

        xp = ComplexMinimizer( com, t.tempDir(), log=self.xp_log )
        xp.run()
        

    def random_complex( self, inp_mirror=None ):

        self.cm = ComplexMinimizer( self.random_complex_remote(),
                                    debug=self.debug )
        self.cm.run( inp_mirror=inp_mirror)

        com = Complex( self.rec, self.lig )

        rt = com.extractLigandMatrix( self.cm.lig )
        com.setLigMatrix( rt )

        return com


class ComplexMinimizer( Xplorer ):
    """
    Rigid-body minimize receptor and ligand of a Complex using soft vdW pot.
    """

    def __init__( self, com, debug=0, **params ):

        self.com = com
        
        self.rec_psf = com.rec().getPsfFile()
        self.lig_psf = com.lig().getPsfFile()

        recCode = com.rec().getPdbCode()
        ligCode = com.lig().getPdbCode()

        self.rec_in = tempfile.mktemp( recCode + ".pdb" )
        self.lig_in = tempfile.mktemp( ligCode + ".pdb" )

        self.lig_out = tempfile.mktemp( "lig_out.pdb" )
        self.rec_out = tempfile.mktemp( "rec_out.pdb" )
        
        self.inp_template = t.projectRoot() +\
                            '/external/xplor/rb_minimize_complex.inp'

        self.param19 = t.projectRoot() + \
                        '/external/xplor/toppar/param19.pro'

        self.result = None

        Xplorer.__init__( self, self.inp_template, debug=debug, **params )

    def prepare( self ):
        
        self.com.rec().writePdb( self.rec_in )
        self.com.lig().writePdb( self.lig_in )

        
    def cleanup( self ):

        Xplorer.cleanup( self )

        if not self.debug:

            t.tryRemove( self.rec_in )
            t.tryRemove( self.lig_in )

            t.tryRemove( self.rec_out )
            t.tryRemove( self.lig_out )


    def finish( self ):

        self.rec = PCRModel( self.com.rec_model.getPsfFile(), self.rec_out )
        self.lig = PCRModel( self.com.lig_model.getPsfFile(), self.lig_out )



#####TEST######

if __name__ == '__main__':

    from Biskit.tools import *
    from Biskit import Trajectory

    rec = Load( '~/interfaces/c15/rec_wet/1AVV.model' )
    lig = Load( '~/interfaces/c15/lig_wet/1SHF.model' )

    inp_copy = absfile( 'test_randomizer_xplor.in')

    cr = ComplexRandomizer( rec, lig, debug=1 )

    cs = [ cr.random_complex( inp_mirror=inp_copy ) for i in range(5) ]

    traj = Trajectory( [ c.model() for c in cs ] )

    traj.ref.writePdb( '~/test.pdb' )
    traj.writeCrd( '~/test.crd' )

    
