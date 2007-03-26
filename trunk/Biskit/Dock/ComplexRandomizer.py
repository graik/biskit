## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

## generate random orientations of receptor and ligand
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
"""
Create Complexes with random orientation from a receptor and ligand structure.
"""
    
from Biskit.Dock.Complex import Complex
import Biskit.mathUtils as ma
import Biskit.molUtils as mol
import Biskit.tools as t
import numpy.oldnumeric.random_array as ra
import numpy.oldnumeric as N
from Biskit import Xplorer, PCRModel

import tempfile

class ComplexRandomizer:
    """
    Create Biskit.Dock.Complex(es) with random orientation
    """

    def __init__( self, mrec, mlig, rec_out=None, lig_out=None, debug=0 ):
        """
        @param mrec: receptor model
        @type  mrec: PCRModel
        @param mlig: ligand model
        @type  mlig: PCRModel
        @param rec_out: rec output(default: None)
        @type  rec_out: str
        @param lig_out: lig output (default: None)
        @type  lig_out: str
        @param debug: 1, keep temporary xplor files (default: 0)
        @type  debug: 1|0
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
        
        @param model: model to center
        @type  model: PDBModel
        
        @return: PDBModel (clone of model)
        @rtype: PDBModel
        """
        r = model.clone()
        r.keep( N.nonzero( N.logical_not( r.maskH2O() ) ) )
        center = r.centerOfMass()
        r.setXyz( r.getXyz() - center )

        return r


    def __max_distance( self, model ):
        """
        largest center to any other atom distance
        
        @param model: model with centered coordinates
        @type  model: PDBModel
        
        @return: largest distance
        @rtype: float
        """
        center = model.centerOfMass()
        dist = N.sqrt( N.sum( ( model.getXyz()-center )**2 , 1 ) )

        return max( dist )


    def __random_translation( self ):
        """
        Random translation on a sphere around 0,0,0 with fixed radius
        The radius is the sum of the (max) radius of receptor and ligand
        
        @return: translation array 3 x 1 of float
        @rtype: array
        """
        radius = (self.d_max_rec + self.d_max_lig) / 2.0
        xyz = ra.random( 3 ) - 0.5

        scale = radius*1.0 / N.sqrt( N.sum( xyz**2 ) )

        return scale * xyz


    def __random_matrix( self ):
        """
        Random rotation matrix.
        
        @return: 4 x 4 array of float, random rotation and translation matrix
        @rtype: array
        """
        r = ma.randomRotation()
##         r = N.array([[1,0,0],[0,1,0],[0,0,1]],'f')
        t = self.__random_translation()

        ## create 3 x 4 matrix: 0:3, 0:3 contains rot; 3,0:3 contains trans
        result = N.concatenate( (r, N.transpose( [ t.tolist() ] )), 1)

        ## make it square
        result = N.concatenate( (result, N.array([[0,0,0,1]], N.Float32)), 0 )

        return result


    def random_complex_remote( self ):
        """
        Create a complex where the recrptor and ligand have random
        orientations but are spaced within contact distance.
        
        @return: rec & lig spaced r_rec + r_lig apart in random orientation
        @rtype: Complex
        """
        return Complex( self.rec, self.lig,
                        ligMatrix= self.__random_matrix() )


    def __minimize_complex( self, com ):
        """
        Use Xplor to rigid body minimize the random complex.

        @param com: random complex
        @type  com: Complex
        """
        xp = ComplexMinimizer( com, t.tempDir(), log=self.xp_log )
        xp.run()


    def random_complex( self, inp_mirror=None ):
        """
        @return: randomized and minimized complex
        @rtype: Complex
        """
        self.cm = ComplexMinimizer( self.random_complex_remote(),
                                    debug=self.debug )
        self.cm.run( inp_mirror=inp_mirror )

        com = Complex( self.rec, self.lig )

        rt = com.extractLigandMatrix( self.cm.lig )
        com.setLigMatrix( rt )

        return com


class ComplexMinimizer( Xplorer ):
    """
    Rigid-body minimize receptor and ligand of a Complex
    using soft vdW pot.
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
        """
        Prepare for calculation. Write input files.
        """
        self.com.rec().writePdb( self.rec_in )
        self.com.lig().writePdb( self.lig_in )


    def cleanup( self ):
        """
        Remove temporary files.
        """
        Xplorer.cleanup( self )

        if not self.debug:

            t.tryRemove( self.rec_in )
            t.tryRemove( self.lig_in )

            t.tryRemove( self.rec_out )
            t.tryRemove( self.lig_out )


    def finish( self ):
        """
        When done, write result to disc.
        """
        self.rec = PCRModel( self.com.rec_model.getPsfFile(), self.rec_out )
        self.lig = PCRModel( self.com.lig_model.getPsfFile(), self.lig_out )



#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case
    
    The test generates 3 random complexes. In interactive mode,
    the 3 complexes are displayed as movie in Pymol. They are
    written out as Amber trajectory if debug=True.
    """

    TAGS = [ BT.EXE, BT.LONG ]

    def prepare(self):
        import tempfile
	self.f_pfb = tempfile.mktemp('_test.pdb')
	self.f_crd = tempfile.mktemp('_test.crd')
	
    def cleanUp(self):
	t.tryRemove( self.f_pfb )
	t.tryRemove( self.f_crd )

    def test_ComplexRandomizer(self):
	"""Dock.ComplexRandomizer test"""
        from Biskit import Trajectory

	if self.local:
	    print "\nLoading Rec and Lig files ...",

        rec_pdb = t.testRoot() + '/rec/1A2P.pdb' 
        lig_pdb = t.testRoot() + '/lig/1A19.pdb' 

        rec_psf = t.testRoot() + '/rec/1A2P.psf' 
        lig_psf = t.testRoot() + '/lig/1A19.psf' 

        rec = PCRModel( rec_psf, rec_pdb )
        lig = PCRModel( lig_psf, lig_pdb )

	if self.local:
	    print "Initializing Randomizer..."

        self.cr = ComplexRandomizer( rec, lig, debug=self.DEBUG )

	if self.local:
	    print "Creating 3 random complexes..."

        cs = [ self.cr.random_complex() for i in range(3) ]

        self.traj = Trajectory( [ c.model() for c in cs ] )

        if self.local:
	    self.display( self.traj )
            globals().update( locals() )

        self.assertEqual( len(self.traj), 3 )
    

    def display(self, traj ):
	"""Display random complexes as trajectory in Pymol.
	Only run in local interactive mode.
	"""
	from Biskit import Pymoler

	print "activate debug switch to get random complexes written to disc!"
	if self.DEBUG:
	    print "writing random complex as trajectory to file..."
	    traj.ref.writePdb( self.f_pfb )
	    traj.writeCrd( self.f_crd )
	    print 'Wrote reference pdb file to: %s' % self.f_pfb
	    print 'Wrote crd file to: %s' % self.f_crd

	self.pm = Pymoler( full=0 )

	mname = self.pm.addMovie( [ traj[i] for i in range(len(traj)) ] )
	self.pm.add('hide all')
	self.pm.add('show cartoon')
	self.pm.add('spectrum')
	self.pm.add('mplay')

	self.pm.run()


if __name__ == '__main__':

    BT.localTest()

