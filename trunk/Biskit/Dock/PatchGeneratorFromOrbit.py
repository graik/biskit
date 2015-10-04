##
## Biskit, a toolkit for the manipulation of macromolecular structures
## (C) 2004-2006 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

import numpy.random as ra
import numpy as N

class PatchGeneratorFromOrbit:
    """
    Generate chunks / patches from a given PDBModel. To generate
    surface patches, give a model with only the surface-accessible
    atoms.

    Patches resembling a cutoff-based binding interface are created by
    placing a point at random above the molecule and selecting the
    atoms that are closest to this point. The constant distance of
    this point from the molecule's center (it's orbit) should be set
    to the center-to-center distance between the molecule and its
    binding partner in the reference complex.

    Note that, due to the selection 'from outer space', the patches
    are often discontinuous and not yet a good random reference for
    real binding interfaces. In a_trajPatchDynamics.py, they are only
    used as seed patch which is then re-centered in the same way as
    the interface patch is re-centered.

    See scripts/analysis/a_trajPatchDynamics.py
    """

    def __init__( self, model, orbit, center=None ):
        """
        model  - PDBModel
        orbit  - float, center-to-center distance of binding partner in complex
        exclude- [ 1|0 ], atom mask of atoms to exclude from patches
        """
        self.model   = model
        self.orbit   = orbit ## or max( model.getXyz() - model.center() )
        self.center  = center or model.center()


    def __distances( self, point, xyz=None ):
        """
        point - 3 x 1 array of float; point of origin
        xyz   - 3 x n array of float; coordinates, if None -- take model atoms
        -> distances of all atoms to given point
        """
        if xyz is None:
            xyz = self.model.getXyz()
        return N.sqrt( N.sum( N.power( xyz - point, 2), 1 ) )


    def random_translations( self, n=1, center=None ):
        """
        n Random translations on a sphere around center with fixed radius.
        The radius must be given as orbit to __init__.
        n      - int, number of random coordinates to generate
        center - 3 array of float
        -> array n x 3 of float
        """
        if center is None:
            center = self.center

        xyz = ra.random( (n,3) ) - 0.5

        scale = self.orbit*1.0 / N.sqrt( N.sum( xyz**2, 1 ) )

        r = N.array( [ scale[i]*xyz[i] for i in range(n) ] )

        return r + center


    def patchAround( self, center, nAtoms ):
        """
        patchAround( float_center, int_nAtoms ) -> mask for self.model
        Create single patch of nAtoms atoms that are closest to center.
        """
        dist = self.__distances( center )
        order = N.argsort( dist )

        r = N.zeros( len( self.model ), 'i' )
        N.put( r, order[:nAtoms], 1 )

        return self.centerPatch( r )


    def centerPatch( self, patch_mask ):
        """
        patch_mask - [ 1|0 ], mask of non-centered patch
        -> [ 1|0 ], mask of patch around geometric center of first patch
        """
        c    = self.model.center( patch_mask )
        dist = self.__distances( c )

        n_atoms= len( N.nonzero( patch_mask ) )
        i_dist = N.argsort( dist )[:n_atoms]

        result = N.zeros( len( patch_mask ) )
        N.put( result, i_dist, 1 )

        return result


    def orderCenters( self, points, origin=None ):
        """
        Order random points by increasing distance to first or to origin.
        points  - n x 3 array of float, random center coordinates
        origin  - 3 array of float
        -> [ int ], indices into  points ordered by increasing distance
        """
        origin = origin
        if origin is None:
            origin = points[0]

        dist = self.__distances( origin, points )

        return N.take( points, N.argsort( dist ) )


    def randomPatches( self, size, n=None, exclude=None,
                       max_overlap=0, exclude_all=None ):
        """
        size - int, number of atoms per patch
        n    - int, number of patches (None -> as many as possible, max 100)
        exclude     - [ 1|0 ], don't touch more than |max_overlap| of these
                      atoms (atom mask)
        max_overlap - int
        exclude_all - [ 1|0 ], don't touch ANY of these atoms
        -> [ [ 1|0 ] ], list of atom masks
        """
        if exclude is None:
            exclude = N.zeros( self.model.lenAtoms(), 'i' )

        if exclude_all is None:
            exclude_all = N.zeros( self.model.lenAtoms(), 'i' )

        n = n or 500

        centers = self.random_translations( n=n, center=self.center )

        ## start from excluded patch (if given) working outwards
        origin = centers[0]

        tabu = exclude_all
        if not N.any( tabu ):
            tabu = exclude
        else:
            origin = self.model.center( mask=tabu )

        centers = self.orderCenters( centers, origin )

        r = []

        for i in range(n):

            m = self.patchAround( centers[i], size )

            if N.sum( m * exclude ) <= max_overlap \
               and N.sum( m * exclude_all ) == 0:

                exclude = exclude + m
                r += [ m ]

        return r


def test( model, center2center, nAtoms=10, exclude=None ):

    from Biskit import Pymoler, PDBModel

    g = PatchGeneratorFromOrbit( model, center2center )

    overlap = int( round( nAtoms / 4.0 ) )

    r = g.randomPatches( nAtoms, 500, max_overlap=overlap, exclude=exclude )

    profile = N.sum( N.array(r) )

    pm  = Pymoler()
    pm.addPdb( model, 'all' )

    ms = [ model.take( N.nonzero(mask) ) for mask in r ]

    pm.addMovie( ms )

    return pm



if __name__ == '__main__':

    from Biskit import PDBDope
    from Biskit.tools import *
    from Biskit.Dock import Complex

    print "Loading"
    m_com = load( testRoot() + '/com/ref.complex' ).model()
    rec = m_com.takeChains([0])
    lig = m_com.takeChains([1])

    ## add molecular surface to components
    doper = PDBDope( rec )
    doper.addSurfaceRacer( probe=1.4 )
    surf_rec = rec.profile2mask( 'MS', 0.0001, 101 )

    doper = PDBDope( lig )
    doper.addSurfaceRacer( probe=1.4 )
    surf_lig = lig.profile2mask( 'MS', 0.0001, 101 )

    ## kick out non-surface
    rec = rec.compress( surf_rec )
    lig = lig.compress( surf_lig )

    com = Complex( rec, lig )

    ## get interface patch
    cont = com.atomContacts( cutoff=6.0 )
    rec_if = N.sum( cont, 1 )
    lig_if = N.sum( cont, 0 )

    ## center distance
    c2c = N.sqrt( N.sum( (rec.center() - lig.center())**2, 0 ) )
    print "Center2Center: ", c2c

    ## get patches and put them into Pymoler for display
    print "Patching"
    excl = N.compress( N.ones( len( rec_if ) ), rec_if )
    pm = test( rec, c2c, nAtoms=len(N.nonzero(rec_if)), exclude=rec_if )


    pm.addPdb( rec.compress( rec_if ), 'rec_interface' )
    pm.addPdb( lig.compress( lig_if ), 'lig_interface' )
    pm.addPdb( com.model(), 'complex')

    ## show everything
    ## the patches are as movie in 'model' 
    pm.show()
