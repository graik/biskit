##
## Biskit, a toolkit for the manipulation of macromolecular structures
## (C) 2004-2006 Raik Gruenberg & Johan Leckner; All rights reserved
##


import random
import numpy as N

class PatchGenerator:
    """
    Generate chunks / patches from a given PDBModel. To generate surface
    patches, give a model with only the surface-accessible atoms.
    """

    def __init__( self, model ):
        self.model = model

    def __distances( self, atom ):
        """distances( int_atom ) -> distances of all atoms to given atom"""
        xyz = self.model.xyz
        return N.sqrt( N.sum( N.power( xyz - xyz[atom], 2), 1 ) )
        

    def patchAround( self, atom, size ):
        """
        patchAround( int_atom, int_nAtoms ) -> mask for self.model
        Create single patch of given size around given atom
        """
        dist = self.__distances( atom )
        order = N.argsort( dist )

        r = N.zeros( len( self.model ), 'i' )
        N.put( r, order[:size], 1 )

        return r


    def distantAtoms( self, n, first=None ):
        """
        Select n atoms more or less equaly distributed (the why-not-algorithm)
        n     - int, number of atoms to select
        first - int, first atom (None -> random )
        -> list of int, atom indices
        """
        random.seed()
        xyz = self.model.xyz

        atoms = [ first or random.randint(0, self.model.lenAtoms()-1) ]
        dist = []

        for i in range(1, n):

            dist += [ self.__distances( atoms[-1] ) ]

            mindist = [ min( N.array(dist)[:,i] )
                        for i in range(len(self.model) ) ]

            atoms += [ N.argmax( mindist ) ]

        return atoms
    

    def orderCenters( self, atoms ):
        """
        Order center atoms by increasing distance to first.
        atoms  - [ int ], atom indices selected by distantAtoms
        -> [ int ], atom indices re-ordered
        """
        dist = self.__distances( atoms[0] ).tolist()
        
        pairs = [(dist, a) for a in atoms ]
        pairs.sort()
        return [ x[1] for x in pairs ]
        

    def randomPatches( self, size, n=None, first_atom=None, exclude=None,
                       max_overlap=0, exclude_all=None ):
        """
        size - int, number of atoms per patch
        n    - int, number of patches (None -> as many as possible, max 100)
        first_atom  - int, center first patch around this atom (None ->random)
        exclude     - [ 1|0 ], don't touch more than |max_overlap| of these
                      atoms (atom mask)
        max_overlap - int
        exclude_all - [ 1|0 ], don't touch ANY of these atoms
        -> [ [ 1|0 ] ], list of atom masks
        """
        exclude     = exclude     or N.zeros( self.model.lenAtoms(), 'i' )
        exclude_all = exclude_all or N.zeros( self.model.lenAtoms(), 'i' )

        n = n or 50

        centers = self.orderCenters( self.distantAtoms( n, first_atom ) )
        
        r = []

        for i in range(n):

            m = self.patchAround( centers[i], size )

            if N.sum( m * exclude ) <= max_overlap \
                   and N.sum( m * exclude_all ) == 0:

                exclude += m
                r += [ m ]

        return r


def test( model ):

    from Biskit import Pymoler
    
    g = PatchGenerator( model )
    r = g.randomPatches( 85, 50, max_overlap=25 )

    profile = N.sum( N.array(r) )

    pm  = Pymoler()
    pm.addPdb( model, 'o' )

    ms = [ model.take( N.flatnonzero(mask) ) for mask in r ]

    pm.addMovie( ms )

    return pm


if __name__ == '__main__':

    from Biskit import PDBDope
    from Biskit.tools import *
    from Biskit.Dock import Complex

    print "Loading"
##     m_com = load( testRoot() + '/com_wet/dry_com.model' )
    m_com = load( testRoot() + '/com/ref.complex').model()
    m_com._resIndex = None   ## HACK, something wrong with legacy _resIndex
    
    doper = PDBDope( m_com )
    doper.addSurfaceRacer()
    
    mask = m_com.profile2mask( 'relAS', 5, 101 )
    m = m_com.compress( mask )

    ## get patches and put them into Pymoler for display
    print "Patching"
    pm = test( m )

    ## show real interface patch
    com = Complex( m_com.takeChains([0]), m_com.takeChains([1]))
    cont = com.atomContacts()
    rec_i = N.flatnonzero( N.sum( cont, 1 ) )
    lig_i = N.flatnonzero( N.sum( cont, 0 ) )

    pm.addPdb( m_com.take( rec_i ), 'rec_interface' )
    pm.addPdb( m_com.takeChains([1]).take( lig_i ), 'lig_interface' )

    ## show everything
    pm.show()
