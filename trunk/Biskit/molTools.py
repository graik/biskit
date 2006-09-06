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
## last $Author$
## last $Date$
## $Revision$
"""
Various structure related calculations.
"""

import Numeric as N
from Biskit import molUtils as molU
import Biskit.tools as T

def hbonds( model ):
    """
    Collect a list with all potential hydrogen bonds in model.

    @param model: PDBModel for which 
    @type  model: PDBModel

    @return: a list of potential hydrogen bonds containing a lists
             with donor index, acceptor index, distance and angle.
    @rtype: [ int, int, float, float ]
    """
    hbond_lst = []
    donors = molU.hbonds['donors']
    accept = molU.hbonds['acceptors']

    ## indices if potential donors
    d_ind = []
    for res , aList in donors.items():  
        for a in aList:
            if a in molU.hydrogenSynonyms.keys():
                aList.append( molU.hydrogenSynonyms[a] )
                
        d_ind += model.filterIndex( residue_name=res, name=aList )
        
    ## indices if potential acceptors
    a_ind = []
    for res , aList in accept.items():               
        a_ind += model.filterIndex( residue_name=res, name=aList )
        
    ## calculate pairwise distances and angles
    for d in d_ind:
        d_xyz  = model.xyz[d]
        d_nr   = model.atoms[d]['residue_number']
        d_cid  = model.atoms[d]['chain_id']
        d_segi = model.atoms[d]['segment_id']

        for a in a_ind:
            a_xyz  = model.xyz[a]
            a_nr   = model.atoms[a]['residue_number']
            a_cid  = model.atoms[a]['chain_id']
            a_segi = model.atoms[a]['segment_id']
            
            dist = N.sqrt( sum( (d_xyz - a_xyz)**2 ) )

            ## don't calculate angles within the same residue and 
            ##  for distances definately are not are h-bonds
            if dist < 3.0 and not\
                  ( d_nr == a_nr and d_cid == a_cid and d_segi == a_segi ):

                ## calculate angle for potenital hbond
                d_xyz_cov = xyzOfNearestCovalentNeighbour( d, model )
                a_xyz_cov = xyzOfNearestCovalentNeighbour( a, model )
                d_vec = d_xyz_cov - d_xyz
                a_vec = a_xyz - a_xyz_cov
                
                d_len = N.sqrt( sum( (d_vec)**2 ) )
                a_len = N.sqrt( sum( (a_vec)**2 ) )
                
                da_dot = N.dot( d_vec, a_vec)
                
                angle = 180 - N.arccos( da_dot / (d_len * a_len) )*180/N.pi
    
                if hbondCheck( angle, dist ):
                    hbond_lst += [[ d, a, dist, angle ]]
                    
    return hbond_lst
        

def hbondCheck( angle, length ):
    """
    A h-bond is longer that 2.4A but shorter that 3.5A, measuring
    from the nitrogen (donor) to the acceptor. The corresponding
    length from the hydrogen (donor) to the acceptor (used here)
    is 1.4 and 2.5A.
    
     - the optimal length is about 2.8A (1.8A measured from the hydrogen).
     - long h-bond (>3.6A) should be quite linear angle >150 
     - short h-bond (<3.6A) could be quite bent, but not more than 90-110

    @param angle: angle of bond to check
    @type  angle: float
    @param length: lenght of bond to check (D-H...A)
    @type  length: float

    @return: if the test is passed the cutoff angle for a
             hydrogen bond of the given length is returned, else None.
    @rtype: float OR None
    """
    ## h-bond has to be within the span 2.4 - 3.6 A
    ## here we will give a quite generous upper bound and lower
    if length < 1.4 or length > 2.9:
        return None

    ## two length angle pairs on a line that will give
    ## the cutoff value for allowed angles
    len_1, ang_1 = 1.5, 90  # 1.6, 90
    len_2, ang_2 = 2.8, 150 # 2.5, 150
    slope = (ang_2-ang_1)/(len_2-len_1)
    intersect = ang_1 - slope*len_1

    cutoff_angle = slope*length + intersect
    if angle < cutoff_angle:
        return None
    else:
        return cutoff_angle
        

def xyzOfNearestCovalentNeighbour( i, model ):
    """
    Closest atom in the same residue as atom with index i

    @param model: PDBModel 
    @type  model: PDBModel
    @param i: atom index 
    @type  i: int

    @return: coordinates of the nearest atom 
    @rtype: [float, float, float]
    """
    resModel = model.filter( residue_number=model.atoms[i]['residue_number'] )
    dist = N.sqrt( N.sum( (resModel.xyz - model.xyz[i])**2 , 1) )

    ## set distance to self to something high
    dist[ N.argmin(dist) ] = 100.
    
    pos_shortest =  N.nonzero( dist == min(dist) )[0]
 
    return resModel.xyz[ pos_shortest ]



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: icmCad values
        @rtype: [float]
        """
        from Biskit import PDBModel
        
        ## Loading PDB...
        m = PDBModel( T.testRoot() + '/lig/1A19.pdb' )
        m = m.compress( m.maskProtein() )

        hb = hbonds( m )

        xyz = xyzOfNearestCovalentNeighbour( 40, m )
        
        if local:
            print '\nThe nearest covalently attached atom to the'
            print '  atom with index 40 has the coordinates:'
            print xyz
    
            print 'Potential h-bonds in model:'
            print '(donor index, acceptor index, distance and angle)'
            for h in hb:
                print h
                
            globals().update( locals() )
                              
        return N.sum(N.ravel(hb[3:5])) + N.sum(xyz)


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: icmCad values
        @rtype:  [float]
        """
        return 2025.8997840075292 + 152.687011719

    
        
if __name__ == '__main__':

    test = Test()

    assert abs( test.run( local=1 ) - test.expected_result() ) < 1e-8




