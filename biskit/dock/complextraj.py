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
Trajectory of two proteins.
"""

from biskit.md import Trajectory, TrajError, EnsembleTraj
from biskit import hist
from biskit.dock import Complex as ProteinComplex

import biskit.core.oldnumeric as N0

import biskit.gnuplot as gnuplot


class ComplexTrajError( TrajError ):
    pass

class ComplexTraj( EnsembleTraj ):
    """
    Multi-copy trajectory of a protein-protein complex.
    Quite basic.

    ComplexTraj keeps the reference PDBModel of the normal Trajectory
    but, additionally, knows which chains in the reference belong to the
    receptor and which chains belong to the ligand. It tries to keep the
    2 lists of chain indices ( cr and cl ) up-to-date, even if chains are
    removed or re-ordered.
    
    Appending of atoms/chains (e.g. with concatAtoms) is currently not
    considered and will most likely leed to additional chains that are
    ignored by the complex-specific functions.
    """

    def __init__(self, traj=None, recChains=[0], n_members=10, rec=None):
        """
        Use::
           ComplexTraj( traj, recChains=[0], n_members=10 )
        
        @param traj: Trajectory (default: 0)
        @type  traj: Trajectory
        @param recChains: list of chain indices defining receptor
                          (default: [0])
        @type  recChains: [int]
        @param n_members: number of ensemble members (for multi-copy MD)
                          (default: 10)
        @type  n_members: int
        @param rec: instead of using recChains, autodetect rec chains
                    by comparing traj.ref with the model given as rec
                    (default: None)
        @type  rec: PDBModel
        """
        EnsembleTraj.__init__( self )

        if traj:
            self.__dict__.update( traj.__dict__ )

            if rec:
                self.cr = self.ref.compareChains( rec )[0]
            else:
                self.cr = N0.sort( recChains ).tolist()

            self.cl = list(range( 0, self.getRef().lenChains()))
            for c in self.cr:
                self.cl.remove(c)

        else:
            self.cr = self.cl = None


    def ligTraj( self ):
        """
        @return: ligand part of this Trajectory (copy)
        @rtype: EnsembleTraj
        """
        return self.takeChains( self.cl, EnsembleTraj )


    def recTraj( self ):
        """
        @return: receptor part of this Trajectory (copy)
        @rtype: EnsembleTraj
        """
        return self.takeChains( self.cr, EnsembleTraj )


    def refRec( self ):
        """
        @return: copy of the receptor part of the reference model
        @rtype: PDBModel
        """
        return self.getRef().takeChains( self.cr )


    def refLig( self ):
        """
        @return: copy of the ligand part of the reference model
        @rtype: PDBModel
        """
        return self.getRef().takeChains( self.cl )


    def refCom( self ):
        """
        @return: Complex
        @rtype: 
        """
        return ProteinComplex( self.refRec(), self.refLig() )


    def __getitem__( self, i ):
        return self.getComplex( i )


    def getComplex( self, index ):
        """
        Use::
          getComplex( frame_index ) -> Complex

        @param index: frame index
        @type  index: int

        @return: Complex
        @rtype: Complex 
        """
        m = self.getPDBModel( index )
        return ProteinComplex( m.takeChains(self.cr), m.takeChains(self.cl) )


    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.

        @param traj: Trajectory
        @type  traj: Trajectory
        
        @raise ComplexTrajError: if given traj is no ComplexTraj.
        """
        if not isinstance( traj, ComplexTraj ):
            raise ComplexTrajError(
                "Cannot replace ComplexTraj by normal Trajectory.")

        Trajectory.replaceContent( self, traj )
        self.cr = traj.cr
        self.cl = traj.cl


    def __translateChainIndices( self, atomIndices, newChainMap ):
        """
        Translate current chain indices into what they would look like in
        a PDBModel containing only the given atoms in the given order.
        
        @param atomIndices: indices of atoms
        @type  atomIndices: [int]
        @param newChainMap: chain map [0000011133333..]
        @type  newChainMap: [int]
        
        @return: { int:int, .. } map current chain indices to new ones
        @rtype: {int:int}
        
        @raise ComplexTrajError: if (parts of) chains are inserted into
                                 each other
        """
        ## todo: looks not very elegant

        oldChainMap = N0.take( self.ref.chainMap(), atomIndices )

        r = {}
        for i in range( len( oldChainMap ) ):
            old, new = oldChainMap[i], newChainMap[i]
            if old in r:
                if r[old] != new:
                    raise ComplexTrajError(
                        "Can't insert different chains into each other.")
            else:
                r[old] = new

        return r


    def takeAtoms( self, indices, returnClass=None ):
        """
        takeAtoms( indices, [returnClass] ) -> ComplexTraj
    
        @param indices: atoms to extract
        @type  indices: [int]
        @param returnClass: return type (default: current class)
        @type  returnClass: class

        @return: ComplexTraj
        @rtype: ComplexTraj
        
        @raise ComplexTrajError: if (parts of) chains are inserted into
                                 each other
        """
        r = Trajectory.takeAtoms( self, indices, returnClass )

        oldToNew = self.__translateChainIndices( indices, r.ref.chainMap() )

        r.cr = [ oldToNew[ c ] for c in self.cr if c in oldToNew ]
        r.cl = [ oldToNew[ c ] for c in self.cl if c in oldToNew ]

        return r


    ## NOT TESTED!!
    def concatAtoms( self,  recChains, *traj ):
        """
        Concatenate 2 trajectories of same (frame) length 'horizontally', i.e.
        for each frame the coordinates of one are appended to the coordinates
        of the other. The ref model of the new trajectory is a 'semi-deep' copy
        of this trajectory's model (see L{Biskit.PDBModel.take()} ).::
          concatAtoms( recChains, traj1 [traj2, traj3..]) -> ComplexTraj

          
        @param recChains: chain indices into traj that will be assigned to the
                          receptor; all remaining chains (if any) go to the
                          ligand
        @type  recChains: [int]
        @param traj: trajectories to concatenate     
        @type  traj: Trajectory

        @return: ComplexTraj
        @rtype: ComplexTraj        
        
        @warning: NOT TESTED!!
        """
        currentLen = self.ref.lenChains()

        recChains.sort()
        cl = [c for c in range(0, traj.ref.lenChains() ) if not c in recChains]

        recChains = [ c + currentLen for c in recChains ]
        cl        = [ c + currentLen for c in cl ]

        r = EnsembleTraj.concatAtoms( self, *traj )
        r.cr = self.cr + recChains
        r.cl = self.cl + cl

        return r


    def atomContacts( self, index, cutoff=4.5, rec_mask=None, lig_mask=None ):
        """
        Use::
          atomContacts( frame_index ) -> array len_rec x len_lig of 0||1

        @param index: frame index
        @type  index: [int]
        @param cutoff: contact cutoff in Angstrom (default: 4.5)
        @type  cutoff: float
        @param rec_mask: atom mask
        @type  rec_mask: [int]
        @param lig_mask: atom mask
        @type  lig_mask: [int]

        @return: contact matrix 
        @rtype: matrix
        """
        return self[ index ].atomContacts( cutoff, rec_mask, lig_mask )


    def averageContacts( self, step=10, cutoff=4.5 ):
        """
        Use::
          averageContacts( step=1, cutoff=4.5 )
        
        @param step: take only each |step|th frame (default: 10)
        @type  step: int
        @param cutoff: distance cutoff in Angstrom (default: 4.5)
        @type  cutoff: float
        
        @return: contact matrix with frequency of each contact in
                 (thinned) traj.
        @rtype: matrix
        """
        r = [ self.atomContacts( i, cutoff=cutoff ) 
             for i in range(0, len(self), step ) ]
        return N0.sum( N0.array( r ) ) / ( 1. * len(r) )


    def plotContactDensity( self, step=1, cutoff=4.5 ):
        """
        Example. plot histogramm of contact density. Somehing wrong??

        @raise ComplexTrajError: if gnuplot program is not installed
        """
        if not gnuplot.installed:
            raise ComplexTrajError('gnuplot program is not installed')
        r = self.averageContacts( step, cutoff )
        r = N0.ravel( r )
        r = N0.compress( r, r )
        gnuplot.plot( hist.density( r, 10 ) )



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    TAGS = [ BT.LONG ]

    def test_ComplexTraj(self):
        """Dock.ComplexTraj test"""

        import biskit.tools as T

        ## there is no complex trajectory in the test folder so will have
        ## to create a fake trajectory with a complex
        f =  [ T.testRoot()+ '/com/1BGS.pdb' ] * 5
        t = Trajectory( f, verbose=self.local )

        t = ComplexTraj( t, recChains=[0] )

        #if self.local:
            #print 'plotting contact density...'
            #t.plotContactDensity( step=2 )

        ## create a fake second chain in the ligand
        for i in range( 1093+98, 1968 ):
            t.ref.atoms['chain_id'][i] = 'B'

        t.ref.chainIndex( force=1, cache=1 )
        t.cl = [1,2]

        r = N0.concatenate((list(range(1093,1191)), list(range(0,1093)), list(range(1191,1968))))

        tt = t.takeAtoms( r )

        contactMat = tt.atomContacts( 1 )
        
        if self.local:
            print('Receptor chains: %s    Ligand chains: %s'%(t.cr, t.cl))
            
        self.assertEqual( N0.sum(N0.ravel(contactMat)), 308 )

if __name__ == '__main__':

    #import biskit.tools as T

    ### there is no complex trajectory in the test folder so will have
    ### to create a fake trajectory with a complex
    #f =  [ T.testRoot()+ '/com/1BGS.pdb' ] * 5
    #t = Trajectory( f, verbose=1 )

    #t = ComplexTraj( t, recChains=[0] )

    ##if self.local:
        ##print 'plotting contact density...'
        ##t.plotContactDensity( step=2 )

    ##for i in range( 1093+98, 1968 ):
    #t.ref.atoms['chain_id'][1093+98:1968] = 'B'
    #t.ref.chainIndex( force=1, cache=1 )
    #t.cl = [1,2]

    #r = N0.concatenate((range(1093,1191), range(0,1093), range(1191,1968)))

    #tt = t.takeAtoms( r )

    #contactMat = tt.atomContacts( 1 )
    BT.localTest()
