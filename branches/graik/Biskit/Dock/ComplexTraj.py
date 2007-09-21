##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$
"""
Trajectory of two proteins.
"""

from Biskit import Trajectory, TrajError, EnsembleTraj, hist
from Complex import Complex as ProteinComplex

import Numeric as N

import Biskit.gnuplot as gnuplot
    

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

    def __init__(self, traj=None, recChains=[0], n_members=10, rec=None ):
        """
        ComplexTraj( traj, recChains=[0], n_members=10 )
        traj       - Trajectory
        recChains  - [int], list of chain indices defining receptor
        n_members  - int, number of ensemble members (for multi-copy MD)
        rec        - PDBModel, instead of using recChains, autodetect rec
                     chains by comparing traj.ref with the model given as rec
        """

        EnsembleTraj.__init__( self )

        if traj:
            self.__dict__.update( traj.__dict__ )

            if rec:
                self.cr = self.ref.compareChains( rec )[0]
            else:
                self.cr = N.sort( recChains ).tolist()

            self.cl = range( 0, self.getRef().lenChains() )
            for c in self.cr:
                self.cl.remove(c)

        else:
            self.cr = self.cl = None


    def version( self ):
        return EnsembleTraj.version(self) + '; ComplexTraj $Revision$'

    def ligTraj( self ):
        """-> EnsembleTraj, ligand part of this Trajectory (copy)"""
        return self.takeChains( self.cl, EnsembleTraj )

    def recTraj( self ):
        """-> EnsembleTraj, receptor part of this Trajectory (copy)"""
        return self.takeChains( self.cr, EnsembleTraj )

    def refRec( self ):
        """-> PDBModel, copy of the receptor part of the reference model """
        return self.getRef().takeChains( self.cr )

    def refLig( self ):
        """-> PDBModel, copy of the ligand part of the reference model """
        return self.getRef().takeChains( self.cl )

    def refCom( self ):
        """-> Complex """
        return ProteinComplex( self.refRec(), self.refLig() )

    def __getitem__( self, i ):
        return self.getComplex( i )

    def getComplex( self, index ):
        """getComplex( frame_index ) -> Complex """
        m = self.getPDBModel( index )
        return ProteinComplex( m.takeChains(self.cr), m.takeChains(self.cl) )

        
    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.
        !! ComplexTrajError, if given traj is no ComplexTraj.
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
        atomIndices  - [int]
        newChainMap  - [0000011133333..]
        -> { int:int, .. } map current chain indices to new ones
        !! ComplexTrajError, if (parts of) chains are inserted into each other
        """
        ## todo: looks not very elegant
        
        oldChainMap = N.take( self.ref.chainMap(), atomIndices )

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
        The 
        indices     - [int], atoms to extract
        returnClass - class, return type (default: current class)
        !! ComplexTrajError, if (parts of) chains are inserted into each other
        """
        r = Trajectory.takeAtoms( self, indices, returnClass )

        oldToNew = self.__translateChainIndices( indices, r.ref.chainMap() )

        r.cr = [ oldToNew[ c ] for c in self.cr if c in oldToNew ]
        r.cl = [ oldToNew[ c ] for c in self.cl if c in oldToNew ]

        return r


    ## NOT TESTED!!
    def concatAtoms( self,  recChains, *traj ):
        """
        concatAtoms( recChains, traj1 [traj2, traj3..]) -> ComplexTraj
        Concatenate 2 trajectories of same (frame) length 'horizontally', i.e.
        for each frame the coordinates of one are appended to the coordinates
        of the other. The ref model of the new trajectory is a 'semi-deep' copy
        of this trajectory's model (see PDBModel.take() ).
        recChains - [int], chain indices into traj that will be assigned to the
                    receptor; all remaining chains (if any) go to the ligand
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
        """atomContacts( frame_index ) -> array len_rec x len_lig of 0||1"""
        return self[ index ].atomContacts( cutoff, rec_mask, lig_mask )


    def averageContacts( self, step=10, cutoff=4.5 ):
        """
        averageContacts( step=1, cutoff=4.5 )
        step   - int, take only each |step|th frame
        cutoff - float, distance cutoff in A
        -> contact matrix with frequency of each contact in (thinned) traj.
        """
        r = [ self.atomContacts( i ) for i in range(0, len(self), step ) ]
        return N.sum( N.array( r ) ) / ( 1. * len(r) )


    def plotContactDensity( self, step=1, cutoff=4.5 ):
        """Example. plot histogramm of contact density. Somehing wrong??"""

        if not gnuplot.installed:
            raise ComplexTrajError, 'gnuplot program is not installed'
        r = self.averageContacts( step, cutoff )
        r = N.ravel( r )
        r = N.compress( r, r )
        gnuplot.plot( hist.density( r, 10 ) )
        

if __name__ == '__main__':

    from tools import *
    
    t = Load( testRoot()+'/com_pc2_00/traj.dat' )

    t = ComplexTraj( t, [0] )

##    t.plotContactDensity( step=2 )
    for i in range( 1093+98, 1968 ):
        t.ref.atoms[i]['chain_id'] = 'C'

    t.cl = [1,2]

    r = N.concatenate((range(1093,1191), range(0,1093), range(1191,1968)))

    tt = t.takeAtoms( r )