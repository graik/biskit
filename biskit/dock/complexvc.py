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
Complex that keeps track of its previous versions (conformations).
"""
import numpy as N
import biskit.mathUtils as MU
from biskit import EHandler
from biskit.dock import Complex as ProtComplex
from biskit.dock import ComplexList


class ComplexVC( ProtComplex ):
    """
    Complex that keeps track of its previous versions. The object behaves
    like a normal L{Biskit.Dock.Complex} but encapsulates a
    L{Biskit.Dock.ComplexList} with older versions/conformations of
    this Complex.

      - com[0] gives the very first version, com[1] the second, etc.
      - com[-1] gives the current version (but as normal Complex, not
        ComplexVC).
      - e.g. com[0, 'fnac'] gives the according info value of the first version
      - e.g. com['fnac']    gives the current info value (as for Dock.Complex)
    """

    def __init__(self, rec_model, lig_model, com_0=None,
                 ligMatrix=None, info={} ):
        """
        Create a new ComplexVC from a previous Complex or ComplexVC
        and a new set of receptor, ligand conformation and transformation.

        @param rec_model: PDBModel/XplorModel, receptor conformation
        @type  rec_model: PDBModel
        @param lig_model: PDBModel/XplorModel, ligand conformation
        @type  lig_model: PDBModel
        @param com_0: Complex /ComplexVC, previous version(s) of this com
        @type  com_0: Complex OR ComplexVC
        @param ligMatrix: transformation matrix of ligand versus receptor
        @type  ligMatrix: 4x4 array
        @param info: info dictionary {'info_key': value, ..}, additional infos
        @type  info: dict
        """
        ProtComplex.__init__(self, rec_model, lig_model, ligMatrix, info )

        if isinstance( com_0, ComplexVC ):
            ## get history from mother complex ...
            self.history = com_0.history

            com_0 = ProtComplex( com_0.rec_model, com_0.lig_model,
                                 com_0.ligandMatrix, com_0.info )

        else:
            ## ... or create a new history
            self.history = ComplexList()

        if com_0 is not None:
            self.history.append( com_0 )

            ## save only differences between old and new conformations
            self.rec_model = self.__syncModel( self.rec_model, com_0.rec_model)
            self.lig_model = self.__syncModel( self.lig_model, com_0.lig_model)


    def __iter__(self):
        """
        __iter__() <==> for k in self
        """
        return iter( self.history + [self] )


    def __len__( self ):
        """
        length of self
        """
        return len( self.history ) + 1


    def getComplex( self, i, copy=0 ):
        """
        Return one version as Complex. -1 returns a copy of the latest
        version. By default the info dictionary remains connected but
        other fields don't. I.e. replacing rec_model in the copy does
        not affect the original complex.

        @param i: index in complex history, -1 returns toComplex()
        @type  i: int
        @param copy: copy info dictionary in case of i==-1 (changes in
                     c.getComplex( -1 ).info will not appear in c.info [0]
        @type  copy: 1|0

        @return: complex
        @rtype: Dock.Complex
        """
        if i == -1:
            return self.toComplex( copy=copy )

        if i >= 0 and i < len( self.history ):
            return self.history[ i ]

        return self.toSimpleList()[i]


    def __getitem__( self, k ):
        """
        Get a Complex OR a info key value for a Complex specified Complex
        OR info dic value.

        @param k: tuple OR int OR str
        @type  k: any OR Complex
        """
        if type(k) == tuple:

            i, key = k[0], k[1]

            if i == len( self.history ):
                return self.info[key]

            return self.getComplex(i).info[key]

        if isinstance(k, (int, N.integer)):

            return self.getComplex( k )

        return self.info[k]


    def __syncModel( self, new_model, old_model ):
        """
        Connect new rec or lig model to old one, to minimize storage.

        @param new_model: PDBModel / XplorModel
        @type  new_model: PDBModel
        @param old_model: PDBModel / XplorModel
        @type  old_model: PDBModel

        @return: PDBModel / XplorModel, new model that only keeps
                 changes relative to old, the old model becomes the
                 source of the new, if possible
        @rtype: PDBModel
        """
        ## try to fix atom order of new_model so that it is identical to old
        if old_model.equals( new_model ) != [1,1]:
            i_new, i_old = new_model.compareAtoms( old_model )

            if len( i_new ) == len( new_model ):
                new_model.keep( i_new )

        ## create result model that only keeps difference of new and old
        if old_model.equals( new_model ) == [1,1]:

            ## stays compatible with XplorModel.__init__ and PDBModel.__init
            r = old_model.__class__( source=old_model )

            r.setXyz( new_model.getXyz() )

            ## check for profiles identical to source and adapt 'changed'
            r.update()  

            if not MU.arrayEqual( r.xyz, old_model.xyz ):
                r.removeProfile( 'relASA', 'ASA_sc', 'ASA_total', 'ASA_bb' )

            return r

        EHandler.warning(
            'ComplexVC: Cannot connect new to old PDBModel.')

        new_model.disconnect()
        return new_model


    def sortHistory( self ):
        """
        Sort by date
        """
        self.history = self.history.sortBy( 'date' )


    def toList( self ):
        """
        @return: all historic complexes plus current one as Complex
        @rtype: ComplexList
        """
        return self.history + [ self.toComplex() ]


    def toSimpleList( self ):
        """
        @return: all historic complexes plus current one as Complex
        @rtype: [Complex]
        """
        return self.history.toList() + [ self.toComplex() ]


    def toComplex( self, copy=0 ):
        """
        Copy of latest version as a normal Complex.

        @param copy: also disconnect info dict (default: 0)
        @type  copy: 1|0

        @return: Complex
        @rtype: Complex
        """
        r = ProtComplex( self.rec_model, self.lig_model,
                         self.ligandMatrix, self.info )
        if not copy:
            r.info = self.info

        return r


    def valuesOf( self, infoKey, default=None ):
        """
        Get info values from all versions of this complex (oldest first).

        @param infoKey: info dic key
        @type  infoKey: str
        @param default: default value, if key is not present
        @type  default: any

        @return: list of values
        @rtype: [any]
        """
        return [ c.get( infoKey, default ) for c in self ]



#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_ComplexVC(self):
        """Dock.ComplexVC test"""
        import time
        import biskit.tools as T

        c = T.load( T.testRoot() + '/com/ref.complex' )

        self.ce= ComplexVC( c.rec_model, c.lig(), c,
                                  info={'comment':'test'} )

        time.sleep( 2 )

        lig = self.ce.lig().transform( MU.randomRotation(), [0,0,0] )
        self.ce2 = ComplexVC( self.ce.rec_model, lig, self.ce,
                                    info={'comment':'test2'})

        if self.local:
            print('\nFound %i versions of the complex: ' % len(self.ce2))
            for x in self.ce2:
                print('\t* ' + x['date'])

            print('Comments: ', self.ce2.valuesOf('comment'))

        self.assertEqual( self.ce2.valuesOf('comment'),
                          [None, 'test', 'test2'])


if __name__ == '__main__':

    BT.localTest()

