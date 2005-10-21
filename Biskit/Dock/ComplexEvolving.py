##
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Author$
## last $Date$

import Biskit.tools as t
import Biskit.mathUtils as MU
from Biskit import EHandler
from Biskit.Dock.Complex import Complex as ProtComplex
from ComplexList import ComplexList, ComplexListError, ConditionError
     

class ComplexEvolving( ProtComplex ):
    """
    Complex that keeps track of its previous versions. The object behaves
    like a normal Biskit.Complex but encapsulates a ComplexList with older
    versions/conformations of this Complex.

    * com[0] gives the very first version, com[1] the second, etc.
    * com[-1] gives the current version (but as normal Complex, not
      ComplexEvolving).
    * e.g. com[0, 'fnac'] gives the according info value of the first version
    * e.g. com['fnac']    gives the current info value (as for Dock.Complex)
    """

    def __init__(self, rec_model, lig_model, com_0=None,
                 ligMatrix=None, info={} ):
        """
        Create a new ComplexEvolving from a previous Complex or ComplexEvolving
        and a new set of receptor, ligand conformation and transformation.
        rec_model  - PDBModel/PCRModel, receptor conformation
        lig_model  - PDBModel/PCRModel, ligand conformation
        com_0      - Complex /ComplexEvolving, previous version(s) of this com
        ligMatrix  - 4x4 array, transformation matrix of ligand versus receptor
        info       - {'info_key': value, ..}, additional infos
        """

        ProtComplex.__init__(self, rec_model, lig_model, ligMatrix, info )

        if isinstance( com_0, ComplexEvolving ):
            ## get history from mother complex ...
            self.history = com_0.history

            com_0 = ProtComplex( com_0.rec_model, com_0.lig_model,
                                 com_0.ligandMatrix, com_0.info )

        else:
            ## ... or create a new history
            self.history = ComplexList()

        if com_0 != None:
            self.history.append( com_0 )

            ## save only differences between old and new conformations
            self.rec_model = self.__syncModel( self.rec_model, com_0.rec_model)
            self.lig_model = self.__syncModel( self.lig_model, com_0.lig_model)


    def version( self ):
        return 'ComplexEvolving $Revision$'

    def __iter__(self):
        """__iter__() <==> for k in self"""
        return iter( self.history + [self] )

    def __len__( self ):
        return len( self.history ) + 1


    def getComplex( self, i, copy=0 ):
        """
        Return one version as Complex. -1 returns a copy of the latest
        version. By default the info dictionary remains connected but
        other fields don't. I.e. replacing rec_model in the copy does
        not affect the original complex.
        i    - int, index in complex history, -1 returns toComplex()
        copy - 1||0, copy info dictionary in case of i==-1 (changes in
               c.getComplex( -1 ).info will not appear in c.info [0]
        -> Dock.Complex
        """
        if i == -1:
            return self.toComplex( copy=copy )
        
        if i >= 0 and i < len( self.history ):
            return self.history[ i ]

        return self.toSimpleList()[i]
        

    def __getitem__( self, k ):
        if type(k) == tuple:

            i, key = k[0], k[1]

            if i == len( self.history ):
                return self.info[key]
            
            return self.getComplex(i).info[key]

        if type(k) == int:

            return self.getComplex( k )
        
        return self.info[k]


    def __syncModel( self, new_model, old_model ):
        """
        Connect new rec or lig model to old one, to minimize storage.
        new_model - PDBModel / PCRModel
        old_model - PDBModel / PCRModel
        -> PDBModel / PCRModel, new model that only keeps changes relative
                                to old,
        the old model becomes the source of the new, if possible
        """
        ## try to fix atom order of new_model so that it is identical to old
        if old_model.equals( new_model ) != [1,1]:
            i_new, i_old = new_model.compareAtoms( old_model )

            if len( i_new ) == len( new_model ):
                new_model.keep( i_new )

        ## create result model that only keeps difference of new and old
        if old_model.equals( new_model ) == [1,1]:

            ## stays compatible with PCRModel.__init__ and PDBModel.__init
            r = old_model.__class__( source=old_model )
            
            r.setXyz( new_model.getXyz() )

            ## check for profiles identical to source and adapt 'changed'
            r.update()  

            if not MU.arrayEqual( r.xyz, old_model.xyz ):
                r.removeProfile( 'relASA', 'ASA_sc', 'ASA_total', 'ASA_bb' )

            return r

        EHandler.warning(
            'ComplexEvolving: Cannot connect new to old PDBModel.')

        new_model.disconnect()
        return new_model


    def sortHistory( self ):
        self.history = self.history.sortBy( 'date' )


    def toList( self ):
        """
        -> ComplexList, all historic complexes plus current one as Complex
        """
        return self.history + [ self.toComplex() ]

    def toSimpleList( self ):
        """
        -> [ Complex ], all historic complexes plus current one as Complex
        """
        return self.history.toList() + [ self.toComplex() ]


    def toComplex( self, copy=0 ):
        """
        Copy of latest version as a normal Complex.
        copy - 1||0, also disconnect info dict [0]
        -> Complex
        """
        r = ProtComplex( self.rec_model, self.lig_model,
                         self.ligandMatrix, self.info )
        if not copy:
            r.info = self.info

        return r

    def valuesOf( self, infoKey, default=None ):
        """
        Get info values from all versions of this complex (oldest first).
        infoKey - str
        default - any, default value, if key is not present
        -> [ any ]
        """
        return [ c.get( infoKey, default ) for c in self ]
        

### TEST ###

import time

if __name__ == '__main__':

##     from Biskit import *
##     from Biskit.Dock import *

    c = t.Load( t.testRoot() + '/com_wet/ref.complex')

    ce= ComplexEvolving( c.rec_model, c.lig(), c, info={'comment':'test'} )

    time.sleep( 2 )

    lig = ce.lig().transform( MU.randomRotation(), [0,0,0] )
    ce2 = ComplexEvolving( ce.rec_model, lig, ce,
                           info={'comment':'test2'})

    for x in ce2:
        print x['date']
