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
This is a helper class for ComplexList.
"""

import biskit
from biskit.errors import BiskitError
from biskit import LocalPath, PDBModel
import biskit.tools as T

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.dock; __package__ = "biskit.dock"

from .complex import Complex

class RegistryError( BiskitError ):
    pass

class ComplexModelRegistry:
    """
    This is a helper class for ComplexList.

    Keep unique copies of the rec and lig models from many Complexes.
    Make sure that 2 Complexes with the same rec_model (same by file
    name and unchanged) always point to the same PDBModel instance.
    """

    def __init__( self ):

        self.rec_f2model = {}
        self.lig_f2model = {}

        self.rec_f2com = {}
        self.lig_f2com = {}

        self.initVersion = biskit.__version__


    def addComplex( self, com ):
        """
        Register Complex with the registry.
        
        @param com: complex
        @type  com: Complex
        """
        com.rec_model,fr = self.__sync_model( com.rec_model, self.rec_f2model )
        com.lig_model,fl = self.__sync_model( com.lig_model, self.lig_f2model )

        ## optimized for speed, that's why a bit awkward
        if fr is not None:
            coms = self.rec_f2com.get( fr, None )
            if coms is None:
                self.rec_f2com[ fr ] = [ com ]
            else:
                coms.append( com )

        if fl is not None:
            coms = self.lig_f2com.get( fl, None )
            if coms is None:
                self.lig_f2com[ fl ] = [ com ]
            else:
                coms.append( com )


    def removeComplex( self, com ):
        """
        Remove a Complex from the registry.
        
        @param com: complex
        @type  com: Complex
        """        
        self.__removeModel(com.rec_model,com, self.rec_f2model, self.rec_f2com)
        self.__removeModel(com.lig_model,com, self.lig_f2model, self.lig_f2com)


    def __removeModel( self, model, com, f2model, f2com ):
        """
        Remove model of a complex from the registry.
        
        @param model: receptor or ligand model
        @type  model: PDBModel
        @param com: complex
        @type  com: Complex
        @param f2model: dictionary mapping files to models
        @type  f2model: {str:PDBModel}
        @param f2com: dictionary mapping files to complexes
        @type  f2com: {str:Complex}
        """        
        ## optimized for speed
        f = model.source
        coms = f2com[ f ]
        coms.remove( com )

        if len( coms ) == 0:
            del f2com[ f ]
            del f2model[ f ]


##     def update( self, otherReg ):
##         """
##         Add content of another registry to this one, make all
##         Complex.rec_model / lig_model to point to this Registry.
##         """
##         ## currently not used
##         for k,v in otherReg.rec_f2model:

##             if not k in self.rec_f2model:
##                 self.rec_f2model[ k ] = v
##             else:

##                 m = self.rec_f2model[ k ]

##                 if otherReg.rec_f2model != m:
##                     for c in otherReg.rec_f2com[ k ]:
##                         c.rec_model = m

##                 for c in otherReg.rec_f2com[ k ]:
##                     if not c in self.rec_f2com[ k ]:
##                         self.rec_f2com[ k ] + [c]

##         for k,v in otherReg.lig_f2model:

##             if not k in self.lig_f2model:
##                 self.lig_f2model[ k ] = v

##             else:
##                 m = self.lig_f2model[ k ]

##                 if otherReg.lig_f2model != m:
##                     for c in otherReg.lig_f2com[ k ]:
##                         c.lig_model = m

##                 for c in otherReg.lig_f2com[ k ]:
##                     if not c in self.lig_f2com[ k ]:
##                         self.lig_f2com[ k ] + [c]


    def getRecModel( self, source ):
        """
        Get receptor model belonging to source.
        
        @param source: path
        @type  source: str
        
        @return: model
        @rtype: PDBModel       
        """
        return self.rec_f2model[source]


    def getLigModel( self, source ):
        """
        Get ligand model belonging to source.
        
        @param source: path
        @type  source: str
        
        @return: model
        @rtype: PDBModel       
        """        
        return self.lig_f2model[source]


    def getModel( self, source):
        """
        Get model belonging to source.
        
        @param source: path
        @type  source: str
        
        @return: model
        @rtype: PDBModel       
        """        
        if source in self.rec_f2model:
            return self.rec_f2model[ source ]

        return self.lig_f2model[ source ]


    def recModels( self ):
        """
        Get a list with all receptor models.

        @return: list of models
        @rtype: [PDBModel]       
        """              
        return list(self.rec_f2model.values())


    def ligModels( self ):
        """
        Get a list with all ligand models.

        @return: list of models
        @rtype: [PDBModel]       
        """                 
        return list(self.lig_f2model.values())


##     def getSubRegistry( self, cl ):
##         """
##         Get fraction of this Registry that is used by given list of
##         complexes.
##         cl - [ Complex ], all complexes must be part of the registry
##         -> ComplexModelRegistry
##         """
##         if isinstance( cl, ComplexList ):
##             return cl.models

##         r = self.__class__()

##         for c in cl:
##             r.addComplex( c )

##         return r


    def getRecComplexes( self, model ):
        """
        Get the complexes of which a given receptor model is a component.

        @param model: LocalPath or PDBModel
        @type  model: object
        
        @return: list of Complexes
        @rtype: [Complex]
        """
        return self.__getComplexes( model, self.rec_f2com )


    def getLigComplexes( self, model ):
        """
        Get the complexes of which a given ligand model is a component.
        
        @param model: LocalPath or PDBModel
        @type  model: object
        
        @return: list of Complexes
        @rtype: [Complex]
        """
        return self.__getComplexes( model, self.lig_f2com )


    def __getComplexes( self, model, f2com ):
        """
        Get the complexes of which a given ligand model is a component.
        
        @param model: LocalPath or PDBModel
        @type  model: object
        @param f2com: dictionary mapping paths to complexes
        @type  f2com: {str:Complex}
        
        @return: list of Complexes
        @rtype: [Complex]

        @raise RegistryError: if model has no source file.
        """
        if isinstance( model, LocalPath):
            f = model

        try:
            if isinstance( model, PDBModel ):
                f = model.source

            return f2com[ f ]

        except KeyError as why:
            raise RegistryError( "Model with source '%r' is unknown" % f )


    def __sync_model( self, m, f2model ):
        """
        Get a shared model instance that is equal to m. If there is no such
        instance in the registry, m is returned and m is added to the
        registry - but only if it has been pickled to disc and hasn't
        changed since.
        
        @param m: PDBModel
        @type  m: PDBModel
        @param f2model: dictionary with the path to the file as key
        @type  f2model: { str:PDBModel }
        
        @return: Model from f2model equivalent to m, otherwise add m to f2model
                 and return m. File name or None for stray model
        @rtype: PDBModel, str
        """
        ## The problem here is to minimize the calls to PDBModel.validSource()
        ## because 60.000 calls to os.path.exists() would take a lot of time
        assert isinstance( m, PDBModel )

        if isinstance( m.source, LocalPath ):
            f = m.source
        else:
            f = None

        if f is not None: ## don't use "if f" because that would call f.__len__
            if m.xyzChanged: # or m.atomsChanged:
                ## don't add it to registry, it will remain a stray model
                return m, f

            found = f2model.get( f, None )
            if found is not None:
                ## the source has already been checked, we can skip the usual..
                ## .. exists()
                return found, f

            ## potentially new entry, check also whether the file exists
            if m.source.exists():
                f2model[ f ] = m

        return m, f


    def __str__( self ):
        s = "Receptor models:"
        for f in self.rec_f2model:
            s += "\n%s: %i complexes" % \
                 (f.formatted(), len( self.rec_f2com[ f ] ) )

        s += "\nLigand models:"
        for f in self.lig_f2model:
            s += "\n%s: %i complexes" % \
                 (f.formatted(), len( self.lig_f2com[ f ] ) )

        return s


    def __repr__( self ):
        return "ComplexModelRegistry\n" + self.__str__()


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_ComplexModelRegistry(self):
        """Dock.ComplexModelRegistry test"""
        from biskit.dock import ComplexList
        
        self.cl = T.load( T.testRoot() +'/dock/hex/complexes.cl' )
        self.cl = self.cl.toList()
        
        self.r = ComplexModelRegistry()

        for c in self.cl[:500]:
            self.r.addComplex( c )

        check = self.r.getLigComplexes( self.r.ligModels()[0] )

        self.assertEqual( len(check), 500 )
    

if __name__ == '__main__':

    BT.localTest()
