## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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
Store and manipulate coordinates and atom information.
"""

import biskit.tools as T
from biskit import molUtils
from biskit import mathUtils
from biskit import match2seq
from biskit import rmsFit
from biskit.core.localpath import LocalPath
from biskit.errors import BiskitError
from biskit.profileCollection import ProfileCollection, ProfileError
from biskit.core.pdbparserFactory import PDBParserFactory
from biskit.core.pdbparseFile import PDBParseFile
from biskit import biounit as BU
from biskit.core import oldnumeric as N0
from biskit.core.scientificIO import PDB as IO
from biskit import EHandler

from biskit.future import Residue

import biskit as B

import numpy as N

import os, sys
import copy
import time
import string
import types
import functools


class PDBProfiles( ProfileCollection ):
    """
    A ProfileCollection that triggers an update() of its parent :class:`PDBModel`
    if an empty or (optionally) missing profile is requested.

    .. seealso:: `biskit.ProfileCollection`
    """

    def __init__(self, model=None, profiles=None, infos=None ):
        """
        :param model: parent model of this ProfileCollection
        :type  model: PDBModel
        :param profiles: dictionary of existing profiles
        :type  profiles: { 'name' : list/array }
        :param infos:    dictionary of existing meta infos
        :type  infos:    { 'name' : { 'date' : ... } }
        """
        ProfileCollection.__init__( self, profiles=profiles, infos=infos )
        self.model = model
    
    def clone(self):
        """
        Include reference to (same) parent model.
        """
        return self.__class__(self.model, copy.deepcopy(self.profiles),
                              copy.deepcopy(self.infos))


## should work but leads to some kind of loop condition
##    def profLength(self, default=0):
##        r = ProfileCollection.profLength(self, default=None)
##        if r is not None:
##            return r
##        
##        if self.model and self.model.xyz is not None:
##            return len(self.model.xyz)
##        
##        return default
        

    def get( self,  name, default=None, update=True, updateMissing=False ):
        """
        Fetch a profile::
          get( profKey, [default] ) -> list or array of values 
        Or:
          get( (profKey, infoKey), [default] ) -> single value of metainfo dict

        This method extends the standard :class:`ProfileCollection.get` by the
        ability to fetch empty (None) or missing profiles from a source
        file or model.

        :param name: profile key or profile and info key
        :type  name: str OR (str, str)
        :param default: default result if no profile is found,
                        if None and no profile is found, attempt update
        :type  default: any
        :param update: update from source before returning empty profile [0]
        :type  update: bool
        :param updateMissing: update from source before reporting missing
                              profile [0]
        :type  updateMissing: bool

        :raise ProfileError: if no profile is found with |name|
        """
        try:
            r = ProfileCollection.get( self, name, default=default)

        except ProfileError as e:
            if updateMissing:
                r = None
            else:
                raise ProfileError(e)

        if r is None and (update or updateMissing):

            ## only read PDB source if this is indeed useful
            if not name in PDBModel.PDB_KEYS and \
               PDBParseFile.supports( self.model.validSource() ):
                return None

            self.model.update( updateMissing=updateMissing )

            ## try again
            r = ProfileCollection.get(self, name )

        return r

class PDBResidueProfiles(PDBProfiles):
    """Work in progress -- give residue profiles a default length"""
    
    def profLength(self, default=0):
        r = ProfileCollection.profLength(self, default=None)
        if r is not None:
            return r
        
        if self.model and self.model._resIndex is not None:
            return len(self.model._resIndex)
        
        return default


class PDBError(BiskitError):
    pass

class PDBIndexError( PDBError):
    """Errors warning of issues with residue or chain index."""
    pass

class PDBModel:
    """
    Store and manipulate coordinates and atom infos stemming from a PDB file.
    Coordinates are stored in the numpy array 'xyz'; the additional atom infos
    from the PDB (name, residue_name, and many more) are efficiently stored in
    a :class:`PDBProfiles` instance 'atoms' which can be used to also associate
    arbitrary other data to the atoms. Moreover, a similar collection
    'residues' can hold data associated to residues (but is initially empty).
    A normal dictionary 'info' accepts any information about the whole model.

    For detailed documentation,
    see http://biskit.pasteur.fr/doc/handling_structures/PDBModel

    @todo:
       * outsource validSource into PDBParserFactory
       * prevent repeated loading of test PDB for each test
    """

    #: keys of all atom profiles that are read directly from the PDB file
    PDB_KEYS = ['name', 'residue_number', 'insertion_code', 'alternate',
                'name_original', 'chain_id', 'occupancy', 'element',
                'segment_id', 'charge', 'residue_name', 'after_ter',
                'serial_number', 'type', 'temperature_factor']

    def __init__( self, source=None, pdbCode=None, noxyz=0, skipRes=None,
                  headPatterns=[] ):
        """
        Examples:

        - `PDBModel()` creates an empty Model to which coordinates (field xyz)
          and PDB records (atom profiles) have still to be added.
        - `PDBModel( file_name )` creates a complete model with coordinates
          and PDB records from file_name (pdb, pdb.gz, or pickled PDBModel)
        - `PDBModel( PDBModel )` creates a copy of the given model
        - `PDBModel( PDBModel, noxyz=1 )` creates a copy without coordinates

        :param source: str, file name of pdb/pdb.gz file OR pickled PDBModel OR
                      PDBModel, template structure to copy atoms/xyz field from
        :type  source: str or PDBModel
        :param pdbCode: PDB code, is extracted from file name otherwise
        :type  pdbCode: str or None
        :param noxyz: 0 (default) || 1, create without coordinates
        :type  noxyz: 0||1
        :param headPatterns: [(putIntoKey, regex)] extract given REMARK values
        :type  headPatterns: [(str, str)]

        :raise PDBError: if file exists but can't be read
        """
        self.source = source
        if type( source ) is str and ( len(source) != 4 or \
                                       os.path.isfile( source ) ):
            self.source = LocalPath( source )

        self.__validSource = 0
        self.fileName = None
        self.pdbCode = pdbCode
        self.xyz = None

        #: save atom-/residue-based values
        self.residues = PDBResidueProfiles( self )
        self.atoms = PDBProfiles( self )

        #: cached atom masks, calculated when first needed
        self.__maskCA = None
        self.__maskBB = None
        self.__maskHeavy = None
        #: cache position of chain breaks (clear when xyz changes)
        self.__chainBreaks = None

        #: starting positions of each residue
        self._resIndex = None
        #: starting positions of each chain
        self._chainIndex = None

        if noxyz:
            ## trick update to leave xyz untouched
            self.xyz = 0

        #: monitor changes of coordinates
        self.xyzChanged = 0

        self.forcePickle = 0

        #: version as of creation of this object
        self.initVersion = B.__version__

        #: to collect further informations
        self.info = { 'date':T.dateSortString() }

        if source != None:
            self.update( skipRes=skipRes, updateMissing=1, force=1,
                         headPatterns=headPatterns )

        if noxyz:
            ## discard coordinates, even when read from PDB file
            self.xyz = None


    def __getstate__(self):
        """
        called before pickling the object. (but also called by deepcopy)
        """
        self.slim()
        self.forcePickle = 0
        return self.__dict__

    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        ## until 2.0.1 PDBModel.atoms contained a list of dictionaries
        ## now we merged this info into PDBModel.aProfiles and renamed
        ## aProfiles back to .atoms
        if 'atoms' in state and type( state['atoms'] ) in [list, type(None)] :
            state['old_atoms'] = state['atoms']
            del state['atoms']

        self.__dict__ = state

        ## backwards compability
        self.__defaults() 

    def __len__(self):
        return self.lenAtoms()


    def __getitem__( self, k ):
        """
        Get atom profile or profile item or CrossView for one atom::
        
          m['prof1']         <==>  m.atoms.get( 'prof1' )         
          m['prof1','info1'] <==>  m.atoms.get( 'prof1','info1' )
          m[10]              <==>  CrossView( m.atoms, 10 )

        :return: profile OR meta infos thereof OR CrossView dict
        :rtype: list OR array OR any OR CrossView
        """
        if type( k ) is str:
            if k in self.atoms:
                return self.atoms.get( k )
            if k in self.residues:
                return self.residues.get( k )
            if k in self.info:
                return self.info[ k ]

            return self.profile( k )

        if type( k ) is tuple:
            return self.profileInfo( k[0] )[ k[1] ]

        return self.atoms[k]


    def __setitem__( self, k, v ):
        """
        Set atom profile or profile item (or meta info)::
        
          m['prof1'] = range(10)  
          is same as        <==> m.atoms.set( 'prof1', range(10) )
          OR                <==> m.residues.set( 'prof1', range(10) )

          m['prof1','info1'] = 'comment'
          is same as        <==> m.atoms.setInfo('prof1',info1='comment')
          OR                <==> m.residues.setInfo('prof1',info1='comment')

          m['version'] = '1.0.0'    <==> m.info['version'] = '1.0.0'
          (but only if 'version' already exists in `m.info`) 

        :return: item
        :rtype: any        
        """
        if type( k ) is str:
            if k in self.atoms:
                return self.atoms.set( k, v )
            if k in self.residues:
                return self.residues.set( k, v )
            if k in self.info:
                self.info[ k ] = v
                return
            if v is not None and len( v ) == self.lenAtoms():
                return self.atoms.set( k, v )
            if v is not None and len( v ) == self.lenResidues():
                return self.residues.set( k, v )
            raise ProfileError('Value cannot clearly be assigned to either atom or '+\
                  'residue profiles')

        if type( k ) is tuple:
            key, infokey = k
            if key in self.residues:
                self.residues[key, infokey] = v
                return
            self.atoms[key, infokey] = v
            return

        raise ProfileError('Cannot interpret %r as profile name or profile info record' % k)


    def __getslice__( self, *arg ):
        """
        Get list of CrossViews::
          m[0:100:5] <==> [ CrossView(m.atoms, i) for i in range(0,100,5) ]
        """
        return self.atoms.__getslice__( *arg )

    def __iter__( self ):
        return self.atoms.iterCrossViews()


    def __repr__( self ):
        code = self.pdbCode or ''
        return  '[%s %s %5i atoms, %4i residues, %2i chains]' % \
                ( self.__class__.__name__, code, self.lenAtoms(lookup=False), 
                  self.lenResidues(), self.lenChains() )

    def __str__( self ):
        return self.__repr__()

    def report( self, prnt=True, plot=False, clipseq=60 ):
        """
        Print (or return) a brief description of this model.

        :param prnt: directly print report to STDOUT (default True)
        :type  prnt: bool
        :param plot: show simple 2-D line plot using gnuplot [False]
        :type  plot: bool
        :param clipseq: clip chain sequences at this number of letters [60]
        :type  clipseq: int
        :return: if prnt==True: None, else: formatted description of this model
        :rtype: None or str
        """
        r = self.__repr__()

        for c in range( self.lenChains() ):
            m = self.takeChains( [c] )
            r += '\n\t* chain %-2i(%s): %s' % ( c, m['chain_id'][0],
                                          T.clipStr( m.sequence(), clipseq ) )

        r += '\n source: ' + repr( self.source )
        r += '\n %2i atom profiles:    %s' % ( len( self.atoms ), 
                                               T.clipStr( repr(list(self.atoms.keys())), 57 )) 
        r += '\n %2i residue profiles: %s' % ( len( self.residues ),
                                               T.clipStr( repr(list(self.residues.keys())), 57 ))
        r += '\n %2i info records:     %s' % ( len( self.info ), 
                                               T.clipStr( repr( list(self.info.keys()) ), 57 ))

        if plot:
            self.plot()

        if prnt:
            print(r)
        else:
            return r

    def plot( self, hetatm=False ):
        """
        Get a quick & dirty overview over the content of a PDBModel. plot
        simply creates a 2-D plot of all x-coordinates versus all y coordinates, 
        colored by chain. This is obviously not publication-quality ;-). 
        Use the Biskit.Pymoler class for real visalization.

        :param hetatm: include hetero & solvent atoms (default False)
        :type  hetatm: bool
        """
        from biskit import gnuplot

        m = self
        if not hetatm:
            mask = self.maskHetatm()
            mask = mask + m.maskSolvent()
            m = self.compress( N0.logical_not( mask ) )

        chains = [ self.takeChains( [i] ) for i in range( m.lenChains())]
        xy = [ list(zip( m.xyz[:,0], m.xyz[:,1] )) for m in chains ]

        gnuplot.plot( *xy )


    def __vintageCompatibility( self ):
        """
        backward compatibility to vintage PDBModels < 2.0.0
        """
        ## convert first generation profile dictionaries into new ProfileCollections
        if 'resProfiles' in self.__dict__:
            self.residues=PDBProfiles( self,
                                       profiles=getattr(self,'resProfiles',{}),
                                       infos=getattr(self,'resProfiles_info',{}) )
            try:
                del self.resProfiles; del self.resProfiles_info
            except: pass

        if 'atomProfiles' in self.__dict__:
            self.atoms=PDBProfiles( self,
                                    profiles=getattr(self,'atomProfiles',{}),
                                    infos=getattr(self,'atomProfiles_info',{}) )
            try:
                del self.atomProfiles; del self.atomProfiles_info
            except: pass

        ## fix old bug: slim() was creating a self.xyx instead of xyz
        if getattr( self, 'xyx', 0 ):
            del self.xyx

        ## first generation source was just a simple string
        if type( self.source ) == str:
            self.source = LocalPath( self.source )

        self.__validSource = getattr( self, '_PDBModel__validSource', 0)

        self.initVersion = getattr( self, 'initVersion', 'old PDBModel')

        self.forcePickle = getattr( self, 'forcePickle', 0 )

        self.info = getattr( self, 'info', { 'date':T.dateSortString() } )

        ## fix previous bug; old PDBModel pickles often have stray terAtoms
        ## records
        if getattr( self, '_PDBParseFile__terAtoms', None) is not None:
            del self._PDBParseFile__terAtoms



    def __defaults(self ):
        """
        backwards compatibility to earlier pickled models
        """
        self.__vintageCompatibility()

        ## if there were not even old profiles...
        if getattr( self, 'atoms', 0) == 0:
            self.atoms = PDBProfiles(self)
        if getattr( self, 'residues', 0) == 0:
            self.residues = PDBResidueProfiles(self)

        ## between release 2.0.1 and 2.1, aProfiles were renamed to atoms
        if getattr( self, 'aProfiles', None) is not None:
            self.atoms = self.aProfiles
            del self.aProfiles
        ## between release 2.0.1 and 2.1, rProfiles were renamed to residues
        if getattr( self, 'rProfiles', None) is not None:
            self.residues = self.rProfiles
            del self.rProfiles

        ## old aProfiles and rProfiles didn't keep a reference to the parent
        if getattr( self.atoms, 'model', 0) == 0:
            self.atoms.model = self
            self.residues.model = self

        ## biskit <= 2.0.1 kept PDB infos in list of dictionaries
        atoms = getattr( self, 'old_atoms', 0)
        if not atoms == 0:
            ## atoms to be fetched from external source
            if atoms is None:
                for k in self.PDB_KEYS:
                    self.atoms.set( k, None, changed=0 )

            else:
                atoms = B.DictList( atoms )
                for k in self.PDB_KEYS:
                    self.atoms.set( k, atoms.valuesOf( k ), 
                                    changed=getattr(self, 'atomsChanged',1) )

            del self.old_atoms
            del self.atomsChanged

        ## biskit <= 2.0.1 kept positions of TER records in a separate list
        ter_atoms = getattr( self, '_PDBModel__terAtoms', 0)
        if ter_atoms:
            mask = N0.zeros( self.atoms.profLength() )
            N0.put( mask, ter_atoms, 1 )
            self.atoms.set('after_ter', mask,
                           comment='rebuilt from old PDBModel.__terAtoms')
        if ter_atoms is not 0:
            del self.__terAtoms

        ## biskit <= 2.0.1 cached a volatile index in __resIndex & __chainIndex
        self._resIndex = getattr( self, '_resIndex', None)
        self._chainIndex=getattr( self, '_chainIndex', None)
        if getattr( self, '__resIndex', None) is not None:
            try:
                del self.__resIndex, self.__chainIndex
            except:
                print('DEBUG ', T.lastError())
                pass

        self.__maskCA = getattr( self, '__maskCA', None )
        self.__maskBB = getattr( self, '__maskBB', None )
        self.__maskHeavy = getattr( self, '__maskHeavy', None )

        ## test cases of biskit < 2.3 still contain Numeric arrays
        if self.xyz is not None and type( self.xyz ) is not N.ndarray:
            self.xyz = N0.array( self.xyz )
        if self._resIndex is not None and type( self._resIndex ) is not N.ndarray:
            self._resIndex = N0.array( self._resIndex )
        if self._chainIndex is not None and type(self._chainIndex) is not N.ndarray:
            self._chainIndex = N0.array( self._chainIndex )

        try:
            del self.caMask, self.bbMask, self.heavyMask
        except:
            pass


    def update( self, skipRes=None, updateMissing=0, force=0, 
                headPatterns=[] ):
        """
        Read coordinates, atoms, fileName, etc. from PDB or
        pickled PDBModel - but only if they are currently empty.
        The atomsChanged and xyzChanged flags are not changed.

        :param skipRes: names of residues to skip if updating from PDB
        :type  skipRes: list of str
        :param updateMissing: 0(default): update only existing profiles
        :type  updateMissing: 0|1
        :param force: ignore invalid source (0) or report error (1)
        :type  force: 0|1
        :param headPatterns: [(putIntoKey, regex)] extract given REMARKS
        :type  headPatterns: [(str, str)]

        :raise PDBError: if file can't be unpickled or read: 
        """
        source = self.validSource()

        if source is None and force:
            raise PDBError( str(self.source) + ' is not a valid source.')

        if source is None:
            return

        parser = PDBParserFactory.getParser( source )
        parser.update(self, source, skipRes=skipRes,
                      updateMissing=updateMissing, force=force,
                      headPatterns=headPatterns )


    def setXyz(self, xyz ):
        """
        Replace coordinates.

        :param xyz: Numpy array ( 3 x N_atoms ) of float
        :type  xyz: array

        :return: array( 3 x N_atoms ) or None, old coordinates
        :rtype: array
        """
        old = self.xyz
        self.xyz = xyz

        self.xyzChanged = self.xyzChanged or \
            not mathUtils.arrayEqual(self.xyz,old )
        return old


    def setSource( self, source ):
        """
        :param source: LocalPath OR PDBModel OR str
        """
        if type( source ) == str and len( source ) != 4:
            self.source = LocalPath( source )
        else:
            self.source = source
        self.__validSource = 0


    def getXyz( self, mask=None ):
        """
        Get coordinates, fetch from source PDB or pickled PDBModel,
        if necessary.

        :param mask: atom mask
        :type  mask: list of int OR array of 1||0

        :return: xyz-coordinates, array( 3 x N_atoms, Float32 )
        :rtype: array 
        """
        if self.xyz is None and self.validSource() is not None:
            self.update( force=1 )

        if self.xyz is None:
            ## empty array that can be concatenated to other xyz arrays
            return N0.zeros( (0,3), N0.Float32 )

        if mask is None:
            return self.xyz

        return N0.compress( mask, self.xyz, 0 )


    def getAtoms( self, mask=None ):
        """
        Get atom CrossViews that can be used like dictionaries.
        Note that the direct manipulation of individual profiles is more 
        efficient than the manipulation of CrossViews (on profiles)!

        :param mask: atom mask
        :type  mask: list of int OR array of 1||0

        :return: list of CrossView dictionaries
        :rtype: [ :class:`ProfileCollection.CrossView` ]
        """
        r = self.atoms.toCrossViews()

        if mask is None:
            return r

        return [ r[i] for i in N0.nonzero( mask ) ]


    def profile( self, name, default=None, update=True, updateMissing=False ):
        """
        Use::
           profile( name, updateMissing=0) -> atom or residue profile

        :param name: name to access profile
        :type  name: str        
        :param default: default result if no profile is found, if None,
        try to update from source and raise error [None]
        :type  default: any
        :param update: update from source before returning empty profile [True]
        :type  update: bool
        :param updateMissing: update from source before reporting missing
                              profile [False]
        :type  updateMissing: 0||1

        :raise ProfileError: if neither atom- nor rProfiles contains |name|
        """
        if updateMissing and not name in self.atoms and \
           not name in self.residues:
            self.update( updateMissing=True )

        if name in self.atoms:
            return self.atoms.get( name, default,
                                   update=update, updateMissing=0)

        if name in self.residues:
            return self.residues.get( name, default,
                                      update=update, updateMissing=0)

        if default is not None:
            return default

        raise ProfileError( 'No profile info found with name '+str(name))


    def profileInfo( self, name, updateMissing=0 ):
        """
        Use:
        
           profileInfo( name ) -> dict with infos about profile

        :param name: name to access profile
        :type  name: str       
        :param updateMissing: update from source before reporting missing \
               profile. Guaranteed infos are:
                    
                    - 'version' (str)
                    - 'comment' (str)
                    - 'changed' (1||0)
                    
        :type  updateMissing: 0|1

        :raise ProfileError: if neither atom - nor rProfiles contains |name|
        """
        if updateMissing and not name in self.atoms and \
           not name in self.residues:
            self.update()

        if name in self.atoms:
            return self.atoms.getInfo( name )

        if name in self.residues:
            return self.residues.getInfo( name )

        raise ProfileError( 'No profile info found with name '+str(name))


    def removeProfile( self, *names ):
        """
        Remove residue or atom profile(s)

        Use:
        
           removeProfile( str_name [,name2, name3] ) -> 1|0,

        :param names: name or list of residue or atom profiles
        :type  names: str OR list of str

        :return: 1 if at least 1 profile has been deleted,
                 0 if none has been found
        :rtype: int
        """
        r = 0

        for n in names:
            if n in self.atoms:
                del self.atoms[n]
                r = 1

            if n in self.residues:
                del self.residues[n]
                r = 1
        return r


    def xyzIsChanged(self):
        """
        Tell if xyz or atoms have been changed compared to source file or
        source object (which can be still in memory).

        :return: xyz field has been changed with respect to source
        :rtype: (1||0, 1||0)
        """
        return self.xyzChanged


    def xyzChangedFromDisc(self):
        """
        Tell whether xyz can currently be reconstructed from a
        source on disc. Same as xyzChanged() unless source is another not yet
        saved PDBModel instance that made changes relative to its own source.

        :return: xyz has been changed
        :rtype: bool
        """
        if self.validSource() is None:
            return True

        if isinstance( self.source, B.PDBModel ):
            return self.xyzIsChanged() or \
                   self.source.xyzChangedFromDisc()

        return self.xyzIsChanged()


    def profileChangedFromDisc(self, pname):
        """
        Check if profile has changed compared to source.

        :return: 1, if profile |pname| can currently not be
                 reconstructed from a source on disc.
        :rtype: int

        :raise ProfileError: if there is no atom or res profile with pname
        """
        if self.validSource() is None:
            return True

        if isinstance( self.source, B.PDBModel ):
            return self.profileInfo( pname )['changed'] or \
                   self.source.profileChangedFromDisc( pname )

        return self.profileInfo( pname )['changed']


    def __slimProfiles(self):
        """
        Remove profiles, that haven't been changed from a direct
        or indirect source on disc
        **AUTOMATICALLY CALLED BEFORE PICKLING and by deepcopy**
        """
        for key in self.residues:

            if not self.profileChangedFromDisc( key ):
                self.residues.set( key, None )

        for key in self.atoms:

            if not self.profileChangedFromDisc( key ):
                self.atoms.set( key,  None )


    def slim( self ):
        """
        Remove xyz array and profiles if they haven't been changed and
        could hence be loaded from the source file (only if there is a source
        file...).
        **AUTOMATICALLY CALLED BEFORE PICKLING**
        **Currently also called by deepcopy via getstate**
        """
        ## remove atoms/coordinates if they are unchanged from an existing
        ## source
        ## override this behaviour with forcePickle
        if not self.forcePickle:

            if not self.xyzChangedFromDisc():
                self.xyz = None

            if type( self.xyz ) is N0.arraytype and self.xyz.dtype.char != 'f':
                self.xyz = self.xyz.astype(N0.Float32)

            self.__slimProfiles()

        self.__maskCA = self.__maskBB = self.__maskHeavy = None
        self.__validSource = 0


    def validSource(self):
        """
        Check for a valid source on disk.

        :return:  str or PDBModel, None if this model has no valid source
        :rtype: str or PDBModel or None
        """
        if self.__validSource == 0:

            if isinstance( self.source, LocalPath ) and self.source.exists():
                self.__validSource = self.source.local()
            else:
                if isinstance( self.source, B.PDBModel ):
                    self.__validSource = self.source
                else:
                    ## risky: the PDB code may not exist!
                    if type( self.source ) is str and len( self.source )==4:
                        self.__validSource = self.source
                    else:
                        self.__validSource = None

        return self.__validSource


    def sourceFile( self ):
        """
        Name of pickled source or PDB file. If this model has another
        PDBModel as source, the request is passed on to this one.

        :return: file name of pickled source or PDB file
        :rtype: str

        :raise PDBError: if there is no valid source
        """
        s = self.validSource()

        if s is None:
            raise PDBError('no valid source')

        if type( s ) == str:
            return s

        return self.source.sourceFile()


    def disconnect( self ):
        """
        Disconnect this model from its source (if any).

        .. note:: 
        
            If this model has an (in-memory) PDBModel instance as source,
            the entries of 'atoms' could still reference the same dictionaries.
        """
        self.update()

        try:
            self.fileName = self.fileName or self.sourceFile()
        except:
            pass

        self.setSource( None )

        self.xyzChanged = 1

        for p in self.residues:
            self.residues.setInfo( p, changed=1 )
        for p in self.atoms:
            self.atoms.setInfo( p, changed=1 )


    def getPdbCode(self):
        """
        Return pdb code of model.

        :return: pdb code
        :rtype: str
        """
        return self.pdbCode

    def setPdbCode(self, code ):
        """
        Set model pdb code.

        :param code: new pdb code
        :type  code: str
        """
        self.pdbCode = code


    def sequence(self, mask=None, xtable=molUtils.xxDic ):
        """
        Amino acid sequence in one letter code.

        :param mask: atom mask, to apply before  (default None)
        :type  mask: list or array
        :param xtable: dict {str:str}, additional residue:single_letter mapping
                       for non-standard residues (default molUtils.xxDic)
                       [currently not used]
        :type  xtable: dict

        :return: 1-letter-code AA sequence (based on first atom of each res).
        :rtype: str
        """
        firstAtm = self.resIndex()
        if mask is not None:
            m_first = N0.zeros( self.lenAtoms() )
            N0.put( m_first, firstAtm, 1 )
            m_first = mask * m_first
            firstAtm = N0.nonzero( m_first )

        l = self.atoms['residue_name']
        l = [ l[i] for i in firstAtm ]

        return ''.join( molUtils.singleAA( l, xtable ) )


    def xplor2amber( self, aatm=True, parm10=False ):
        """
        Rename atoms so that tleap from Amber can read the PDB.
        If HIS residues contain atoms named HE2 or/and HD2, the residue
        name is changed to HIE or HID or HIP, respectively. Disulfide bonds
        are not yet identified - CYS -> CYX renaming must be done manually
        (see AmberParmBuilder for an example). 
        Internally amber uses H atom names ala HD21 while (old) standard pdb 
        files use 1HD2. By default, ambpdb produces 'standard' pdb atom names 
        but it can output the less ambiguous amber names with switch -aatm.

        :param change: change this model's atoms directly (default:1)
        :type  change: 1|0
        :param aatm: use, for example, HG23 instead of 3HG2 (default:1)
        :type  aatm: 1|0
        :param parm10: adapt nucleic acid atom names to 2010 Amber forcefield
        :type  parm10: 1|0

        :return: [ {..} ], list of atom dictionaries
        :rtype: list of atom dictionaries
        """
        numbers = list(map( str, list(range(10)) ))

        ## nucleic acid atom names have changed in parm10;
        if parm10:  ## this evidently is a bit of a hack. Should be revisited.
            def __parm10rename( a ):
                if "'1" in a: return a.replace( "'1", "'" )
                if "'2" in a: return a.replace( "'2", "''" )
                if a == 'O1P': return 'OP1'
                if a == 'O2P': return 'OP2'
                if a == 'H5T': return "HO5'"
                if a == 'H3T': return "HO3'"
                return a
            self.atoms['name'] = list(map( __parm10rename, self.atoms['name'] ))

        resI = self.resIndex().tolist()
        resI = N0.concatenate( (resI, [len(self)] ) )

        names    = self.atoms['name']
        resnames = self.atoms['residue_name']

        for i in range( len(resI)-1 ):
            first = resI[i]
            last  = resI[i+1]

            res = resnames[first]

            for j in range(first, last):

                if aatm:
                    a = names[j]

                    if len(a)>2 and a[0] in numbers:
                        names[j] = a[1:] + a[0]

            if res == 'HIS':
                anames = names[first:last]

                if 'HE2' in anames: resnames[first:last]= ['HIE'] *(last-first)
                if 'HD1' in anames: resnames[first:last]= ['HID'] *(last-first)
                if 'HE2' in anames and 'HD1' in anames:
                    resnames[first:last] = ['HIP'] *(last-first)
            

    def renameAmberRes( self ):
        """
        Rename special residue names from Amber back into standard names
        (i.e CYX S{->} CYS )
        """
        l = self.atoms['residue_name']

        for i in range(len(l)):
            if l[i] == 'CYX':
                l[i] = 'CYS'
            if l[i] in  ['HIE','HID','HIP']:
                l[i] = 'HIS'

    def writePdb( self, fname, ter=1, amber=0, original=0, left=0, wrap=0,
                  headlines=None, taillines=None):
        """
        Save model as PDB file.

        :param fname: name of new file
        :type  fname: str
        :param ter: Option of how to treat the terminal record:
        
                    * 0 - don't write any TER statements
                    * 1 - restore original TER statements (doesn't work, \
                          if preceeding atom has been deleted) [default]
                    * 2 - put TER between all detected chains
                    * 3 - as 2 but also detect and split discontinuous chains
        
        :type  ter: int
        :param amber: amber formatted atom names
                      (implies ter=3, left=1, wrap=0) (default 0)
        :type  amber: 1||0
        :param original: revert atom names to the ones parsed in from PDB
                         (default 0)
        :type  original: 1||0
        :param left: left-align atom names (as in amber pdbs)(default 0)
        :type  left: 1||0
        :param wrap: write e.g. 'NH12' as '2NH1' (default 0)
        :type  wrap: 1||0
        :param headlines: [( str, dict or str)], list of record / data tuples::
                          e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
        :type  headlines: list of tuples 
        :param taillines: same as headlines but appended at the end of file
        :type  taillines: list of tuples 
        """
        try:
            f = IO.PDBFile( fname, mode='w' )

            numbers = list(map( str, list(range(10)) ))

            if amber:
                __resnames = copy.copy(self.atoms['residue_name'])
                __anames   = copy.copy(self.atoms['name'])
                self.xplor2amber()
                ter = ter if ter!=1 else 3 ## adjust unless changed from default
                wrap = 0
                left = 1

            if ter == 2 or ter == 3:
                ## tolist to workaround numeric bug
                terIndex = N0.array( self.chainIndex( breaks=(ter==3) ).tolist()[1:] )
            if ter == 1:
                terIndex = N0.nonzero( self.atoms['after_ter'] )

            if headlines:
                for l in headlines:
                    f.writeLine( l[0], l[1] )

            if taillines:
                for l in taillines:
                    f.writeLine( l[0], l[1] )

            i = -1
            for a in self.atoms.toDicts():
                i += 1

                ## fetch coordinates Vector
                a['position'] = self.xyz[ i ]

                aname = a['name']

                if original and not amber:
                    aname = a['name_original']

                ## PDBFile prints atom names 1 column too far left
                if wrap and len(aname) == 4 and aname[0] in numbers:
                    aname = aname[1:] + aname[0]
                if not left and len(aname) < 4:
                    aname = ' ' + aname.strip()

                a['name'] = aname

                ## write line
                f.writeLine( a['type'], a )

                ## write TER line with details from previous atom
                if (ter>0 and  i+1 in terIndex):
                    f.writeLine('TER', a )

            f.close()

            if amber:
                self.atoms['residue_name'] = __resnames
                self.atoms['name'] = __anames

        except:
            EHandler.error( "Error writing "+fname )


    def saveAs(self, path):
        """
        Pickle this PDBModel to a file, set the 'source' field to
        this file name and mark atoms, xyz, and profiles  as unchanged.
        Normal pickling of the object will only dump those data that can not
        be reconstructed from the source of this model (if any).
        saveAs creates a 'new source' without further dependencies.

        :param path: target file name
        :type  path: str OR LocalPath instance
        """
        try:
            self.update()

            ## pickle all atoms, coordinates, profiles, even if unchanged
            self.forcePickle = 1

            self.setSource( path )

            self.xyzChanged = 0

            for p in self.residues:
                self.residues.setInfo( p, changed=0 )
            for p in self.atoms:
                self.atoms.setInfo( p, changed=0 )

            T.dump( self, str(path) )

        except IOError as err:
            raise PDBError("Can't open %s for writing." % T.absfile(str(path)))


    def maskF(self, atomFunction, numpy=1 ):
        """
        Create list whith result of atomFunction( atom ) for each
        atom. (Depending on the return value of atomFunction, the
        result is not necessarily a mask of 0 and 1. Creating masks
        should be just the most common usage).

        Note:

        This method is slow compared to maskFrom because the dictionaries
        that are given to the atomFunction have to be created from aProfiles
        on the fly. If performance matters, better combine the result from
        several maskFrom calls, e.g. instead of::

          r = m.maskF( lambda a: a['name']=='CA' and a['residue_name']=='ALA' )

        use::

          r = m.maskFrom( 'name', 'CA' ) * m.maskFrom('residue_name', 'ALA')

        :param atomFunction: function( dict_from_aProfiles.toDict() ),
                             true || false (Condition)
        :type  atomFunction: 1||0
        :param numpy: 1(default)||0, convert result to Numpy array of int
        :type  numpy: int

        :return: Numpy array( [0,1,1,0,0,0,1,0,..], Int) or list
        :rtype: array or list
        """
        try:
            result = list(map( atomFunction, self.atoms.toDicts() ))
        except:

            ## fall-back solution: assign 0 to all entries that raise
            ## exception
            EHandler.warning("mask(): Error while mapping funtion "+
                             "to all atoms.")
            result = []

            for a in self.atoms.iterDicts():
                try:
                    result.append( atomFunction( a ) )
                ## put 0 if something goes wrong
                except :
                    EHandler.warning("mask(): Error while save-mapping ")
                    result.append(0)

        if numpy:
            return N0.array( result )
        return result


    def maskFrom( self, key, cond ):
        """
        Create an atom mask from the values of a specific profile.
        Example, the following three statements are equivalent:

          >>> mask = m.maskFrom( 'name', 'CA' )
          >>> mask = m.maskFrom( 'name', lambda a: a == 'CA' )
          >>> mask = N0.array( [ a == 'CA' for a in m.atoms['name'] ] )

        However, the same can be also achieved with standard numpy operators:
        
          >>> mask = numpy.array(m.atoms['name']) == 'CA'

        :param key: the name of the profile to use
        :type  key: str
        :param cond: either a function accepting a single value or a value or
                     an iterable of values (to allow several alternatives)
        :type  cond: function OR any OR [ any ]
        :return: array or list of indices where condition is met
        :rtype: list or array of int
        """

        if type( cond ) is types.FunctionType:
            return N0.array( list(map( cond, self.atoms[ key ] )) )

        ## several allowed values given
        elif type( cond ) in [ list, tuple ]:
            return N0.array( [ x in cond for x in self.atoms[key] ] )

        ## one allowed value given
        ## Numeric splits lists of str into 2-D char arrays, 'O' prevents that
        else:
            return N0.array( self.atoms[key] ) == cond


    def maskCA( self, force=0 ):
        """
        Short cut for mask of all CA atoms.

        :param force: force calculation even if cached mask is available
        :type  force: 0||1

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        if self.__maskCA is None or force:
            self.__maskCA = self.maskFrom( 'name', 'CA' )

        return self.__maskCA


    def maskBB( self, force=0, solvent=0 ):
        """
        Short cut for mask of all backbone atoms. Supports standard protein
        and DNA atom names. Any residues classified as solvent (water, ions)
        are filtered out.

        :param force: force calculation even if cached mask is available
        :type  force: 0||1
        :param solvent: include solvent residues (default: false)
        :type  solvent: 1||0

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        if self.__maskBB is None or force or solvent:
            mask = self.maskFrom( 'name', ['CA', 'C', 'N', 'O', 'H','OXT',
                                           "P","O5'","C5'","C4'","C3'","O3'"])
            if not solvent:
                mask = N0.logical_not(self.maskSolvent()) * mask

                self.__maskBB = mask  ## cache
            else:
                return mask  ## don't cache
            
        return self.__maskBB


    def maskHeavy( self, force=0 ):
        """
        Short cut for mask of all heavy atoms. ('element' <> H)

        :param force: force calculation even if cached mask is available
        :type  force: 0||1

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        if self.__maskHeavy is None or force:
            self.__maskHeavy = self.maskFrom( 'element', lambda a: a != 'H' )

        return self.__maskHeavy

    def maskH( self ):
        """
        Short cut for mask of hydrogens. ('element' == H)

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return N0.logical_not( self.maskHeavy() )


    def _maskCB( self ):
        """
        Short cut for mask of all CB I{and} CA of GLY.

        :return: mask of all CB and CA of GLY
        :rtype: array
        """
        f = lambda a: a['name'] == 'CB' or\
          a['residue_name'] == 'GLY' and a['name'] == 'CA' 

        return self.maskF( f )


    def maskCB( self ):
        """
        Short cut for mask of all CB I{and} CA of GLY.

        :return: mask of all CB plus CA of GLY
        :rtype: array
        """
        m_cb = self.maskFrom( 'name', 'CB' )
        m_g  = self.maskFrom( 'residue_name', 'GLY' )
        m_ca = self.maskCA()

        return m_cb + (m_g * m_ca)


    def maskH2O( self ):
        """
        Short cut for mask of all atoms in residues named TIP3, HOH and  WAT

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom( 'residue_name', ['TIP3','HOH','WAT'] )

    def maskSolvent( self ):
        """
        Short cut for mask of all atoms in residues named
        TIP3, HOH, WAT, Na+, Cl-, CA, ZN

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom('residue_name', ['TIP3','HOH','WAT','Na+', 'Cl-',
                                              'CA', 'ZN'])

    def maskHetatm( self ):
        """
        Short cut for mask of all HETATM 

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom( 'type', 'HETATM' )

    def maskProtein( self, standard=0 ):
        """
        Short cut for mask containing all atoms of amino acids.

        :param standard: only standard residue names (not CYX, NME,..)
                         (default 0)
        :type  standard: 0|1

        :return: array( 1 x N_atoms ) of 0||1,
                 mask of all protein atoms (based on residue name)
        :rtype: array
        """
        d = molUtils.aaDic
        if standard:
            d = molUtils.aaDicStandard

        names = list(map( str.upper, list(d.keys()) ))
        return N0.array(
            [ n.upper() in names for n in self.atoms['residue_name'] ] )


    def maskDNA( self ):
        """
        Short cut for mask of all atoms in DNA (based on residue name).

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom( 'residue_name', ['DA','DC','DG','DT'] )

    def maskRNA( self ):
        """
        Short cut for mask of all atoms in RNA (based on residue name).

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom( 'residue_name', ['A','C','G','U'] )

    def maskNA( self ):
        """
        Short cut for mask of all atoms in DNA or RNA
        (based on residue name).

        :return: array( 1 x N_atoms ) of 0||1
        :rtype: array
        """
        return self.maskFrom( 'residue_name',
                              ['A','C','G','U','T','DA','DC','DG','DT'] )

    def indicesFrom( self, key, cond ):
        """
        Get atom indices conforming condition applied to an atom profile.
        Corresponds to::
          >>> numpy.nonzero( m.maskFrom( key, cond) )

        :param key: the name of the profile to use
        :type  key: str
        :param cond: either a function accepting a single value or a value or
                     an iterable of values
        :type  cond: function OR any OR [any]
        :return: array of indices where condition is met
        :rtype : array of int
        """
        return N0.nonzero( self.maskFrom( key, cond) )


    def indices( self, what ):
        """
        Get atom indices conforming condition. This is a convenience method
        to 'normalize' different kind of selections (masks, atom names,
        indices, functions) to indices as they are e.g. required by
        :class:`PDBModel.take`. 

        :param what: Selection::
             - function applied to each atom entry,
                e.g. lambda a: a['residue_name']=='GLY'
             - list of str, allowed atom names
             - list of int, allowed atom indices OR mask with only 1 and 0
             - int, single allowed atom index
        :type  what: function OR list of str or int OR int

        :return: N_atoms x 1 (0||1 )
        :rtype: Numeric array

        :raise PDBError: if what is neither of above
        """
        ## lambda funcion
        if type( what ) is types.FunctionType:
            return N0.nonzero( self.maskF( what) )

        if type( what ) is list or type( what ) is N0.arraytype:

            ## atom names
            if type( what[0] ) == str:
                return self.indicesFrom( 'name', what )

            if isinstance( what[0] , int) or \
               (isinstance(what, N.ndarray) and what.dtype in [int, bool]):
                ## mask
                if len( what ) == self.lenAtoms() and max( what ) < 2:
                    return N0.nonzero( what )
                ## list of indices
                else:
                    return what

        ## single index
        if isinstance( what , (int, N.integer)):
            return N0.array( [what], N0.Int )

        raise PDBError("PDBModel.indices(): Could not interpret condition ")


    def mask( self, what ):
        """
        Get atom mask. This is a convenience method to 'normalize'
        different kind of selections (masks, atom names, indices,
        functions) to a mask as it is e.g. required by :class:`PDBModel.compress`.

        :param what: Selection::
                     - function applied to each atom entry,
                        e.g. lambda a: a['residue_name']=='GLY'
                     - list of str, allowed atom names
                     - list of int, allowed atom indices OR mask with
                       only 1 and 0
                     - int, single allowed atom index
        :type  what: function OR list of str or int OR int

        :return: N_atoms x 1 (0||1 )
        :rtype: Numeric array

        :raise PDBError: if what is neither of above
        """
        ## lambda funcion
        if type( what ) == types.FunctionType:
            return self.maskF( what )

        if type( what ) == list or type( what ) is N0.arraytype:

            ## atom names
            if type( what[0] ) == str:
                return self.maskFrom( 'name', what)

            if isinstance( what[0] , int) or \
               (isinstance(what, N.ndarray) and what.dtype in [int, bool]):
                ## mask
                if len( what ) == self.lenAtoms() and max( what ) < 2:
                    return what
                ## list of indices
                else:
                    r = N0.zeros( self.lenAtoms(),N0.Int )
                    N0.put( r, what, 1 )
                    return r

        ## single index
        if isinstance( what , int):
            return self.mask( [what] )

        raise PDBError("PDBModel.mask(): Could not interpret condition ")


    def index2map( self, index, len_i ):
        """
        Create a map of len_i length, giving the residue(/chain) numer of
        each atom, from list of residue(/chain) starting positions.

        :param index: list of starting positions, e.g. [0, 3, 8]
        :type  index: [ int ] or array of int
        :param len_i: length of target map, e.g. 10
        :type  len_i: int

        :return: list mapping atom positions to residue(/chain) number,
                 e.g. [0,0,0, 1,1,1,1,1, 2,2] from above example
        :rtype: array of int (and of len_i length)
        """
        index = N0.concatenate( (index, [len_i]) )
        delta = index[1:] - index[:-1] 
        ## Numeric: delta = N0.take( index, range(1, len(index) ) ) - index[:-1]
        return N0.repeat( list(range(len(delta))), delta.astype( N0.Int32) )


    def map2index( self, imap ):
        """
        Identify the starting positions of each residue(/chain) from a map
        giving the residue(/chain) number of each atom.

        :param imap: something like [0,0,0,1,1,1,1,1,2,2,2,...]
        :type  imap: [ int ]

        :return: list of starting positions, e.g. [0, 3, 8, ...] in above ex.
        :rtype: array of int
        """
        try:
            imap = N0.concatenate( (imap, [imap[-1]] ) )
            delta = imap[1:] - imap[:-1] 
            # Numeric: delta = N0.take( imap, range(1, len(imap) ) ) - imap[:-1]
            r = N0.nonzero( delta ) + 1
            return N0.concatenate( ( [0], r ) )

        except IndexError:
            ## handle empty imap parameter
            return N0.zeros(0)


    def extendMask( self, mask, index, len_i ):
        """
        Translate a mask that is defined,e.g., on residues(/chains) to a mask
        that is defined on atoms.

        :param mask : mask marking positions in the list of residues or chains
        :type  mask : [ bool ] or array of bool or of 1||0
        :param index: starting positions of all residues or chains
        :type  index: [ int ] or array of int
        :param len_i: length of target mask
        :type  len_i: int

        :return: mask that blows up the residue / chain mask to an atom mask
        :rtype: array of bool
        """
        index = N0.concatenate( (index, [len_i]) )
        delta = index[1:] - index[:-1] 

        return N0.repeat( mask, delta.astype( N0.Int32 ) )


    def extendIndex( self, i, index, len_i ):
        """
        Translate a list of positions that is defined, e.g., on residues
        (/chains) to a list of atom positions AND also return the starting
        position of each residue (/chain) in the new sub-list of atoms.

        :param i : positions in higher level list of residues or chains
        :type  i : [ int ] or array of int
        :param index: atomic starting positions of all residues or chains
        :type  index: [ int ] or array of int
        :param len_i: length of atom index (total number of atoms)
        :type  len_i: int

        :return: (ri, rindex) - atom positions & new index
        :rtype:  array of int, array of int
        """
        ## catch invalid indices
        i = self.__convert_negative_indices( i, len( index ) )
        if len(i)==0:
            return i, N0.array( [], int )
        if max( i ) >= len( index ) or min( i ) < 0:
            raise PDBError("invalid indices")

        ## last atom of each residue / chain
        stop = N0.concatenate( (index[1:], [len_i]) ) - 1

        ifrom = N0.take( index, i )
        ito   = N0.take( stop, i )

        ## number of atoms in each of the new residues 
        rangelen = ito - ifrom + 1

        rindex   = N0.concatenate( ([0], N0.cumsum( rangelen[:-1] )) )

        ## (1) repeat position of first atom in each residue as often as there
        ## are atoms in this residue. (2) add a range array so that numbers
        ## are increasing from each atom to the next but (3) reset the added
        ## range to 0 at each residue starting position (-delta).

        ri    = N0.repeat( ifrom,  rangelen )
        delta = N0.repeat( rindex, rangelen )

        ri =  ri + N0.arange( len(ri), dtype=N0.Int32 ) - delta

        return ri, rindex


    def atom2resMask( self, atomMask ):
        """
        Mask (set 0) residues for which all atoms are masked (0) in atomMask.

        :param atomMask: list/array of int, 1 x N_atoms
        :type  atomMask: list/array of int

        :return: 1 x N_residues (0||1 )
        :rtype: array of int
        """
        res_indices = self.atom2resIndices( N0.nonzero( atomMask) )
        r = N0.zeros( self.lenResidues() )
        N0.put( r, res_indices, 1 )
        return r


    def atom2resIndices( self, indices ):
        """
        Get list of indices of residues for which any atom is in indices.

        Note: in the current implementation, the resulting residues are
        returned in their old order, regardless of the order of input
        positions.

        :param indices: list of atom indices
        :type  indices: list of int

        :return: indices of residues
        :rtype: list of int
        """
        new_resmap = N0.take( self.resMap(), indices )
        resIndex   = self.map2index( new_resmap )

        return N0.take( new_resmap, resIndex )


    def res2atomMask( self, resMask ):
        """
        Convert residue mask to atom mask.

        :param resMask: list/array of int, 1 x N_residues
        :type  resMask: list/array of int

        :return: 1 x N_atoms
        :rtype: array of int
        """
        return self.extendMask( resMask, self.resIndex(), self.lenAtoms() )

    def __convert_negative_indices( self, indices, length ):
        """
        Replace negative indices by their positive equivalent.
        :return: modified copy of indices (or unchanged indices itself)
        :rtype: array of int
        """
        if len(indices)==0 or min( indices ) >= 0:
            return indices

        indices = copy.copy( N0.array( indices ) )

        negatives = N.flatnonzero( indices < 0 )

        a = N0.zeros( len( indices ) )
        N0.put( a, negatives, length )

        indices += a

        #for i in negatives:
            #indices[ i ] = length + indices[i] ## substract from end

        return indices


    def res2atomIndices( self, indices ):
        """
        Convert residue indices to atom indices.

        :param indices: list/array of residue indices
        :type  indices: list/array of int

        :return: array of atom positions
        :rtype: array of int
        """
        if max( indices ) > self.lenResidues() or min( indices ) < 0:
            raise PDBError("invalid residue indices")

        return self.extendIndex( indices, self.resIndex(), self.lenAtoms() )[0]


    def atom2chainIndices( self, indices, breaks=0 ):
        """
        Convert atom indices to chain indices. Each chain is only
        returned once.

        :param indices: list of atom indices
        :type  indices: list of int
        :param breaks: look for chain breaks in backbone coordinates (def. 0)
        :type  breaks: 0||1

        :return: chains any atom which is in indices
        :rtype: list of int
        """
        new_map = N0.take( self.chainMap( breaks=breaks ), indices )
        index   = self.map2index( new_map )

        return N0.take( new_map, index )


    def atom2chainMask( self, atomMask, breaks=0 ):
        """
        Mask (set to 0) chains for which all atoms are masked (0) in atomMask.
        Put another way: Mark all chains that contain any atom that is marked
        '1' in atomMask.

        :param atomMask: list/array of int, 1 x N_atoms
        :type  atomMask: list/array of int

        :return: 1 x N_residues (0||1 )
        :rtype: array of int
        """
        indices = self.atom2chainIndices( N0.nonzero( atomMask), breaks=breaks )
        r = N0.zeros( self.lenChains(breaks=breaks) )
        N0.put( r, indices, 1 )
        return r


    def chain2atomMask( self, chainMask, breaks=0 ):
        """
        Convert chain mask to atom mask.

        :param chainMask: list/array of int, 1 x N_chains
        :type  chainMask: list/array of int
        :param breaks: look for chain breaks in backbone coordinates (def. 0)
        :type  breaks: 0||1

        :return: 1 x N_atoms
        :rtype: array of int
        """
        return self.extendMask( chainMask, self.chainIndex( breaks=breaks ),
                                self.lenAtoms() )


    def chain2atomIndices( self, indices, breaks=0 ):
        """
        Convert chain indices into atom indices.

        :param indices: list/array of chain indices
        :type  indices: list/array of int

        :return: array of atom positions, new chain index
        :rtype: array of int
        """
        if max( N0.absolute(indices) ) > self.lenChains( breaks=breaks ):
            raise PDBError("invalid chain indices")

        return self.extendIndex( indices, self.chainIndex( breaks=breaks ),
                                 self.lenAtoms() )[0]


    def res2atomProfile( self, p ):
        """
        Get an atom profile where each atom has the value its residue has
        in the residue profile.

        :param p: name of existing residue profile OR ...
                  [ any ], list of lenResidues() length
        :type  p: str

        :return: [ any ] OR array, atom profile
        :rtype: list or array
        """
        if type( p ) is str:
            p = self.residues.get( p )

        isArray = isinstance( p, N0.arraytype )

        resMap = self.resMap()

        r = [ p[ resMap[a] ] for a in range( len(resMap) ) ]

        if isArray:
            r = N0.array( r )

        return r


    def atom2resProfile( self, p, f=None ):
        """
        Get a residue profile where each residue has the value that its first
        atom has in the atom profile.
        :param p: name of existing atom profile OR ...
                  [ any ], list of lenAtoms() length
        :type  p: str
        :param f: function to calculate single residue from many atom values 
                  f( [atom_value1, atom_value2,...] ) -> res_value
                  (default None, simply take value of first atom in each res.)
        :type  f: func

        :return: [ any ] OR array, residue profile
        :rtype: list or array
        """
        if type( p ) is str:
            p = self.atoms.get( p )

        isArray = isinstance( p, N0.arraytype )

        if not f:
            r = N0.take( p, self.resIndex() )
        else:
            r = [ f( values ) for values in self.profile2resList( p ) ]
            r = N0.array( r )

        if not isArray:
            return r.tolist()
        return r


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2mask( str_profname, [cutoff_min, cutoff_max=None])

        :param cutoff_min: low value cutoff (all values >= cutoff_min)
        :type  cutoff_min: float
        :param cutoff_max: high value cutoff (all values < cutoff_max)
        :type  cutoff_max: float

        :return: mask len( profile(profName) ) x 1||0
        :rtype: array

        :raise ProfileError: if no profile is found with name profName
        """
        if profName in self.atoms:
            return self.atoms.profile2mask( profName, cutoff_min, cutoff_max)
        return self.residues.profile2mask( profName, cutoff_min, cutoff_max)


    def profile2atomMask( self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2atomMask( str_profname, [cutoff_min, cutoff_max=None])
        Same as :class:`profile2mask`, but converts residue mask to atom mask.

        :param cutoff_min: low value cutoff
        :type  cutoff_min: float
        :param cutoff_max: high value cutoff
        :type  cutoff_max: float

        :return: mask N_atoms x 1|0
        :rtype: array

        :raise ProfileError: if no profile is found with name profName
        """
        r = self.profile2mask( profName, cutoff_min, cutoff_max )

        if len( r ) == self.lenResidues():
            r = self.res2atomMask( r )

        return r


    def profile2resList( self, p ):
        """
        Group the profile values of each residue's atoms into a separate list.
        :param p: name of existing atom profile OR ...
                  [ any ], list of lenAtoms() length

        :return: a list (one entry per residue) of lists (one entry per resatom)
        :rtype: [ [ any ] ]
        """
        if type( p ) is str:
            p = self.atoms.get( p )

        rI = self.resIndex()       # starting atom of each residue
        rE = self.resEndIndex()    # ending atom of each residue

        r = [ p[ rI[res] : rE[res]+1 ] for res in range( self.lenResidues() ) ]
        return r

    def mergeChains( self, c1, id='', segid='', rmOxt=True,
                     renumberAtoms=False, renumberResidues=True):
        """
        Merge two adjacent chains. This merely removes all internal markers
        for a chain boundary. Atom content or coordinates are not modified.

        PDBModel tracks chain boundaries in an internal _chainIndex. However,
        there are cases when this chainIndex needs to be re-built and new 
        chain boundaries are then infered from jumps in chain- or segment 
        labelling or residue numbering. mergeChains automatically
        re-assigns PDB chain- and segment IDs as well as residue numbering
        to prepare for this situation.

        :param c1   : first of the two chains to be merged
        :type  c1   : int
        :param id   : chain ID of the new chain (default: ID of first chain)
        :type  id   : str
        :param segid: ew chain's segid (default: SEGID of first chain)
        :type  segid: str
        :param renumberAtoms: rewrite PDB serial numbering of the adjacent
                              chain to be consequtive to the last atom of the
                              first chain (default: False)
        :type  renumberAtoms: bool
        :param renumberResidues: shift PDB residue numbering so that the first
                                 residue of the adjacent chain follows the
                                 previous residue. Other than for atom
                                 numbering, later jumps in residue numbering
                                 are preserved. (default: True)
        :type  renumberResidues: bool
        """
        c1 = self.__convert_negative_indices( [c1], self.lenChains() )[0]

        oldI = self.chainIndex()
        assert len(oldI) > c1 + 1, 'no adjacent chain to be merged'

        ## remove chain boundary from chainIndex
        self._chainIndex = N0.concatenate( (oldI[:c1+1], oldI[c1+2:] ) ) 

        ## starting and ending position of (old) second chain
        i_start= oldI[ c1 ]
        if len( oldI ) > c1+2:
            i_next = oldI[ c1+2 ]
        else:
            i_next = len( self )
        i_scar = oldI[ c1+1 ]     ## (old) starting position of second chain
        n_atoms= i_next - i_start

        ## remove trace of PDB TER statement (if any)
        self['after_ter'][ i_scar ] = 0

        ## harmonize chain ID
        id = id or self['chain_id'][i_scar-1]
        self['chain_id'][ i_start : i_next ] = [ id ] * n_atoms

        ## harmonize segID
        segid = segid or self['segment_id'][i_scar-1]
        self['segment_id'][ i_start : i_next ] = [ segid ] * n_atoms

        ## harmonize PDB residue numbering (by *shifting* current numbers)
        if renumberResidues:
            first = self['residue_number'][ i_scar-1 ] + 1
            delta = first - self['residue_number'][ i_scar ]

            profile = self['residue_number'][i_scar:i_next] + delta
            self['residue_number'][i_scar:i_next] = profile

        ## harmonize PDB atom numbering (by *rewriting* current numbers)
        ## usually though, atoms already have consequtive numbering through the PDB
        if renumberAtoms:
            n = i_next - i_scar
            first = self['serial_number'][ i_scar-1 ] + 1 

            self['serial_number'][i_scar:i_next] = N0.arange(first, first + n)

        ## remove OXT and OT2 if requested
        if rmOxt:
            ## overkill: we actually would only need to look into last residue
            anames = N0.array( self.atoms['name'][i_start:i_scar] )
            i_oxt = N.flatnonzero( N0.logical_or( anames=='OXT', anames=='OT2' ))
            if len( i_oxt ) > 0:
                self.remove( i_oxt )


    def mergeResidues( self, r1, name='', residue_number=None, 
                       chain_id='', segment_id='',
                       renumberAtoms=False ):
        """
        Merge two adjacent residues. Duplicate atoms are labelled with
        alternate codes 'A' (first occurrence) to 'B' or later.
        :param r1: first of the two residues to be merged
        :type  r1: int
        :param name: name of the new residue (default: name of first residue)
        :type  name: str
        """
        r1 = self.__convert_negative_indices( [r1], self.lenResidues() )[0]

        oldI = self.resIndex()
        assert len(oldI) > r1 + 1, 'no adjacent residue to be merged'

        ## remove residue boundary from residue Index
        self._resIndex = N0.concatenate( (oldI[:r1+1], oldI[r1+2:] ) ) 

        ## starting and ending position of new fused and (old) second residue
        i_start= oldI[ r1 ]
        if len( oldI ) > r1+2:
            i_next = oldI[ r1+2 ]
        else:
            i_next = len( self )
        i_scar = oldI[ r1+1 ]     ## (old) starting position of second residue
        n_atoms= i_next - i_start

        ## move PDB TER statement (if any) to end of fused residue
        if i_next < len( self ):
            self['after_ter'][ i_next ] = self['after_ter'][ i_scar ]
        self['after_ter'][ i_scar ] = 0

        ## harmonize residue name
        name = name or self['residue_name'][i_scar-1]
        self['residue_name'][ i_start : i_next ] = [ name ] * n_atoms

        ## harmonize chain ID
        id = chain_id or self['chain_id'][i_scar-1]
        self['chain_id'][ i_start : i_next ] = [ id ] * n_atoms

        ## harmonize segID
        segid = segment_id or self['segment_id'][i_scar-1]
        self['segment_id'][ i_start : i_next ] = [ segid ] * n_atoms

        ## harmonize PDB residue numbering
        residue_number = residue_number or self['residue_number'][i_scar-1]
        self['residue_number'][i_start:i_next] = [residue_number] * n_atoms

        ## harmonize PDB atom numbering (by *rewriting* current numbers)
        if renumberAtoms:
            n = i_next - i_scar
            first = self['serial_number'][ i_scar-1 ] + 1 

            self['serial_number'][i_scar:i_next] = N0.arange(first, first + n)

        ## shift chain boundary (if any) to end of fused residue
        ## unless it's the end of the model or there is already a boundary there
        if i_scar in self.chainIndex():
            i = N.flatnonzero( self._chainIndex == i_scar )[0]
            if (not i_next in self._chainIndex) and (i_next != len(self)):
                self._chainIndex[ i ] = i_next
            else:
                self._chainIndex = N0.concatenate( self._chainIndex[:i],
                                                  self._chainIndex[i+1:] )

        ## mark duplicate atoms in the 'alternate' field of the new residue
        r = Residue( self, r1 )
        r.labelDuplicateAtoms()


    def concat( self, *models, **kw ):
        """
        Concatenate atoms, coordinates and profiles. source and fileName
        are lost, so are profiles that are not available in all models.
        model0.concat( model1 [, model2, ..]) -> single PDBModel.

        :param models:   models to concatenate
        :type  models:   one or more PDBModel instances
        :param newRes:   treat beginning of second model as new residue (True)
        :type  newRes:   bool
        :param newChain: treat beginning of second model as new chain (True)
        :type  newChain: bool

        Note: info records of given models are lost.
        """
        newRes  = kw.get('newRes', True)
        newChain= kw.get('newChain', True)

        if len( models ) == 0:
            return self

        m = models[0]

        r = self.__class__()

        self.update()  ## trigger update if xyz or any profile is None

        r.setXyz( N0.concatenate( ( self.getXyz(), m.getXyz() )  ) )

        r.setPdbCode( self.pdbCode )

        r.atoms = self.atoms.concat( m.atoms, )  
        r.residues = self.residues.concat( m.residues, ) 

        r.residues.model = r
        r.atoms.model = r

        append_I = m.resIndex() + self.lenAtoms()
        r._resIndex  = N0.concatenate((self.resIndex(), append_I ))

        append_I = m.chainIndex() +self.lenAtoms()
        r._chainIndex =N0.concatenate((self.chainIndex( singleRes=1 ), append_I))

        ## remove traces of residue or chain breaks
        if not newChain:
            r.mergeChains( self.lenChains() - 1 )

        if not newRes:
            r.mergeResidues( self.lenResidues() -1 )

        r.info = copy.deepcopy( self.info )

## leads to bug 3611835 
##        try:
##            k = max(self.biounit.keys())+1
##            r.residues['biomol'][self.lenResidues():] += k
##            r.biounit = self.biounit.append(m.biounit)
##            r.biounit.model = r
##        except AttributeError:
##            pass

        return r.concat( *models[1:] )


##     def removeChainBreaks( self, chains, breaks=False ):
##         """
##         Remove chain boundaries *before* given chain indices.
##         Example:
##            removeChainBreaks( [1,3] ) --> removes the first and third chain
##              break but keeps the second, e.g. this joins first and second chain
##              but also second and third chain.
##         Coordinates are not modified. removeChainBreaks( [0] ) doesn't make
##         sense.
##         :param chains: [ int ], chain breaks
##         """
##         if 0 in chains:
##             raise PDBError, 'cannot remove chain break 0'

##         cindex = self.chainIndex( breaks=breaks )

##         ## simple removal of terminal OXT and TER label, make it more robust!
##         remove = []
##         for i in chains:
##             lastatom = cindex[i] - 1
##             if self[ lastatom ]['name'] in ['OXT', 'OT2']:
##                 remove += [ lastatom ]
##             self['after_ter'][lastatom+1] = 0 

##         self.remove( remove )

##         ## update chain index
##         cindex = self.chainIndex( breaks=breaks )

##         mask = N0.ones( len( cindex ) )
##         N0.put( mask, chains, 0 )

##         self._chainIndex = N0.compress( mask, cindex )


    def take( self, i, rindex=None, cindex=None,
              *initArgs, **initKw ):
        """
        Extract a PDBModel with a subset of atoms:
        
          take( atomIndices ) -> PDBModel
          
        All other PDBModel methods that extract portions of the model (e.g.
        compress, takeChains, takeResidues, keep, clone, remove) are ultimately
        using `take()` at their core.
          
        Note:
        take employs fast numpy vector mapping methods to re-calculate
        the residue and chain index of the result model. The methods generally
        work but there is one scenario were this mechanism can fail: If take
        is used to create repetitions of residues or chains directly next to
        each other, these residues or chains can get accidentally merged. For
        this reason, calling methods can optionally pre-calculate and provide
        a correct version of the new residue or chain index (which will then
        be used as is).

        :param i: atomIndices, positions to take in the order to take
        :type  i: list/array of int
        
        :param rindex: optional residue index for result model after extraction
        :type  rindex: array of int

        :param cindex: optional chain index for result model after extraction
        :type  cindex: array of int
        
        :param initArgs: any number of additional arguments for constructor of \
                         result model
        
        :param initKw: any additional keyword arguments for constructure of \
                        result model

        :return: new PDBModel or sub-class
        :rtype: PDBModel
        """
        r = self.__class__( *initArgs, **initKw )

        ## the easy part: extract coordinates and atoms
        r.xyz = N0.take( self.getXyz(), i )
        r.xyzChanged = self.xyzChanged or not N.array_equal(r.xyz,self.xyz)

        r.atoms = self.atoms.take( i, r )

        ## more tricky: rescue residue borders and extract residue profiles
        new_resmap   = N0.take( self.resMap(), i )
        if rindex is not None:
            r._resIndex = rindex
        else:
            ## this can fail if residues are repeated in the selection
            r._resIndex = self.map2index( new_resmap )

        i_res     = N0.take( new_resmap, r._resIndex )
        r.residues = self.residues.take( i_res, r )

        ## now the same with chain borders (and later profiles)
        if cindex is not None:
            r._chainIndex = cindex
        else:
            ## this can fail if chains are repeated next to each other in i
            new_chainmap   = N0.take( self.chainMap(), i )
            r._chainIndex = self.map2index( new_chainmap )

        ## copy non-sequential infos
        r.info = copy.deepcopy( self.info )
        r.pdbCode = self.pdbCode
        r.fileName = self.fileName
        r.source = self.source

        ## copy the biounit
##        try:
##            r.biounit = self.biounit.take(i)
##            r.biounit.model = r
##        except AttributeError:
##            pass

        return r



    def keep( self, i ):
        """
        Replace atoms,coordinates,profiles of this(!) model with sub-set.
        (in-place version of N0.take() )

        :param i: atom positions to be kept
        :type  i: list or array of int
        """
        if len(i)==self.lenAtoms() and max(i)<2:
            EHandler.warning('dont use PDBModel.keep() with mask.', trace=0) 

        r = self.take( i )

        self.xyz = r.xyz
        self.xyzChanged = r.xyzChanged

        self._resIndex  = r._resIndex
        self._chainIndex= r._chainIndex

        self.atoms = r.atoms
        self.residues = r.residues

        self.info = r.info

        self.__maskCA = self.__maskBB = self.__maskHeavy = None
        self.__chainBreaks = None


    def clone( self ):
        """
        Clone PDBModel.

        :return: PDBModel / subclass, copy of this model,
                 see comments to numpy.take()
        :rtype: PDBModel
        """
        return self.take( self.atomRange() )


    def compress( self, mask, *initArgs, **initKw ):
        """
        Compress PDBmodel using mask.
        
            compress( mask ) -> PDBModel

        :param mask: array( 1 x N_atoms of 1 or 0 ):
        
                     * 1 .. keep this atom
        
        :type  mask: array

        :return: compressed PDBModel using mask
        :rtype: PDBModel
        """
        return self.take( N0.nonzero( mask ), *initArgs, **initKw )


    def remove( self, what ):
        """        
        Convenience access to the 3 different remove methods.
        The mask used to remove atoms is returned. This mask can be used
        to apply the same change to another array of same dimension as
        the old(!) xyz and atoms.

        :param what: Decription of what to remove:
        
              - function( atom_dict ) -> 1 || 0    (1..remove) OR
              - list of int [4, 5, 6, 200, 201..], indices of atoms to remove
              - list of int [11111100001101011100..N_atoms], mask (1..remove)
              - int, remove atom with this index
        
        :type  what: list of int or int

        :return: array(1 x N_atoms_old) of 0||1, mask used to compress the
                 atoms and xyz arrays. 
        :rtype: array

        :raise PDBError: if what is neither of above
        """
        mask = N0.logical_not( self.mask( what ) )
        self.keep( N0.nonzero(mask) )
        return mask


    def takeResidues( self, i ):
        """
        Copy the given residues into a new model.
        
        :param i: residue indices
        :type  i: [ int ]
        :return: PDBModel with given residues in given order
        :rtype: PDBModel
        """
        ##i, index = self.res2atomIndices( i )
        i, index = self.extendIndex( i, self.resIndex(), self.lenAtoms() )
        return self.take( i, rindex=index )


    def takeChains( self, chains, breaks=0, force=0 ):
        """
        Get copy of this model with only the given chains.
        
        Note, there is one very special scenario where chain boundaries can get
        lost: If breaks=1 (chain positions are based on normal chain boundaries
        as well as structure-based chain break detection) AND one or more
        chains are extracted several times next to each other, for example
        chains=[0, 1, 1, 2], then the repeated chain will be merged. So in
        the given example, the new model would have chainLength()==3. This case
        is tested for and a PDBIndexError is raised. Override with force=1 and
        proceed at your own risk. Which, in this case, simply means you should
        re-calculate the chain index after takeChains(). Example::
        
           >>> repeat = model.takeChains( [0,0,0], breaks=1, force=1 )
           >>> repeat.chainIndex( force=1, cache=1 )
           
        This works because the new model will have back-jumps in residue 
        numbering.

        :param chains: list of chains, e.g. [0,2] for first and third
        :type  chains: list of int
        :param breaks: split chains at chain breaks (default 0)
        :type  breaks: 0|1
        :param maxDist: (if breaks=1) chain break threshold in Angstrom
        :type  maxDist: float
        :param force: override check for chain repeats (only for breaks==1)
        :type  force: bool

        :return: PDBModel consisting of the given chains in the given order
        :rtype: PDBModel
        """
        ## i, index = self.chain2atomIndices( chains, breaks=breaks )
        i, index = self.extendIndex( chains, self.chainIndex( breaks=breaks ),
                                     self.lenAtoms() )
        if not breaks or len(chains)==0:
            return self.take( i, cindex=index )

        ## test for repeats:
        if not force:
            chains = N0.array( chains, int )
            delta = N0.concatenate( (chains[1:], [chains[-1]+1]) ) - chains
            if not N.all( delta != 0 ):
                raise PDBIndexError('Chain boundaries cannot be preserved for repeats.' +\
                   "Use 'force=1' to override, then re-calculate chainIndex().")
            
        ## Give up on repeat treatement: 
        ## more important: the new model's chain index should NOT include breaks
        return self.take( i )


    def addChainFromSegid(self, verbose=1):
        """
        Takes the last letter of the segment ID and adds it as chain ID.
        """
        chain_ids   = self.atoms['chain_id']
        segment_ids = self.atoms['segment_id']

        for i in self.atomRange():

            try:
                chain_ids[i] = segment_ids[i][-1]
            except:
                if verbose:
                    EHandler.warning("addChainId(): Problem with atom "+str(self[i]))


    def addChainId( self, first_id=None, keep_old=0, breaks=0 ):
        """
        Assign consecutive chain identifiers A - Z to all atoms.

        :param first_id: str (A - Z), first letter instead of 'A'
        :type  first_id: str 
        :param keep_old: don't override existing chain IDs (default 0)
        :type  keep_old: 1|0
        :param breaks: consider chain break as start of new chain (default 0)
        :type  breaks: 1|0
        """
        ids = self.atoms['chain_id']

        old_chains = []
        if keep_old:
            old_chains = N0.take( ids, self.resIndex() )
            old_chains = mathUtils.nonredundant( old_chains )
            if '' in old_chains: old_chains.remove('')

        letters = string.ascii_uppercase
        if first_id:
            letters = letters[ letters.index( first_id ): ]
        letters = mathUtils.difference( letters, old_chains )

        chainMap = self.chainMap( breaks=breaks )

        try:
            for i in self.atomRange():

                if not (keep_old and ids[i] in old_chains):
                    ids[i] = letters[ chainMap[i] ]

        except IndexError:
            raise PDBError('Too many chains, running out of letters.')


    def renumberResidues( self, mask=None, start=1, addChainId=1 ):
        """
        Make all residue numbers consecutive and remove any insertion
        code letters. Note that a backward jump in residue numbering
        (among other things) is interpreted as end of chain by
        chainMap() and chainIndex() when a PDB file is loaded.

        :param mask: [ 0||1 x N_atoms ] atom mask to apply BEFORE
        :type  mask: list of int
        :param start: starting number (default 1)
        :type  start: int
        :param addChainId: add chain IDs if they are missing
        :type  addChainId: 1|0
        """
        if addChainId:
            self.addChainId( keep_old=1, breaks=0 )

        i = start
        for res in self.resList( mask ):

            for a in res:
                a['residue_number'] = i
                a['insertion_code'] = ''

            i += 1


    def atomRange( self ):
        """
        >>> m.atomRange() == range( m.lenAtoms() )
        
        :return: integer range for lenght of this model
        :rtype: [ int ]
        """
        return list(range( self.lenAtoms()))


    def lenAtoms( self, lookup=True ):
        """
        Number of atoms in model.

        :return: number of atoms
        :rtype: int
        """
        if not self.xyz is None:
            return len( self.xyz )

        if len( self.atoms ) > 0:
            r = self.atoms.profLength( default=-1 )
            if r != -1:
                return r

        if self.source is None and not lookup:
            return 0

        return len( self.getXyz() )


    def lenResidues( self ):
        """
        Number of residues in model.

        :return: total number of residues
        :rtype: int
        """
##         if self._resIndex is None:
##             return 0

        return len( self.resIndex() )


    def lenChains( self, breaks=0, maxDist=None, singleRes=0, solvent=0 ):
        """
        Number of chains in model.

        :param breaks: detect chain breaks from backbone atom distances (def 0)
        :type  breaks: 0||1
        :param maxDist: maximal distance between consequtive residues
                        [ None ] .. defaults to twice the average distance
        :type  maxDist: float
        :param singleRes: allow chains consisting of single residues (def 0)
        :type  singleRes: 1||0
        :param solvent: also check solvent residues for "chain breaks" (def 0)
        :type  solvent: 1||0

        :return: total number of chains
        :rtype: int
        """
        try:
            return len( self.chainIndex( breaks=breaks, maxDist=maxDist,
                                         singleRes=singleRes, solvent=solvent))
        except IndexError:  ## empty residue map
            return 0


    def resList( self, mask=None ):
        """
        Return list of lists of atom pseudo dictionaries per residue,
        which allows to iterate over residues and atoms of residues.

        :param mask: [ 0||1 x N_atoms ] atom mask to apply BEFORE
        :type  mask: 

        :return: a list (one per residue) of lists (one per atom) of dictionaries
            ::
            
                [ [ CrossView{'name':'N', ' residue_name':'LEU', ..},          
                    CrossView{'name':'CA', 'residue_name':'LEU', ..} ],   
                
                  [ CrossView{'name':'CA', 'residue_name':'GLY', ..}, .. ] 
                ]
            
        :rtype: [ [ `biskit.ProfileCollection.CrossView` ] ]    
        """
        ri = N0.concatenate( (self.resIndex( mask=mask ), [self.lenAtoms()] ) )
        resLen = len( ri ) - 1
        atoms = self.getAtoms()
        if mask is not None:
            atoms = N0.compress( mask, atoms ).tolist()

        return [ atoms[ ri[res] : ri[res+1] ] for res in range( resLen ) ] 


    def resModels( self, i=None ):
        """
        Creates one new PDBModel for each residue in the parent PDBModel.
        
        :param i: range of residue positions (default: all residues)
        :type  i: [ int ] or array( int )

        :return: list of PDBModels, one for each residue
        :rtype: [ `PDBModel` ]
        """
        ri = self.resIndex()
        re = self.resEndIndex()
        if i is None:
            i = N0.arange( len(ri) )
##        xyz = self.getXyz()

        result = [ self.take(N0.arange(ri[x],re[x]+1)) for x in i ]
        return result


    def resMapOriginal(self, mask=None):
        """
        Generate list to map from any atom to its ORIGINAL(!) PDB
        residue number.

        :param mask: [00111101011100111...] consider atom: yes or no
                     len(mask) == N_atoms
        :type  mask: list of int (1||0)

        :return: list all [000111111333344444..] with residue number
                 for each atom
        :rtype: list of int
        """
        ## by default take all atoms
        if mask is None: mask = N0.ones( self.lenAtoms() , N0.Int )

        ## apply mask to this list
        return N0.compress( mask, self.atoms['residue_number'] )


    def __inferResIndex( self ):
        """
        Determine residue borders.

        :return: starting position of each residue
        :rtype:  list of int
        """
        result = []

        if self.lenAtoms() == 0:
            return N0.array( result, N0.Int )

        lastResNumber = -100
        lastResName   = ''
        index = -1
        lastAlt = 'x'
        lastSegid = -1

        res_nrs = self.atoms['residue_number']
        res_nam = self.atoms['residue_name']
        ins_cod = self.atoms['insertion_code']
        seg_id  = self.atoms['segment_id']

        ## create residue numbering for selected atoms
        for i in range( self.lenAtoms() ):

            if res_nrs[i] != lastResNumber or \
               res_nam[i] != lastResName   or \
               seg_id[i]  != lastSegid or \
               ins_cod[i] != lastAlt:

                ## start of new residue
                lastResNumber = res_nrs[i]
                lastResName   = res_nam[i]
                lastAlt       = ins_cod[i]
                lastSegid     = seg_id[i]
                index += 1

                result.append( i )

        return N0.array(result, N0.Int)        


    def resIndex( self, mask=None, force=0, cache=1 ):
        """
        Get the position of the each residue's first atom.

        :param force: re-calculate even if cached result is available (def 0)
        :type  force: 1||0
        :param cache: cache the result if new (def 1)
        :type  cache: 1||0
        :param mask: atom mask to apply before (i.e. result indices refer to \
                     compressed model)
        :type  mask: list of int (1||0)

        :return: index of the first atom of each residue
        :rtype: list of int
        """

        if self._resIndex is not None and not force and mask is None:
            return self._resIndex

        r = self.__inferResIndex()

        if mask is not None:
            m = self.index2map( r, len( mask ) )
            m = N0.compress( mask, m )
            r = self.map2index( m )

        if mask is None and cache:
            self._resIndex = r

        return r


    def resMap(  self, force=0, cache=1 ):
        """
        Get list to map from any atom to a continuous residue numbering
        (starting with 0). A new residue is assumed to start whenever the
        'residue_number' or the 'residue_name' record changes between 2
        atoms.

        See :class:`resList()` for an example of how to use the residue map.

        :param force: recalculate map even if cached one is available (def 0)
        :type  force: 0||1
        :param cache: cache new map (def 1)
        :type  cache: 0||1

        :return: array [00011111122223333..], residue index for each atom
        :rtype:  list of int
        """
        return self.index2map( self.resIndex( force=force,cache=cache ),
                               self.lenAtoms() )

    def resEndIndex( self ):
        """
        Get the position of the each residue's last atom.
        
        :return: index of the last atom of each residue
        :rtype: list of int
        """
        r = self.resIndex()
        return N0.concatenate( (r[1:], [self.lenAtoms()]) ) - 1 

    def __inferChainIndex( self ):

        result = []

        if self.lenAtoms() == 0:
            return N0.array( result, N0.Int )

        lastResidue = -100
        lastChainID = None
        lastSegID = None

        chn_ids = self.atoms['chain_id']
        seg_ids = self.atoms['segment_id']
        res_nrs = self.atoms['residue_number']
        ter_atm = self.atoms['after_ter']

        for i in self.atomRange():

            if chn_ids[i] != lastChainID or \
               seg_ids[i] != lastSegID   or \
               res_nrs[i] <  lastResidue or \
               ter_atm[i]:

                result.append( i )

            lastResidue = res_nrs[i]
            lastChainID = chn_ids[i]
            lastSegID   = seg_ids[i]

        return N0.array( result, N0.Int )


    def __filterSingleResChains( self, chainindex, ignore_resnumbers=0 ):
        """
        Join chains containing single residues with identical name into
        one chain. Typically required for waters or ions if they are
        separated by TER or picked up by the chain break detection
        
        :param check_resnumbers: (def 1)
        :type  check_resnumbers: 1||0
        """
        # residue name of first atom of each chain
        res_names = N0.take( self.atoms['residue_name'], chainindex )
        # residue number of first atom of each chain
        res_nmbrs = N0.take( self.atoms['residue_number'], chainindex )
        # chain id of first atom of each chain
        chain_ids = N0.take( self.atoms['chain_id'], chainindex )
        #segid of first atom of each chain
        seg_ids   = N0.take( self.atoms['segment_id'], chainindex )

        res_names = N0.concatenate( (['-1'], res_names) )
        chain_ids = N0.concatenate( (['-1'], chain_ids) )
        seg_ids   = N0.concatenate( (['-1'], seg_ids ) )
        res_nmbrs = N0.concatenate( ([-100], res_nmbrs) )

        delta     = res_nmbrs[1:] - res_nmbrs[:-1] 
        same_name = res_names[1:] == res_names[:-1]
        same_chain= chain_ids[1:] == chain_ids[:-1]
        same_seg  =   seg_ids[1:] ==   seg_ids[:-1]

        if ignore_resnumbers:
            delta = N0.ones( len(delta), N0.Int )

        is_single = (delta==1) \
                  * same_name * same_chain * same_seg

        return N0.compress( N0.logical_not(is_single), chainindex)


    def chainIndex( self, breaks=0, maxDist=None, force=0, cache=0,
                    singleRes=0, solvent=0 ):
        """
        Get indices of first atom of each chain.

        :param breaks: split chains at chain breaks (def 0)
        :type  breaks: 1||0
        :param maxDist: (if breaks=1) chain break threshold in Angstrom
        :type  maxDist: float
        :param force: re-analyze residue numbering, chain and segids to
                        find chain boundaries, use with care! (def 0)
        :type  force: 1||0
        :param cache: cache new index even if it was derrived from
                        non-default parameters (def 0)
                        **Note:** a simple m.chainIndex() will always cache
        :type  cache: 1||0
        :param singleRes: allow chains consisting of single residues (def 0)
                          Otherwise group consecutive residues with identical
                          name into one chain.
        :type  singleRes: 1||0
        :param solvent: also check solvent residues for "chain breaks" 
                        (default: false)
        :type  solvent: 1||0

        :return: array (1 x N_chains) of int
        :rtype: list of int
        """
        ## fast track
        if not (breaks or force or maxDist or solvent) \
           and self._chainIndex is not None:
            return self._chainIndex

        r = self._chainIndex

        if r is None or force:
            r = self.__inferChainIndex()

        if breaks:
            break_pos = self.chainBreaks( breaks_only=1, maxDist=maxDist,
                                          solvent=solvent, force=force )
            break_pos = break_pos + 1  ## chainBreaks reports last atom of each chain
            r = mathUtils.union( break_pos, r )
            r.sort()

        ## filter out chains consisting only of a single residue
        if len(r)>0 and not singleRes:
##            r = self.__filterSingleResChains( r, ignore_resnumbers=breaks )
            r = self.__filterSingleResChains( r, ignore_resnumbers=False)
            
        ## cache the result if it has been computed with default parameters
        if not(breaks or force or maxDist or singleRes or solvent) or cache:
            self._chainIndex = r

        return N0.array( r, N0.Int )

    def chainEndIndex( self, breaks=0, solvent=0 ):
        """
        Get the position of the each residue's last atom.
        
        :return: index of the last atom of each residue
        :rtype: list of int
        """
        r = self.chainIndex( breaks=breaks, solvent=solvent )
        return N0.concatenate( (r[1:], [self.lenAtoms()]) ) - 1 


    def chainMap( self, breaks=0, maxDist=None ):
        """
        Get chain index of each atom. A new chain is started between 2 atoms if
        the chain_id or segment_id changes, the residue numbering jumps back or
        a TER record was found.

        :param breaks: split chains at chain breaks (def 0)
        :type  breaks: 1||0
        :param maxDist: (if breaks=1) chain break threshold in Angstrom
        :type  maxDist: float

        :return: array 1 x N_atoms of int, e.g. [000000011111111111122222...]
        :rtype: list of int
        """
        return self.index2map( self.chainIndex( breaks=breaks, maxDist=maxDist ),
                               self.lenAtoms() )


    def chainBreaks( self, breaks_only=1, maxDist=None, force=0, solvent=0,
                     z=6. ):
        """
        Identify discontinuities in the molecule's backbone. By default,
        breaks are identified from the distribution of distances between the
        last backbone atom of a residue and the first backbone atom of the
        next residue. The median distance and standard deviation are
        determined iteratively and outliers (i.e. breaks) are identified
        as any pairs of residues with a distance that is more than z standard
        deviations (default 10) above the median. This heuristics can be
        overriden by specifiying a hard distance cutoff (maxDist).

        :param breaks_only: don't report ends of regular chains (def 1)
        :type  breaks_only: 1|0
        :param maxDist: maximal distance between consequtive residues
                        [ None ] .. defaults median + z * standard dev.
        :type  maxDist: float
        :param z      : z-score for outlier distances between residues (def 6.)
        :type  z      : float
        :param solvent: also check selected solvent residues (buggy!) (def 0)
        :type  solvent: 1||0
        :param force: force re-calculation, do not use cached positions (def 0)
        :type  force: 1||0

        :return: atom indices of last atom **before** a probable chain break
        :rtype: list of int
        """
        if self.__chainBreaks is not None and not force and \
           maxDist is None and breaks_only and not solvent and z==6.:
            r = self.__chainBreaks

        else:

            i_bb = N0.nonzero( self.maskBB( solvent=solvent ) )
            ## outlier detection only works with more than 2,
            ##  hard cutoff works with more than 1
            if len(i_bb) < 2:
                r = []
 
            else:
                bb   = self.take( i_bb )
                bb_ri= bb.resIndex()
                bb_re= bb.resEndIndex()
##                xyz = [ bb.xyz[ bb_ri[i] : bb_ri[i+1] ] for i in range(len(bb_ri)-1) ]
##                xyz +=[ bb.xyz[ bb_ri[-1]: len(bb) ] ]
##                centroid = N0.array([ N0.average( x ) for x in xyz ])

                last = N0.take( bb.xyz, bb_re )[:-1]
                first= N0.take( bb.xyz, bb_ri )[1:]
                
##                dist = N0.sqrt( N0.sum( N0.power(centroid[:-1]-centroid[1:],2), 
##                                      axis=1 ) )
                dist = N0.sqrt( N0.sum( N0.power(last-first,2), axis=1 ) )
                
                outliers, median, sd = mathUtils.outliers( dist, z=z, it=5 )
                
                ## get distances above mean
                cutoff = maxDist or median + z * sd
                r = N0.nonzero( N0.greater( dist, cutoff ) )
                
            if len(r) > 0:
                
##                ## can probably be simplified with self.resEndIndex()
##                ri = self.resIndex()
##                ri_to_e = {}
##                for i in range( len(ri)-1 ):
##                    ri_to_e[ ri[i] ] = ri[ i+1 ]-1
##    
##                ## map back to the original atom indices
##                r = [ ri_to_e[ i_bb[ bb_ri[i] ] ] for i in r ]

                ## this replacement didn't work out (fails PDBCleaner testcase):
                re = self.resEndIndex()
                r = [ re[i] for i in r ]


            if breaks_only:
                ri = self.chainIndex( breaks=0, solvent=solvent )
                r = [ x for x in r if not x+1 in ri ]

                if maxDist is None and not solvent and z==6.:
                    self.__chainBreaks = r

        return N0.array( r, int )


    def removeRes( self, what ):
        """
        Remove all atoms with a certain residue name.

        :param what: indices or name(s) of residue to be removed
        :type  what: str OR [ str ] OR int OR [ int ]
        """
        if not isinstance( what, list ) or isinstance( what, N.ndarray):
            what = T.toList( what )

        if type( what[0] ) is str:
            return self.remove( self.maskFrom( 'residue_name', what) )

        if type( what[0] ) is int:
            return self.remove( self.res2atomIndices( what ) )

        return False


    def rms( self, other, mask=None, mask_fit=None, fit=1, n_it=1 ):
        """
        Rmsd between two PDBModels.

        :param other: other model to compare this one with
        :type  other: PDBModel
        :param mask: atom mask for rmsd calculation
        :type  mask: list of int
        :param mask_fit: atom mask for superposition (default: same as mask)
        :type  mask_fit: list of int
        :param fit: superimpose first (default 1)
        :type  fit: 1||0
        :param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence, kicking out outliers on the way
        :type  n_it: int

        :return: rms in Angstrom
        :rtype: float
        """
        x, y = self.getXyz(), other.getXyz()

        if mask_fit is None: mask_fit = mask

        if fit:

            fx, fy = x, y

            if mask_fit is not None:
                fx = N0.compress( mask_fit, x, 0 )
                fy = N0.compress( mask_fit, y, 0 )
            ## find transformation for best match
            r, t = rmsFit.match( fx, fy, n_iterations=n_it )[0]

            ## transform coordinates
            y = N0.dot(y, N0.transpose(r)) + t

        if mask is not None:
            x = N0.compress( mask, x, 0 )
            y = N0.compress( mask, y, 0 )

        ## calculate row distances
        d = N0.sqrt(N0.sum(N0.power(x - y, 2), 1))

        return N0.sqrt( N0.average(d**2) )


    def transformation( self, refModel, mask=None, n_it=1,
                        z=2, eps_rmsd=0.5, eps_stdv=0.05,
                        profname='rms_outlier'):
        """
        Get the transformation matrix which least-square fits this model
        onto the other model.

        :param refModel: reference PDBModel
        :type  refModel: PDBModel
        :param mask: atom mask for superposition
        :type  mask: list of int
        :param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence
        :type  n_it: int        
        :param z: number of standard deviations for outlier definition
                  (default 2)
        :type  z: float
        :param eps_rmsd: tolerance in rmsd (default 0.5)
        :type  eps_rmsd: float
        :param eps_stdv: tolerance in standard deviations (default 0.05)
        :type  eps_stdv: float
        :param profname: name of new atom profile getting outlier flag
        :type  profname: str
        :return: array(3 x 3), array(3 x 1) - rotation and translation matrices
        :rtype: array, array
        """
        x, y = refModel.getXyz(), self.getXyz()

        if mask is not None:
            x = N0.compress( mask, x, 0 )
            y = N0.compress( mask, y, 0 )
            outlier_mask = N0.zeros( N0.sum(mask)  )

        else:
            outlier_mask = N0.zeros( len(self) )

        r, iter_trace = rmsFit.match( x, y, n_iterations=n_it, z=z,
                                      eps_rmsd=eps_rmsd, eps_stdv=eps_stdv)

        N0.put( outlier_mask, iter_trace[-1][-1], 1 )

        if n_it != 1:
            self.atoms.set( profname, outlier_mask, mask=mask,
                            default=1,
                            comment='outliers in last iterative fitting',
                            n_iterations=len( iter_trace ) )
        return r


    def transform( self, *rt ):
        """
        Transform coordinates of PDBModel.

        :param rt: rotational and translation array:
                   array( 4 x 4 ) OR array(3 x 3), array(3 x 1)
        :type  rt: array OR array, array

        :return: PDBModel with transformed coordinates
        :rtype: PDBModel
        """
        ## got result tuple from transformation() without unpacking
        if len( rt ) == 1 and type( rt[0] ) is tuple:
            rt = rt[0]

        if len(rt) == 2:
            r, t = rt[0], rt[1]
        else:
            rt = rt[0]
            r, t = (rt[0:3,0:3], rt[0:3, 3])

        result = self.clone()
        result.setXyz( N0.dot( self.getXyz(), N0.transpose(r) ) + t  )

        return result


    def fit( self, refModel, mask=None, n_it=1,
             z=2, eps_rmsd=0.5, eps_stdv=0.05,
             profname='rms_outlier'):
        """
        Least-square fit this model onto refMode

        :param refModel: reference PDBModel
        :type  refModel: PDBModel
        :param mask: atom mask for superposition
        :type  mask: list of int (1||0)
        :param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence
        :type  n_it: int        
        :param z: number of standard deviations for outlier definition
                  (default 2)
        :type  z: float
        :param eps_rmsd: tolerance in rmsd (default 0.5)
        :type  eps_rmsd: float
        :param eps_stdv: tolerance in standard deviations (default 0.05)
        :type  eps_stdv: float
        :param profname: name of new atom profile containing outlier flag
        :type  profname: str

        :return: PDBModel with transformed coordinates
        :rtype: PDBModel
        """
        return self.transform(
            self.transformation( refModel, mask, n_it, eps_rmsd=eps_rmsd,
                                 eps_stdv=eps_stdv, profname=profname ) )


    def magicFit( self, refModel, mask=None ):
        """
        Superimpose this model onto a ref. model with similar atom content.
        magicFit( refModel [, mask ] ) -> PDBModel (or subclass )

        :param refModel: reference PDBModel
        :type  refModel: PDBModel
        :param mask: atom mask to use for the fit
        :type  mask: list of int (1||0)

        :return: fitted PDBModel or sub-class
        :rtype: PDBModel
        """
        if mask is not None:
            m_this = self.compress( mask )
        else:
            m_this = self

        i_this, i_ref = m_this.compareAtoms( refModel )

        m_this = m_this.take( i_this )
        m_ref  = refModel.take(i_ref )

        ## find transformation for best match
        r,t = rmsFit.findTransformation( m_ref.getXyz(), m_this.getXyz() )

        result = self.transform( r, t )

        return result

    def structureFit( self, refModel, mask=None ):
        """
        Structure-align this model onto a reference model using the external
        TM-Align program (which needs to be installed).
        
        structureFit( refModel [, mask] ) -> PDBModel (or subclass)

        The result model has additional TM-Align statistics in its info record:
        r = m.structureFit( ref )
        r.info['tm_score'] -> TM-Align score
        the other keys are: 'tm_rmsd', 'tm_len', 'tm_id'
        
        .. seealso:: `biskit.TMAlign`

        :param refModel: reference PDBModel
        :type  refModel: PDBModel
        :param mask: atom mask to use for the fit
        :type  mask: list of int (1||0)

        :return: fitted PDBModel or sub-class
        :rtype: PDBModel
        """
        from biskit.exe import tmalign
        
        if mask is not None:
            m_this = self.compress( mask )
        else:
            m_this = self

        tm = tmalign.TMAlign( m_this, refModel )
        r = tm.run()

        return tm.applyTransformation( self )


    def centered( self, mask=None ):
        """
        Get model with centered coordinates.

        :param mask: atom mask applied before calculating the center
        :type  mask: list of int (1||0)

        :return: model with centered coordinates
        :rtype: PDBModel
        """
        r = self.clone()
        if mask is None: mask = N0.ones( len(self) )

        avg = N0.average( N0.compress( mask, r.getXyz(), 0 ) )

        r.setXyz( r.getXyz() - avg )

        return r


    def center( self, mask=None ):
        """
        Geometric centar of model.

        :param mask: atom mask applied before calculating the center
        :type  mask: list of int (1||0)

        :return: xyz coordinates of center
        :rtype: (float, float, float)
        """
        if mask is None:
            return N0.average( self.getXyz() )

        return N0.average( N0.compress( mask, self.getXyz(), axis=0 ) )


    def centerOfMass( self ):
        """
        Center of mass of PDBModel.

        :return: array(Float32)
        :rtype: (float, float, float)
        """
        M = self.masses()
        return mathUtils.wMean( self.getXyz(), M )


    def masses( self ):
        """
        Collect the molecular weight of all atoms in PDBModel.

        :return: 1-D array with mass of every atom in 1/12 of C12 mass.
        :rtype: array of floats

        :raise PDBError: if the model contains elements of unknown mass
        """
        try:
            M = [ molUtils.atomMasses[e] for e in self.atoms['element'] ]

        except KeyError as why:
            raise PDBError('Cannot find mass for '+str(why))

        return N0.array( M )


    def mass( self ):
        """
        Molecular weight of PDBModel.

        :return: total mass in 1/12 of C12 mass
        :rtype: float

        :raise PDBError: if the model contains elements of unknown mass
        """
        return N0.sum( self.masses() )


    def residusMaximus( self, atomValues, mask=None ):
        """
        Take list of value per atom, return list where all atoms of any
        residue are set to the highest value of any atom in that residue.
        (after applying mask)

        :param atomValues: values per atom
        :type  atomValues: list 
        :param mask: atom mask
        :type  mask: list of int (1||0)

        :return: array with values set to the maximal intra-residue value
        :rtype: array of float
        """
        if mask is None:
            mask = N0.ones( len(atomValues) )

        ## eliminate all values that do not belong to the selected atoms
        masked = atomValues * mask

        result = []

        ## set all atoms of each residue to uniform value
        for res in range( 0, self.lenResidues() ):

            ## get atom entries for this residue
            resAtoms = N0.compress( N0.equal( self.resMap(), res ), masked )

            ## get maximum value
            masterValue = max( resAtoms )

            result += resAtoms * 0.0 + masterValue

        return N0.array( result, N0.Float32 )


    def argsort( self, cmpfunc=None ):
        """
        Prepare sorting atoms within residues according to comparison function.

        :param cmpfunc: old style function(m.atoms[i], m.atoms[j]) -> -1, 0, +1
        :type  cmpfunc: function
        
        :param key: new style sort key function(m.atoms[i]) -> sortable
        :type key: function

        :return: suggested position of each atom in re-sorted model 
                 ( e.g. [2,1,4,6,5,0,..] )
        :rtype: list of int
        """
        ## cmp vanished in python 3.x (but still available in past.builtins)
        cmp = lambda x, y: (x > y) - (x < y)
        
        ## by default sort alphabetically by atom name
        cmpfunc = cmpfunc or ( lambda a1, a2: cmp(a1['name'],a2['name']) )

        result = []
        pos = 0

        ## get sort list for each residue
        for resAtoms in self.resList():

            resIndex = list(range(0, len( resAtoms )))

            ## convert atom-based function into index-based function
            f_cmp = lambda i,j,atoms=resAtoms : cmpfunc( atoms[i], atoms[j] )
            
            ## convert py 2.x cmp to py 3.x key method
            f_key = functools.cmp_to_key(f_cmp)

            ## get sortMap for this residue (numbering starts with 0)
            resIndex.sort( key=f_key )

            ## add first residue atom's position in self.atoms to each entry 
            resIndex = list(map( lambda i, delta=pos: i + delta, resIndex ))

            ## concatenate to result list
            result += resIndex

            pos += len( resAtoms )

        return result


    def sort( self, sortArg=None ):
        """
        Apply a given sort list to the atoms of this model.

        :param sortArg: comparison function
        :type  sortArg: function

        :return: copy of this model with re-sorted atoms (see numpy.take() )
        :rtype: PDBModel
        """
        sortArg = sortArg or self.argsort()

        r = self.take( sortArg )

        return r


    def unsort( self, sortList ):
        """
        Undo a previous sorting on the model itself (no copy).

        :param sortList: sort list used for previous sorting.
        :type  sortList: list of int

        :return: the (back)sort list used ( to undo the undo...)
        :rtype: list of int

        :raise PDBError: if sorting changed atom number
        """
        ## prepare sort functions (py 2.x / py 3.x)
        cmp = lambda x, y: (x > y) - (x < y)
        f_key = functools.cmp_to_key(lambda i,j, l=sortList: cmp(l[i],l[j]))
        
        ## get new sort list that reverts the old one
        backSortList = list(range(0, self.lenAtoms()))
        backSortList.sort( key=f_key )

        ## get re-sorted PDBModel copy
        self.keep( backSortList )

        return backSortList


    def atomNames(self, start=None, stop=None):
        """
        Return a list of atom names from start to stop RESIDUE index

        :param start: index of first residue 
        :type  start: int
        :param stop: index of last residue
        :type  stop: int

        :return: ['C','CA','CB' .... ]
        :rtype: list of str
        """
        ## By default return list of all atoms
        start = start or 0

        if stop is None:    ## don't use "stop = stop or ..", stop might be 0!
            stop = self.lenResidues()-1

        ## first get atom indexes of residues
        i = self.resIndex()[start]

        if stop == self.lenResidues()-1:
            j = self.lenAtoms()
        else:
            j = self.resIndex()[stop+1]

        return self.atoms['name'][i:j]


    def __testDict_and( self, dic, condition ):
        """
        Test if **all** key-value pairs of condition are matched in dic

        :param condition: {..}, key-value pairs to be matched
        :type  condition: dictionary
        :param dic: {..}, dictionary to be tested
        :type  dic: dictionary

        :return: 1|0, 1 if all key-value pairs of condition are matched in dic
        :rtype: 1|0
        """
        for k,v in list(condition.items()):
            if dic.get( k, None ) not in v:
                return 0
        return 1


    def __testDict_or( self, dic, condition ):
        """
        Test if **any** key-value pairs of condition are matched in dic

        :param condition: {..}, key-value pairs to be matched
        :type  condition: dictionary
        :param dic: {..}, dictionary to be tested
        :type  dic: dictionary

        :return: 1|0, 1 if any key-value pairs of condition are matched in dic
        :rtype: 1|0
        """
        for k,v in list(condition.items()):
            if dic.get( k, None ) in v:
                return 1
        return 0


    def filterIndex( self, mode=0, **kw ):
        """
        Get atom positions that match a combination of key=values.
        E.g. filter( chain_id='A', name=['CA','CB'] ) -> index

        :param mode: 0 combine with AND (default), 1 combine with OR
        :type  mode: 0||1
        :param kw: combination of atom dictionary keys and
                   values/list of values that will be used to filter
        :type  kw: filter options, see example

        :return: sort list
        :rtype: list of int
        """
        ## cache to minimize function lookup
        atoms = self.atoms.toDicts()
        if mode == 0:
            f_test = self.__testDict_and
        else:
            f_test = self.__testDict_or

        for k in kw:
            kw[ k ] = T.toList( kw[ k ] )

        r = [ i for i in range(self.lenAtoms()) if f_test( atoms[i], kw ) ]
        return r


    def filter( self, mode=0, **kw):
        """
        Extract atoms that match a combination of key=values.
        E.g. filter( chain_id='A', name=['CA','CB'] ) -> PDBModel

        :param mode: 0 combine with AND (default), 1 combine with OR
        :type  mode: 0||1
        :param kw: combination of atom dictionary keys and
                   values/list of values that will be used to filter
        :type  kw: filter options, see example

        :return: filterd PDBModel
        :rtype: PDBModel
        """
        return self.take( self.filterIndex( mode=mode, **kw ) )


    def equals(self, ref, start=None, stop=None):
        """
        Compares the residue and atom sequence in the given range.
        Coordinates are not checked, other profiles are not checked.

        :param start: index of first residue 
        :type  start: int
        :param stop: index of last residue
        :type  stop: int 

        :return: [ 1||0, 1||0 ],
                 first position sequence identity 0|1,
                 second positio atom identity 0|1
        :rtype: list if int
        """
        ## By default compare all residues
        start = start or 0
        if stop is None:    ## don't use stop = stop or .. stop might be 0!
            stop = self.lenResidues()

        ## set to 1 when identical
        seqID, atmID = 0, 0

        ## compare sequences
        if self.sequence()[start:stop] == ref.sequence()[start:stop]:
            seqID = 1

            ##  then compare atoms
            if self.atomNames(start,stop-1) == ref.atomNames(start,stop-1):
                atmID = 1

        return [seqID, atmID]


    def compareAtoms( self, ref ):
        """
        Get list of atom indices for this and reference model that converts
        both into 2 models with identical residue and atom content.

        E.g.
         >>> m2 = m1.sort()    ## m2 has now different atom order
         >>> i2, i1 = m2.compareAtoms( m1 )
         >>> m1 = m1.take( i1 ); m2 = m2.take( i2 )
         >>> m1.atomNames() == m2.atomNames()  ## m2 has again same atom order

        :return: indices, indices_ref
        :rtype: ([int], [int])
        """
        ## compare sequences
        seqMask, seqMask_ref = match2seq.compareModels(self, ref)

        ## get list of matching RESIDUES
        equal = N0.nonzero(seqMask)
        equal_ref = N0.nonzero(seqMask_ref)

        result, result_ref = [], []
        
        rI = self.resIndex()
        rIref = ref.resIndex()

        ## check that all atoms are equal in matching residues
        for i in range(0, len(equal)):

            ## atom name lists for current residue
            aa = self.atomNames( equal[i],equal[i] )
            aa_ref = ref.atomNames( equal_ref[i],equal_ref[i] )

            ## starting atom of current residue
            ind = rI[ equal[i] ]
            ind_ref = rIref[ equal_ref[i] ]

            for j in range( len(aa_ref) ):

                try:
                    ##shortcut for mostly equal models
                    if aa_ref[j] == aa[j]:      ## throws IndexError
                        result     += [ind + j]
                        result_ref += [ind_ref + j]
                        continue

                except IndexError:
                    pass

                try:
                    pos = aa.index( aa_ref[j] ) ## throws ValueError

                    result     += [ind + pos]
                    result_ref += [ind_ref + j]

                except ValueError:
                    pass

        return result, result_ref

    def unequalAtoms( self, ref, i=None, iref=None ):
        """
        Identify atoms that are not matching between two models. 
        This method returns somewhat of the opposite of compareAtoms().
        
        Not matching means: (1) residue is missing, (2) missing atom within a
        residue, (3) atom name is different. Differences in coordinates or 
        other atom profiles are NOT evaluated and will be ignored.
        
        (not speed-optimized)
        
        :param ref: reference model to compare to
        :type  ref: PDBModel
        :param i:   pre-computed positions that are equal in this model \
                     (first value returned by compareAtoms() )
        :type  i:   array( int ) or [ int ]
        :param iref: pre-computed positions that are equal in ref model \
                     (first value returned by compareAtoms() )
        :type  i:    array( int ) or [ int ]
        
        :return: missmatching atoms of self, missmatching atoms of ref
        :rtype: array(int), array(int)
        """
        if i is None or iref is None:
            i, iref = self.compareAtoms( ref )
        mask_self = N0.ones( len(self), int )
        mask_ref  = N0.ones( len(ref ), int )
        
        N0.put( mask_self, i, 0 )
        N0.put( mask_ref, iref, 0 )
        
        return N0.nonzero( mask_self ), N0.nonzero( mask_ref )
    

    def reportAtoms( self, i=None, n=None ):
        """
        :param i: optional list of atom positions to report (default: all)
        :type i:  [ int ]
        :return: formatted string with atom and residue names similar to PDB
        :rtype:  str
        """
        m = self
        n = n or len(m)
        
        if i is not None:
            m = self.take( i )
        
        s = '%(serial_number)4i %(name)5s %(residue_name)3s %(chain_id)1s '+\
            '%(residue_number)3i'

        atm = [ s % a for a in m ]
        r = '\n'.join( atm[:n] )
        if n < len( m ):
            r += ' ...'
        return r

    def __chainFraction( self, chain, ref ):
        """
        Look how well a given chain matches a continuous stretch of residues
        in ref.

        :param chain: chain index
        :type  chain: int
        :param ref: reference PDBModel
        :type  ref: PDBModel

        :return: average relative length of matching chain fragments
        :rtype: float
        """
        m = self.takeChains([chain])
        
        if len(m) == 0:
            return 0

        m_cast = m.take( m.compareAtoms( ref )[0] )

        f = 1. * len( m_cast )
        if f > 0:
            f = f / m_cast.lenChains( breaks=1, maxDist=5.)

        f = f / len( m )

        return f


    def compareChains( self, ref, breaks=0, fractLimit=0.2 ):
        """
        Get list of corresponding chain indices for this and reference model.
        Use takeChains() to create two models with identical chain content and
        order from the result of this function.

        :param ref: reference PDBModel
        :type  ref: PDBModel
        :param breaks: look for chain breaks in backbone coordinates
        :type  breaks: 1||0
        :param fractLimit: 
        :type  fractLimit: float

        :return: chainIndices, chainIndices_ref
        :rtype: ([int], [int])
        """
        i, i_ref = self.compareAtoms( ref )

        c0  = self.atom2chainIndices( i, breaks=breaks )
        c_r = ref.atom2chainIndices( i_ref, breaks=breaks )

        ## dirty hack to throw out single matching residues
        c0 = [ c for c in c0 \
               if self.__chainFraction( c, ref ) > fractLimit  ]
        c_r= [ c for c in c_r \
               if ref.__chainFraction( c, self ) > fractLimit  ]

        return c0, c_r

    def biomodel(self, assembly = 0):
        """
        Return the 'biologically relevant assembly' of this model
        according to the information in the PDB's BIOMT record (captured in
        info['BIOMT']). 
        
        This removes redundant chains and performs symmetry operations to
        complete multimeric structures. Some PDBs define several alternative
        biological units: usually (0) the author-defined one and (1)
        software-defined -- see :class:`lenBiounits`.

        Note: The BIOMT data are currently not updated during take/compress
        calls which may change chain indices and content. This method is
        therefore best run on an original PDB record before any other
        modifications are performed.
        
        :param assembly: assembly index (default: 0 .. author-determined unit)
        :type  assembly: int

        :return: PDBModel; biologically relevant assembly
        """
        try:
            biounit = BU.BioUnit(self, self.info['BIOMT'])
            r = biounit.makeMultimer(assembly)
        except AttributeError:
            r = self
        return r

    def lenBiounits (self):
        """
        Number of biological assemblies defined in PDB BIOMT record, if any.
        
        :return: number of alternative biological assemblies defined in
                 PDB header
        :rtype: int
        """
        try:
            biounit = BU.BioUnit(self, self.info['BIOMT'])            
            r = len(list(biounit.keys()))
        except AttributeError:
            r = 0
        return r

    def atomkey( self, compress=True ):
        """
        Create a string key encoding the atom content of this model independent
        of the order in which atoms appear within residues. Atom names are
        simply sorted alphabetically within residues and then concatenated.
        
        :param compress: compress key with zlib (default: true)
        :type  compress: bool
        :return: key formed from sorted atom content of model
        :rtype: str
        """
        import zlib
        rindex = N0.concatenate( (self.resIndex(), [len(self)] ) )
        r = ''
        if len(rindex) == 0:
            return r
        
        anames = self.atoms['name'] ## cache for faster access
        
        for i in range(len(rindex)-1):
            a = anames[rindex[i]:rindex[i+1]]
            a.sort()
            a = ''.join(a)
            r = r + a
            
        if compress:
            return zlib.compress(r)
            
        return r
  


#############
##  TESTING        
#############

import biskit.test as BT

class _TestData(object):
    MODEL = None

class Test(BT.BiskitTest):
    """Test class """

    #: load test PDB once into class rather than 3 times into every instance
    ## for some reason doesn't actually work
    MODEL = None

    def prepare(self):
        import tempfile

        ## loading output file from X-plor
        self.MODEL = _TestData.MODEL or B.PDBModel( T.testRoot()+'/com/1BGS.pdb')
        _TestData.MODEL = self.MODEL
        self.m = self.MODEL
        self.fout_pdb = tempfile.mktemp( '_test1.pdb' )
        self.fout1 = tempfile.mktemp( '_test1.model' )
        self.fout2 = tempfile.mktemp( '_test2.model' )

    def cleanUp( self ):
        T.tryRemove( self.fout1 )
        T.tryRemove( self.fout2 )
        T.tryRemove( self.fout_pdb )

    def test_removeRes(self):
        """PDBModel.removeRes test"""
        t = time.time()
        self._m = self.m.clone()
        self._m.removeRes(['TIP3', 'HOH'])
        self.assertEqual( len(self._m), 1968)
        self.assertAlmostEqual( self._m.mass(), 21325.90004, 3 )
        if self.local: print("removeRes: ", time.time() - t)

    def test_chainMethods(self):
        """PDBModel chain methods test"""
        ## X-plor doesn't write chainIds, so during the simulation
        ## we store them in the last letter of the segId. Here we
        ## restore the chainId.
        self._m = self.m.clone()
        self._m.addChainFromSegid()

        ## start positions of all chains
        chainIdx = self._m.chainIndex().tolist()

        ## print some chain info
        if self.local:
            print('The molecule consists of %i chains'% self.m.lenChains())
            print('\tChainId \tFirst atom')
            for i in chainIdx:
                print('\t%s \t\t%i'%(self._m.atoms['chain_id'][i], int(i)))

        ## iterate over all chains
        for c in range( 0, len( chainIdx ) ):

            if self.local:
                print("chain ", c, " starts with ", end=' ') 
                print(self._m.atoms['residue_name'][ chainIdx[c] ], end=' ')

                print(" and has sequence: ")


            ## mask out atoms of all other chains
            chainMask  = N0.equal( self._m.chainMap( breaks=1 ), c )
            if self.local:
                print(self._m.sequence( chainMask ))

        self.assertEqual( self._m.lenChains(), 4)


    def test_sorting(self):
        """PDBModel sorting test"""
        if self.local:
            print("sorting atoms alphabetically...")
        m2 = self.m.compress( self.m.maskProtein() )
        sort = m2.argsort()
        m2 = m2.sort( sort )

        self.assertAlmostEqual( N0.sum( m2.centerOfMass() ),  23.1032009125,2)


    def test_chainBreaks(self):
        """PDBModel chain break handling and writing test"""
        self.m4 = B.PDBModel( T.testRoot()+'/com/1BGS_original.pdb')
        self.assertEqual( self.m4.lenChains(), 9 )
        self.assertEqual( self.m4.lenChains( breaks=1 ), 9 )
        self.assertEqual( self.m4.lenChains( breaks=1, singleRes=1, solvent=1), 
                          9 )
        self.m4.writePdb( self.fout_pdb, ter=2 )
        
    def test_chainBreaks2(self):
        """PDBModel more complicated chain break detection"""
        self.m5 = B.PDBModel( T.testRoot()+'/pdbclean/foldx_citche.pdb')
        breaks = self.m5.chainBreaks()
        self.assertEqual( len(breaks), 2 )
        
        ## limitation of the method: same model but now with capping residues
        ## filling the gap
        self.m6 = B.PDBModel( T.testRoot()+'/pdbclean/citche_capped.pdb')
        breaks = self.m6.chainBreaks()
        self.assertEqual( len(breaks), 1 )

    def test_chainSingleResidues( self ):
        """PDBModel single residue chain test"""
        self.m5 = B.PDBModel( T.testRoot() + '/amber/1HPT_0.pdb' )
        self.assertTrue( self.m5.lenChains() < 10, 'single residue chains' )


    def test_rename(self):
        """PDBModel renameAmberRes tests"""
        self.m3 = B.PDBModel( T.testRoot()+'/amber/leap/1HPT_dry.pdb')

        n_cyx = self.m3.atoms['residue_name'].count('CYX')
        n_hid = self.m3.atoms['residue_name'].count('HID')
        n_hip = self.m3.atoms['residue_name'].count('HIP')
        n_hie = self.m3.atoms['residue_name'].count('HIE')
        n_hix = n_hid + n_hie + n_hip

        self.m3.renameAmberRes()

        self.assertEqual(n_cyx, self.m3.atoms['residue_name'].count('CYS'))
        self.assertEqual(n_hix, self.m3.atoms['residue_name'].count('HIS'))

    def test_xplor2amber(self):
        """PDBModel xplor2amber test"""
        ## test is simply back-converting a PDB comming from 'ambpdb -aatm'
        ## a better test would be starting from an actually xplor-generated PDB
        m1 = B.PDBModel( T.testRoot() +'/amber/1HPT_0dry.pdb' )
        m1.xplor2amber( aatm=True )

        m2 = B.PDBModel( T.testRoot() + '/amber/1HPT_0dry_amberformat.pdb' )
        self.assertEqual( m1.atomNames(), m2.atomNames() )


    def test_report(self):
        """PDBModel report&plot test"""
        self.report_output = self.m.report( prnt=self.local,
                                            plot=(self.local or self.VERBOSITY > 2) )


    def test_sourceHandling(self):
        """PDBModel source / disconnection tests"""
        self._m = self.m.clone()
        anames = self.m.atoms['name']
        xyz0   = self.m.getXyz()[0]

        self._m.slim()

        ## _m2 uses _m1 as source
        self._m2 = B.PDBModel( self._m )
        l1 = self._m2.atoms['name']
        self.assertEqual( l1, anames )

        ## remove unchanged profiles and coordinates
        self._m2.slim()

        ## fetch them again from source (of source)
        l2 = self._m2.atoms['name']
        self.assertEqual( l2, anames )

        ## disconnect _m from PDB file source
        self._m.saveAs( self.fout2 )

        self._m2.slim()
        self._m.slim()

        ## this should now trigger the reloading of fout2
        self.assertEqual( self._m2.atoms['name'], anames )
        self.assertTrue( N.all( self._m2.getXyz()[0] == xyz0) )

        ## after disconnection, slim() should not have any effect
        self._m2.disconnect()
        self._m2.slim()
        self.assertTrue( self._m2.atoms.profiles['name'] is not None )

    def test_mergeChains( self ):
        """PDBModel.mergeChains test"""
        m = self.m.takeChains( [0] )
        res_numbers = m['residue_number']
        atm_numbers = m['serial_number']
        chain_ids   = m['chain_id']

        m1 = m.takeResidues( list(range(3)) )
        m2 = m.takeResidues( list(range(3, m.lenResidues())) )

        m2.renumberResidues()
        m2['chain_id'] = len(m2) * ['X']

        self.m1 = m1
        self.m2 = m2

        self.r = m1.concat( m2 )
        r = self.r
        self.assertTrue( r.lenChains() == m.lenChains() + 1 )

        r.mergeChains( 0 )
        self.r = r
        self.assertTrue( r.lenChains() == m.lenChains() )
        self.assertTrue( N.all( N0.array(r['chain_id']) == chain_ids ) )
        self.assertTrue( N.all( N0.array(r['residue_number']) == res_numbers ) )

    def test_mergeResidues( self ):
        """PDBModel.mergeResidues test"""
        m = self.m.clone()
        gg_position = m.sequence().find( 'GG' )

        len_r = m.lenResidues()
        len_a = m.lenAtoms()

        r_gly = Residue( m, gg_position )

        m.mergeResidues( gg_position )
        r_gly.reset()

        self.assertEqual( m.lenResidues(), len_r - 1 )
        self.assertEqual( m.lenAtoms(), len_a )
        self.assertEqual( len( r_gly ), 2 * 5 )

    def test_getset(self):
        """PDBModel.__get/set__ test"""
        self.assertEqual( self.m[10]['occupancy'], 1.0 )
        self.assertEqual( self.m['chain_id', 'changed'], 0 )
        self.assertEqual( len(self.m['chain_id'] ), len( self.m ) )
        self.assertTrue( type( self.m['date']) is str )
        self.m['resname'] = self.m.atom2resProfile('residue_name')
        self.assertEqual( len( self.m['resname'] ), self.m.lenResidues() )

        self.m.info['tested'] = False
        self.m['tested'] = True
        self.assertTrue( self.m.info['tested'] )

        self.m['serial_number', 'default'] = 1
        self.assertTrue( self.m.atoms['serial_number','default'] == 1 )

        self.m['resname', 'changed'] = 0
        self.assertFalse( self.m.residues.isChanged( 'resname' ) )

        self.m['index'] = list(range( len( self.m)))
        self.assertTrue( self.m['index'][-1] == len( self.m ) - 1 )

    def test_slice(self):
        """PDBModel.__slice__ test"""
        self.assertTrue( len( self.m[0:100:20]  ) == 5 )

    def test_various(self):
        """PDBModel various tests"""
        m = PDBModel()
        self.assertEqual( type( m.getXyz() ), N.ndarray )

    def test_compareChains(self):
        """PDBModel.compareChains test"""
        m = self.m.clone()
        m2 = PDBModel()
        ## extract first 100 atoms of each chain
        for i in range(m.lenChains()):
            m2 = m2.concat( m.takeChains([i]).take( list(range(100)) ) )

        m3 = m2.takeChains( [2,3,0,1] )  ## re-order chains

        i, iref = m2.compareChains( m3 )



class TestExe( BT.BiskitTest ):
    """PDBModel tests that rely on external applications"""

    TAGS = [BT.EXE]

    def test_structureFit( self ):
        """PDBModel.structureFit test"""
        m = T.load( T.testRoot( 'tmalign/1huy_citrine.model' ) )
        ref = T.load( T.testRoot( 'tmalign/1zgp_dsred_dimer.model' ) )
        ref = ref.takeChains( [0] )

        r = m.structureFit( ref )
        diff = r.centerOfMass() - ref.centerOfMass()

        if self.local:
            print('center of mass deviation: \n%r' % diff)

        self.assertEqual( r.info['tm_rmsd'], 1.76 )
        self.assertTrue( N.all( N0.absolute(diff) < 1 ),
                      'superposition failed: %r' % diff)


def clock( s, ns=globals() ):  ## pragma: no cover
    import cProfile

    locals().update( ns )

    r = cProfile.run( s, 'report.out' )

    ## Analyzing
    import pstats
    p = pstats.Stats('report.out')
    p.strip_dirs()

    ## long steps and methods calling them
    p.sort_stats('cumulative').print_stats(20)
    p.print_callers( 20 )

    return r

if __name__ == '__main__':

    BT.localTest()   
