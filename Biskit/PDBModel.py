## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

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
## last $Author$
## last $Date$
## $Revision$

"""
Store and manipulate coordinates and atom information.
"""

import tools as T
import molUtils
import mathUtils
import match2seq
import rmsFit
from LocalPath import LocalPath
from Errors import BiskitError
from Biskit import EHandler
from ProfileCollection import ProfileCollection, ProfileError
from PDBParserFactory import PDBParserFactory
from PDBParseFile import PDBParseFile
import Biskit as B

import numpy.oldnumeric as N
import numpy.oldnumeric.mlab as MLab
import os, sys
import copy
import time
import string
import types
import Scientific.IO.PDB as IO


class PDBProfiles( ProfileCollection ):
    """
    A ProfileCollection that triggers an update() of its parent L{PDBModel}
    if an empty or (optionally) missing profile is requested.

    @see L{ProfileCollection}
    """

    def __init__(self, model=None, profiles=None, infos=None ):
        """
        @param model: parent model of this ProfileCollection
        @type  model: PDBModel
        @param profiles: dictionary of existing profiles
        @type  profiles: { 'name' : list/array }
        @param infos:    dictionary of existing meta infos
        @type  infos:    { 'name' : { 'date' : ... } }
        """
        ProfileCollection.__init__( self, profiles=profiles, infos=infos )
        self.model = model

    def version( self ):
        return ProfileCollection.version(self) + '; PDBModel $Revision$'


    def get( self,  name, default=None, update=True, updateMissing=False ):
        """
        Fetch a profile::
          get( profKey, [default] ) -> list or array of values 
        Or:
          get( (profKey, infoKey), [default] ) -> single value of metainfo dict
        
        This method extends the standard L{ProfileCollection.get} by the
        ability to fetch empty (None) or missing profiles from a source
        file or model.

        @param name: profile key or profile and info key
        @type  name: str OR (str, str)
        @param default: default result if no profile is found,
                        if None and no profile is found, attempt update
        @type  default: any
        @param update: update from source before returning empty profile [0]
        @type  update: bool
        @param updateMissing: update from source before reporting missing
                              profile [0]
        @type  updateMissing: bool
        
        @raise ProfileError: if no profile is found with |name|
        """
        try:
            r = ProfileCollection.get( self, name, default=default)

        except ProfileError, e:
            if updateMissing:
                r = None
            else:
                raise ProfileError, e

        if r is None and (update or updateMissing):

            ## only read PDB source if this is indeed useful
            if not name in PDBModel.PDB_KEYS and \
               PDBParseFile.supports( self.model.validSource() ):
                return None

            self.model.update( updateMissing=updateMissing )

            ## try again
            r = ProfileCollection.get(self, name )

        return r


class PDBError(BiskitError):
    pass

class PDBModel:
    """ 
    Store and manipulate coordinates and atom infos stemming from a PDB
    file. Coordinates are stored in the Numeric array 'xyz'; the additional
    atom infos from the PDB (name, residue_name, and many more) are stored in
    a L{PDBProfiles} instance 'aProfiles' which can be used to also associate
    arbitrary other data to the atoms. Moreover, a field 'rProfiles' can hold
    data associated to residues. A normal dictionary 'info' is there to accept
    any information about the whole model.

    Creating a PDBModel
    ===================
 
         1. from a PDB file::

           >>> m = PDBModel('~/my/structure.pdb')

         2. from a pickled PDBModel::

           >>> m = PDBModel('~/my/pickled.model')

          Or::
           >>> import Biskit.tools as T
           >>> m = T.Load('~/my/pickled.model')

         3. from another PDBModel in memory::

           >>> m = PDBModel( sourcemodel )

      The file or object passed to the constructor becomes the L{source}
      of the PDBModel. See the section on saving PDBModels below!


    Working with profiles
    =====================

      See L{ProfileCollection} and L{CrossView}!
      
      For example:

        >>> m.aProfiles['element']

      ...will give you a list of length m.lenAtoms() with the Element character
      of each atom.
      
        >>> m.aProfiles
      and
        >>> m.rProfiles
        
      ... will give you an overview over the available atom profiles and their
      associated meta data.

      ProfileCollection (the mother class of aProfiles and rProfiles)
      offers many convenient methods to access, modify, iterate
      and even plot profiles. 
      
      For convenience, we have defined some shortcuts so that many profile
      methods can be directly accessed from the model. In the following
      examples, m is again a PDBModel instance::

        >>> m['name']            ## fetch the atom profile 'name'
        ['H1', 'H2', 'N', 'H3', 'CA', 'CB', 'C', 'O', 'N', 'H', '...

        >>> m['residue_name']   ## fetch the **atom** profile 'residue_name'
        ['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'GLN', ...

      Note that 'residue_name' still has one value per atom because it comes
      from the parsing of the PDB file. If you want to have a real residue
      profile with the residue names, you can create a copy::

        >>> m['resname'] = m.atom2resProfile('residue_name')

      Which is equivalent to the direct call::

        >>> m.rProfiles.set('resname', m.atom2resProfile( 'residue_name' ))

      Now let's do some simple calculation. We want to create a profile
      that lists the distance of each atom from the molecule's center
      of mass. This is how it works::

        >>> diff = m.xyz - m.centerOfMass()
        >>> dist = N.sqrt( N.sum( diff**2, axis=1 ) )
        >>> m['from_center'] = dist

      However, the shortcuts don't offer you all the options of the real
      ProfileCollection methods. So you may rather want to create the
      profile like this::

        >>> m.aProfiles.set('from_center', dist, 'my_bday'='01/01/1900',
                            comment='distance in A from center of mass' )

      You can then access the meta infos of your profile like this::
        >>> m['from_center', 'my_bday']
        '01/01/1900'

      Sometimes you don't actually have a value for each atom. So suppose,
      you extracted only C-alpha distances from your dist profile. You
      can still use them as atom profile if you give a mask designating
      the positions for which you actually have values::

        >>> dist = N.compress( dist, m.maskCA() )  ## extract CA distances
        >>> m.aProfiles.set('ca_from_center', dist, mask=m.maskCA(),
                            default=-1, comment='CA distances' )

      The new profile will have the CA distances associated to CA atoms
      and -1 in all other positions.

    Working with atoms
    ==================

      You can fetch a virtual dictionary with the values of all atom
      profiles for one atom::
      
        >>> m[0]
        CrossView{'name': 'H1', 'residue_number': 1, 'insertion_code': '',
        ...

      Which is the same as::
        >>> m.aProfiles[0]
      
      Changes to this dictionary change the underlying profiles at the
      position of this atom. Vice versa, changes to the profiles are
      mirrored in the 'dictionary'.

      Iteration over a model m will give you one CrossView after the other::
        >>> for atom in m:
        ...    if atom['element']=='C' and atom['temperature_factor'] > 70:
        ...        atom['charge'] = 1.

      However, working with CrossViews is usually much slower than directly
      working with the underlying profiles. So the above example
      would be much more efficiently done like this::

        >>> m['charge'] = N.greater( m['temperature_factor'], 70 ) * \
        ...                  m.maskFrom('element', 'C')
      
      See the description of L{CrossView} for a detailed discussion of
      performance issues.

    Manipulating PDBModels
    ======================

      Methods are provided to extract, remove, re-arrange, or concatenate
      atoms, chains or residues in ways similar to the handling of Numeric
      arrays. The major methods for this are L{take}, L{compress}, L{concat},
      L{takeResidues} and L{takeChains}. All of these methods operate on
      coordinates, atom- and residue profiles simultaneously so that
      everything always stays synchronized.

      Atoms, residues and chains can be selected by their indices (positions)
      or with masks (of 1 and 0 or True and False). And there are many methods
      to produce such masks. The following example extracts all CA atoms from
      a PDBModel and puts them into a new model 'm_ca'::

        >>> mask_ca = m.maskFrom('name', 'CA')
        >>> m_ca = m.compress( mask_ca )

      We have a shortcut for this::
        >>> m_ca = m.compress( m.maskCA() )

      Masks are normal Numpy arrays and can be combined with logical_and,
      logical_or (simple * and + does the same) and other methods. Masks
      calculated from atoms can be converted to residue or chain masks and
      back with atom2resMask, atom2chainMask, res2atomMask, etc. methods.

      Many PDBModel methods (called take-something) operate
      with indices rather than masks.  Indices are lists of (atom,
      residue, or chain) positions. A mask can easily be converted
      into a list of indices with the Numeric.nonzero() function. So::

        >>> m_ca = m.compress( mask_ca )
      Is equivalent to::
        >>> m_ca = m.take( N.nonzero( mask_ca ) )
      
      Like masks, indices can be converted between atom, residue, and
      chain level with atom2resIndices, res2atomIndices, and so on. The
      following example builds a new model m2 from only the first and last
      chain of an existing model m::

        >>> m2 = m.takeChains( [0] ).concat( m.takeChains( [-1] ) )
    
      However more efficient and easier would be::
        >>> m2 = m.takeChains( [0,-1] )

      Which is the same as::
      
        >>> m2 = m.take( m.chain2atomIndices( [0,-1] ) )

      numpy.put() is used to convert a list of indices into a mask::

        >>> cacb_pos = m.indicesFrom('name', ['CA','CB'] )
        >>> mask = N.zeros( len(m) )
        >>> N.put( mask, cacb_pos, 1 )
        >>>   ## Test:
        >>> N.all( mask == N.nonzero( cacb_pos ) )
        True

    Comparing PDBModels
    ===================

      You can compare and normalize the atom content of similar models
      with  L{compareAtoms}, calculate an (iterative) RMSD with
      L{rms}, and superposition two equivalent (L{fit}) or almost equivalent
      models (L{magicFit}) -- again with or without iterative removal of
      variable regions.

    Surfaces, Energies, Conservations...
    ====================================

      L{PDBDope} adds surfaces, conservation profiles, energies and other
      information from external programs to PDBModels::

        >>> d = PDBDope( m )
        >>> d.addConservation()  ## needs Hmmer and Pfam installed
        >>> d.addSurfaceRacer()  ## access. and mol. surface from surfrace
        >>> d.addDensity()       ## simple atomic density
        >>> d.addFoldX()         ## puts fold-X scores into m.info
        >>> d.addSecondaryStructure()  ## from DSSP
        
    Saving a PDBModel
    =================

      A normal PDB file can be created from the model at any time::
        >>> m.writePdb('~/my/modified.pdb')
      There are several options for dealing with TER records and
      atom names.

      However, writing and reading PDB files is rather slow and
      additional information in rProfiles, aProfiles, or info is lost.
      We therefore recommend to pickle PDBModels as python objects
      to disc. This can be done in two ways::

        >>> T.Dump( m, '~/my/repickled.model' )
      Or::
        >>> m.saveAs( '~/my/repickled.model' )

      There is a subtle but important difference between the two
      methods: The model remembers its source (a PDB file or a
      PDBModel in memory or a pickled PDBModel on disc) and keeps
      track of whether xyz, or (which) profiles have been changed with
      respect to this source. Normal pickling with T.Dump will
      automatically trigger a L{slim()} method. slim() resets the xyz
      array or any profiles to 'None' if they have not been changed
      since being read from a source on disc (a source in memory
      doesn't count). 

      However, the source connection can also become a little
      confusing, especially if you move or change the file with the
      source model. To avoid any such confusion you should use the
      second method for saving models::

        >>> m.saveAs( '~/my/repickled.model' )
      
      saveAs produces an independent stand-alone PDBModel on
      disc, which can also be copied around between different
      machines (python pickles are, at least in theory, platform
      independent).

    Relation to other Biskit classes
    ================================

      A L{Dock.Complex} can be created from two PDBModels and offers
      methods for contact evaluation, scoring and other things::

        >>> c = Complex( rec_model, lig_model )

      The receptor and ligand PDBModels of a complex are accessible
      as::
        >>> m1 = c.rec()
        >>> m2 = c.lig()

      A L{Trajectory} adds a time dimension to a PDBModel and can be
      created from a list of PDBModels, PDBs or an Amber crd file. The
      underlying PDBModel of a trajectory t is accessed as::
        >>> m = t.ref
      And a PDBModel for a certain time point can be created by indexing::
        >>> model_10 = t[10]
    """

    #: keys of all atom profiles that are read directly from the PDB file
    PDB_KEYS = ['name', 'residue_number', 'insertion_code', 'alternate',
                'name_original', 'chain_id', 'occupancy', 'element',
                'segment_id', 'charge', 'residue_name', 'after_ter',
                'serial_number', 'type', 'temperature_factor']

    def __init__( self, source=None, pdbCode=None, noxyz=0, skipRes=None ):
        """
        - PDBModel() creates an empty Model to which coordinates (field xyz)
          and PDB infos (field atoms) have still to be added.
        - PDBModel( file_name ) creates a complete model with coordinates
          and PDB infos taken from file_name (pdb, pdb.gz, pickled PDBModel)
        - PDBModel( PDBModel ) creates a copy of the given model
        - PDBModel( PDBModel, noxyz=1 ) creates a copy without coordinates

        @param source: str, file name of pdb/pdb.gz file OR pickled PDBModel OR
                      PDBModel, template structure to copy atoms/xyz field from
        @type  source: str or PDBModel
        @param pdbCode: PDB code, is extracted from file name otherwise
        @type  pdbCode: str or None
        @param noxyz: 0 (default) || 1, create without coordinates
        @type  noxyz: 0||1

        @raise PDBError: if file exists but can't be read
        """
        self.source = source
        if type( self.source ) is str:
            self.source = LocalPath( self.source )

        self.__validSource = 0
        self.fileName = None
        self.pdbCode = pdbCode
        self.xyz = None

        #: save atom-/residue-based values
        self.rProfiles = PDBProfiles( self )
        self.aProfiles = PDBProfiles( self )

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
        self.initVersion = self.version()

        #: to collect further informations
        self.info = { 'date':T.dateSortString() }

        if source <> None:

            self.update( skipRes=skipRes, updateMissing=1 )

        if noxyz:
            ## discard coordinates, even when read from PDB file
            self.xyz = None


    def version( self ):
        return 'PDBModel $Revision$'


    def __getstate__(self):
        """
        called before pickling the object.
        """
        self.slim()
        self.forcePickle = 0
        return self.__dict__

    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        self.__dict__ = state
        ## backwards compability
        self.__defaults() 

    def __len__(self):
        return self.lenAtoms()


    def __getitem__( self, k ):
        """
        Get atom profile or profile item or CrossView for one atom::
          m['prof1']         <==>  m.aProfiles.get( 'prof1' )         
          m['prof1','info1'] <==>  m.aProfiles.get( 'prof1','info1' )
          m[10]              <==>  CrossView( m.aProfiles, 10 )

        @return: profile OR meta infos thereof OR CrossView dict
        @rtype: list OR array OR any OR CrossView
        """
        if type( k ) is str:
            if k in self.aProfiles:
                return self.aProfiles.get( k )
            if k in self.rProfiles:
                return self.rProfiles.get( k )
            if k in self.info:
                return self.info[ k ]

            return self.profile( k )
        
        if type( k ) is tuple:
            return self.profileInfo( k[0] )[ k[1] ]
        
        return self.aProfiles[k]
    
    
    def __setitem__( self, k, v ):
        """
        Set atom profile or profile item (or meta info)::
          m['prof1'] = range(10)    <==> m.aProfiles.set( 'prof1', range(10) )
            OR                      <==> m.rProfiles.set( 'prof1', range(10) )
          
          m['prof1','info1]='comment'
                             <==> m.aProfiles.setInfo('prof1',info1='comment')
            OR               <==> m.rProfiles.setInfo('prof1',info1='comment')
        
          m['version'] = '1.0.0'    <==> m.info['version'] = '1.0.0'
            but only if 'version' already exists in m.info 

        @return: item
        @rtype: any        
        """
        if type( k ) is str:
            if v is not None and len( v ) == self.lenAtoms():
                return self.aProfiles.set( k, v )
            if v is not None and len( v ) == self.lenResidues():
                return self.rProfiles.set( k, v )
            if k in self.aProfiles:
                return self.aProfiles.set( k, v )
            if k in self.rProfiles:
                return self.rProfiles.set( k, v )
            if k in self.info:
                self.info[ k ] = v
            raise ProfileError, \
                  'Value cannot clearly be assigned to either atom or '+\
                  'residue profiles'
                
        if type( k ) is tuple:
            key, infokey = k
            if key in self.rProfiles:
                return self.rProfiles.setInfo( (k, infokey), v )
            return self.aProfiles.setInfo( (k, infokey), v )
        
        raise ProfileError, 'Cannot interpret %r as profile name or profile info record' % k

    
    def __getslice__( self, *arg ):
        """
        Get list of CrossViews::
          m[0:100:5] <==> [ CrossView(m.aProfiles, i) for i in range(0,100,5) ]
        """
        return self.aProfiles.__getslice__( *arg )
    
    def __iter__( self ):
        return self.aProfiles.iterCrossViews()

    
    def __defaults(self ):
        """
        backwards compatibility to earlier pickled models
        """
        ## convert old profile dictionaries into new ProfileCollections
        if 'resProfiles' in self.__dict__:
            self.rProfiles=PDBProfiles( self,
                                   profiles=getattr(self,'resProfiles',{}),
                                   infos=getattr(self,'resProfiles_info',{}) )
            try:
                del self.resProfiles; del self.resProfiles_info
            except: pass

        if 'atomProfiles' in self.__dict__:
            self.aProfiles=PDBProfiles( self,
                                   profiles=getattr(self,'atomProfiles',{}),
                                   infos=getattr(self,'atomProfiles_info',{}) )
            try:
                del self.atomProfiles; del self.atomProfiles_info
            except: pass

        ## if there are not even old profiles...
        if getattr( self, 'aProfiles', None) is None:
            self.aProfiles = PDBProfiles(self)
        if getattr( self, 'rProfiles', None) is None:
            self.rProfiles = PDBProfiles(self)
            
        self.forcePickle = getattr( self, 'forcePickle', 0 )

        self.info = getattr( self, 'info', { 'date':T.dateSortString() } )

        ## rescue TER positions from previous version
        ter_atoms = getattr( self, '_PDBModel__terAtoms', 0)
        if ter_atoms:
            mask = N.zeros( self.aProfiles.profLength() )
            N.put( mask, ter_atoms, 1 )
            self.aProfiles.set('after_ter', mask,
                               comment='rebuilt from old PDBModel.__terAtoms')

        ## fix previous bug: slim() was creating a self.xyx instead of xyz
        if getattr( self, 'xyx', 0 ):
            del self.xyx

        if type( self.source ) == str:
            self.source = LocalPath( self.source )

        self.__validSource = getattr( self, '_PDBModel__validSource', 0)

        self.initVersion = getattr( self, 'initVersion', 'old PDBModel')

        # biskit 2.0.1 -> biskit 2.1
        if getattr( self.aProfiles, 'model', 0) is 0:
            self.aProfiles.model = self
            self.rProfiles.model = self

        ## rescue atom dictionaries
        atoms = getattr( self, 'atoms', 0)
        if not atoms is 0:
            ## atoms to be fetched from external source
            if atoms is None:
                for k in self.PDB_KEYS:
                    self.aProfiles.set( k, None, changed=0 )
            
            else:
                atoms = B.DictList( atoms )
                for k in self.PDB_KEYS:
                    self.aProfiles.set( k, atoms.valuesOf( k ), 
                                    changed=getattr(self, 'atomsChanged',1) )

            del self.atoms

        self._resIndex = getattr( self, '_resIndex', None)
        self._chainIndex=getattr( self, '_chainIndex', None)

        self.__maskCA = getattr( self, '__maskCA', None )
        self.__maskBB = getattr( self, '__maskBB', None )
        self.__maskHeavy = getattr( self, '__maskHeavy', None )

        try:
            del self.caMask, self.bbMask, self.heavyMask
        except:
            pass


    def update( self, skipRes=None, updateMissing=0, force=0 ):
        """
        Read coordinates, atoms, fileName, etc. from PDB or
        pickled PDBModel - but only if they are currently empty.
        The atomsChanged and xyzChanged flags are not changed.
        
        @param skipRes: names of residues to skip if updating from PDB
        @type  skipRes: list of str
        @param updateMissing: 0(default): update only existing profiles
        @type  updateMissing: 0|1
        @param force: ignore invalid source (0) or report error (1)
        @type  force: 0|1

        @raise PDBError: if file can't be unpickled or read: 
        """
        source = self.validSource()

        if source is None and force:
            raise PDBError( str(self.source) + ' is not a valid source.')

        if source is None:
            return

        parser = PDBParserFactory.getParser( source )
        parser.update(self, source, skipRes=skipRes,
                      updateMissing=updateMissing, force=force )
        
        
    def setXyz(self, xyz ):
        """
        Replace coordinates.

        @param xyz: Numpy array ( 3 x N_atoms ) of float
        @type  xyz: array

        @return: N.array( 3 x N_atoms ) or None, old coordinates
        @rtype: array
        """
        old = self.xyz
        self.xyz = xyz

        self.xyzChanged = self.xyzChanged or \
                          not mathUtils.arrayEqual(self.xyz,old )
        return old


    def setSource( self, source ):
        """
        @param source: LocalPath OR PDBModel OR str
        """
        if type( source ) == str:
            self.source = LocalPath( source )
        else:
            self.source = source
        self.__validSource = 0


    def getXyz( self, mask=None ):
        """
        Get coordinates, fetch from source PDB or pickled PDBModel,
        if necessary.

        @param mask: atom mask
        @type  mask: list of int OR array of 1||0

        @return: xyz-coordinates, N.array( 3 x N_atoms, N.Float32 )
        @rtype: array 
        """
        if self.xyz is None:
            self.update( force=1 )

        if self.xyz is None:
            return N.array( [], N.Float32 )

        if mask is None:
            return self.xyz

        return N.compress( mask, self.xyz, 0 )


    def getAtoms( self, mask=None ):
        """
        Get atom CrossViews that can be used like dictionaries.

        @param mask: atom mask
        @type  mask: list of int OR array of 1||0

        @return: list of CrossView dictionaries
        @rtype: [ L{ProfileCollection.CrossView} ]
        """
        r = self.aProfiles.toCrossViews()

        if mask is None:
            return r

        return [ r[i] for i in N.nonzero( mask ) ]
        
        
    def profile( self, name, default=None, update=True, updateMissing=False ):
        """
        Use::
           profile( name, updateMissing=0) -> atom or residue profile

        @param name: name to access profile
        @type  name: str        
        @param default: default result if no profile is found,
                        if None, raise exception or try updating
        @type  default: any
        @param updateMissing: update from source before reporting missing
                              profile
        @type  updateMissing: 0||1

        @raise ProfileError: if neither atom- nor rProfiles contains |name|
        """
        if updateMissing and not name in self.aProfiles and \
               not name in self.rProfiles:
            self.update( updateMissing=True )

        if name in self.aProfiles:
            return self.aProfiles.get( name, default,
                                       update=update, updateMissing=0)

        if name in self.rProfiles:
            return self.rProfiles.get( name, default,
                                       update=update, updateMissing=0)

        raise ProfileError( 'No profile info found with name '+str(name))
    

    def profileInfo( self, name, updateMissing=0 ):
        """
        Use::
           profileInfo( name ) -> dict with infos about profile

        @param name: name to access profile
        @type  name: str       
        @param updateMissing:update from source before reporting missing profile::
                           Guaranteed infos: 'version'->str,
                                             'comment'->str,
                                             'changed'->1||0
        @type  updateMissing: 0|1

        @raise ProfileError: if neither atom - nor rProfiles
                                   contains |name|
        """
        if updateMissing and not name in self.aProfiles and \
               not name in self.rProfiles:
            self.update()

        if name in self.aProfiles:
            return self.aProfiles.getInfo( name )

        if name in self.rProfiles:
            return self.rProfiles.getInfo( name )

        raise ProfileError( 'No profile info found with name '+str(name))


    def removeProfile( self, *names ):
        """
        Remove residue or atom profile(s)
        
        Use::
           removeProfile( str_name [,name2, name3] ) -> 1|0,

        @param names: name or list of residue or atom profiles
        @type  names: str OR list of str

        @return: 1 if at least 1 profile has been deleted,
                 0 if none has been found
        @rtype: int
        """
        r = 0

        for n in names:
            if n in self.aProfiles:
                del self.aProfiles[n]
                r = 1

            if n in self.rProfiles:
                del self.rProfiles[n]
                r = 1
        return r


    def xyzIsChanged(self):
        """
        Tell if xyz or atoms have been changed compared to source file or
        source object (which can be still in memory).

        @return: xyz field has been changed with respect to source
        @rtype: (1||0, 1||0)
        """
        return self.xyzChanged


    def xyzChangedFromDisc(self):
        """
        Tell whether xyz can currently be reconstructed from a
        source on disc. Same as xyzChanged() unless source is another not yet
        saved PDBModel instance that made changes relative to its own source.

        @return: xyz has been changed
        @rtype: bool
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

        @return: 1, if profile |pname| can currently not be
                 reconstructed from a source on disc.
        @rtype: int

        @raise ProfileError: if there is no atom or res profile with pname
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
        B{AUTOMATICALLY CALLED BEFORE PICKLING}
        """
        for key in self.rProfiles:

            if not self.profileChangedFromDisc( key ):
                self.rProfiles.set( key, None )

            else:

                if type( self.rProfiles[ key ] ) is N.arraytype:
                    self.rProfiles.set( key, self.rProfiles[key].tolist() )

        for key in self.aProfiles:

            if not self.profileChangedFromDisc( key ):
                self.aProfiles.set( key,  None )

            else:
                
                if type( self.aProfiles[ key ] ) is N.arraytype:
                    self.aProfiles.set( key, self.aProfiles[key].tolist() )


    def slim( self ):
        """
        Remove xyz array and profiles if they haven't been changed and
        could hence be loaded from the source file (only if there is a source
        file...).
        B{AUTOMATICALLY CALLED BEFORE PICKLING}
        """
        ## remove atoms/coordinates if they are unchanged from an existing
        ## source
        ## override this behaviour with forcePickle
        if not self.forcePickle:

            if not self.xyzChangedFromDisc():
                self.xyz = None

            if type( self.xyz ) is N.arraytype and self.xyz.dtype.char != 'f':
                self.xyz = self.xyz.astype(N.Float32)

            self.__slimProfiles()

        self.__maskCA = self.__maskBB = self.__maskHeavy = None
        self.__validSource = 0


    def validSource(self):
        """
        Check for a valid source on disk.

        @return:  str or PDBModel, None if this model has no valid source
        @rtype: str or PDBModel or None
        """
        if self.__validSource == 0:

            if isinstance( self.source, LocalPath ) and self.source.exists():
                self.__validSource = self.source.local()
            else:
                if isinstance( self.source, B.PDBModel ):
                    self.__validSource = self.source
                else:
                    self.__validSource = None

        return self.__validSource


    def sourceFile( self ):
        """
        Name of pickled source or PDB file. If this model has another
        PDBModel as source, the request is passed on to this one.

        @return: file name of pickled source or PDB file
        @rtype: str

        @raise PDBError: if there is no valid source
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

        @note: If this model has an (in-memory) PDBModel instance as source,
        the entries of 'atoms' could still reference the same dictionaries.
        """
        self.update()

        try:
            self.fileName = self.fileName or self.sourceFile()
        except:
            pass

        self.setSource( None )

        self.xyzChanged = 1

        for p in self.rProfiles:
            self.rProfiles.setInfo( p, changed=1 )
        for p in self.aProfiles:
            self.aProfiles.setInfo( p, changed=1 )


    def getPdbCode(self):
        """
        Return pdb code of model.

        @return: pdb code
        @rtype: str
        """
        return self.pdbCode

    def setPdbCode(self, code ):
        """
        Set model pdb code.

        @param code: new pdb code
        @type  code: str
        """
        self.pdbCode = code


    def sequence(self, mask=None, xtable=molUtils.xxDic ):
        """
        Amino acid sequence in one letter code.

        @param mask: atom mask, to apply before  (default None)
        @type  mask: list or array
        @param xtable: dict {str:str}, additional residue:single_letter mapping
                       for non-standard residues (default molUtils.xxDic)
                       [currently not used]
        @type  xtable: dict

        @return: 1-letter-code AA sequence (based on first atom of each res).
        @rtype: str
        """
        firstAtm = self.resIndex()
        if mask is not None:
            m_first = N.zeros( self.lenAtoms() )
            N.put( m_first, firstAtm, 1 )
            m_first = mask * m_first
            firstAtm = N.nonzero( m_first )

        l = self.aProfiles['residue_name']
        l = [ l[i] for i in firstAtm ]

        return ''.join( molUtils.singleAA( l, xtable ) )


    def xplor2amber( self, aatm=1 ):
        """
        Rename atoms so that tleap from Amber can read the PDB.
        If HIS residues contain atoms named HE2 or/and HD2, the residue
        name is changed to HIE or HID or HIP, respectively. Disulfide bonds
        are not yet identified - CYS -> CYX renaming must be done manually
        (see AmberParmBuilder for an example). 
        Internally amber uses H atom names ala HD21 while standard pdb files
        use 1HD2. By default, ambpdb produces 'standard' pdb atom names but
        it gives the less ambiguous amber names with switch -aatm.

        @param change: change this model's atoms directly (default:1)
        @type  change: 1|0
        @param aatm: use, for example, HG23 instead of 3HG2 (default:1)
        @type  aatm: 1|0

        @return: [ {..} ], list of atom dictionaries
        @rtype: list of atom dictionaries
        """
        numbers = map( str, range(10) )

        resI = self.resIndex()

        names    = self.aProfiles['name']
        resnames = self.aProfiles['residue_name']

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
        l = self.aProfiles['residue_name']

        for i in range(len(l)):
            if l[i] == 'CYX':
                l[i] = 'CYS'
            if l[i] in  ['HIE','HID','HIP']:
                l[i] = 'HIS'
            
    def writePdb( self, fname, ter=1, amber=0, original=0, left=0, wrap=0,
                  headlines=None, taillines=None):
        """
        Save model as PDB file.

        @param fname: name of new file
        @type  fname: str
        @param ter: Option of how to treat the terminal record::
                    0, don't write any TER statements
                    1, restore original TER statements (doesn't work,
                         if preceeding atom has been deleted) [default]
                    2, put TER between all detected chains
                    3, as 2 but also detect and split discontinuous chains
        @type  ter: 0, 1, 2 or 3
        @param amber: amber formatted atom names
                      (implies ter=3, left=1, wrap=0) (default 0)
        @type  amber: 1||0
        @param original: revert atom names to the ones parsed in from PDB
                         (default 0)
        @type  original: 1||0
        @param left: left-align atom names (as in amber pdbs)(default 0)
        @type  left: 1||0
        @param wrap: write e.g. 'NH12' as '2NH1' (default 0)
        @type  wrap: 1||0
        @param headlines: [( str, dict or str)], list of record / data tuples::
                          e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
        @type  headlines: list of tuples 
        @param taillines: same as headlines but appended at the end of file
        @type  taillines: list of tuples 
        """
        try:
            f = IO.PDBFile( fname, mode='w' )

            numbers = map( str, range(10) )

            if amber:
                __resnames = copy.copy(self.aProfiles['residue_name'])
                __anames   = copy.copy(self.aProfiles['name'])
                self.xplor2amber()
                ter = 3
                wrap = 0
                left = 1

            if ter == 2 or ter == 3:
                ## tolist to workaround numeric bug
                terIndex = N.array( self.chainIndex( breaks=(ter==3) ).tolist()[1:] )
            if ter == 1:
                terIndex = N.nonzero( self.aProfiles['after_ter'] )

            if headlines:
                for l in headlines:
                    f.writeLine( l[0], l[1] )

            if taillines:
                for l in taillines:
                    f.writeLine( l[0], l[1] )

            i = -1
            for a in self.aProfiles.toDicts():
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
                self.aProfiles['residue_name'] = __resnames
                self.aProfiles['name'] = __anames

        except:
            EHandler.error( "Error writing "+fname )


##     def returnPdb( self, out=None, ter=1, headlines=None, taillines=None):
##         """
##         Restore PDB file from (possibly transformed) coordinates and pdb line
##         dictionaries in self.atoms. This is an older version of writePdb that
##         returns a list of PDB lines instead of writing to a file.

##         @param out: stdout or None, if None a list is returned
##         @type  out: stdout or None
##         @param ter: Option of how to treat the retminal record::
##                     0, don't write any TER statements
##                     1, restore original TER statements (doesn't work,
##                          if preceeding atom has been deleted)
##                     2, put TER between all detected chains
##         @type  ter: 0, 1 or 2
##         @param headlines: [( str, dict or str)], list of record / data tuples::
##                           e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
##         @type  headlines: list of tuples
##         @param taillines: same as headlines but appended at the end of file
##         @type  taillines: list of tuples

##         @return: [ str ], lines of a PDB file
##         @rtype: list of strings
##         """
##         try:
##             pdb_lst=[]
##             i=0

##             if ter == 2:
##                 terIndex = N.array( self.chainIndex()[1:] ) - 1
##             if ter == 1:
##                 terIndex = self.__terAtoms

##             if headlines:
##                 for l in headlines:
##                     pdb_lst += [ '%6s%s\n'%(l[0],l[1]) ]

##             if taillines:
##                 for l in taillines:
##                     pdb_lst += [ '%6s%s\n'%(l[0],l[1]) ]

##             for a in self.getAtoms():

##                 ## fetch coordinates Vector
##                 a['positionX'], a['positionY'], a['positionZ'] = self.xyz[ i ]

##                 ## switch original and short atom name
##                 a['name'], a['name_original'] = a['name_original'], a['name']

##                 ## write line
##                 atom_line = 'ATOM  %(serial_number)5i %(name)-4s %(residue_name)3s %(chain_id)1s%(residue_number)4i%(insertion_code)1s   %(positionX)8.3f%(positionY)8.3f%(positionZ)8.3f%(occupancy)6.2f%(temperature_factor)6.2f      %(segment_id)-4s%(element)2s\n'
##                 pdb_lst += [ atom_line%a ]

##                 ## switch back
##                 a['name'], a['name_original'] = a['name_original'], a['name']
##                 del( a['positionX'], a['positionY'], a['positionZ'] )

##                 ## write TER line with details from previous atom
##                 if (ter and  i in terIndex):
##                     ter_line = 'TER   %(serial_number)5i      %(residue_name)3s %(chain_id)1s%(residue_number)4i%(insertion_code)1s\n'
##                     pdb_lst += [ ter_line % self.getAtoms()[i] ]

##                 i += 1
##             ## write to stdout or return list 
##             if out == 'stdout':
##                 sys.stdout.writelines( pdb_lst )
##             else:
##                 return pdb_lst

##         except:
##             EHandler.error( "PDBModel.returnPdb(): " )


    def saveAs(self, path):
        """
        Pickle this PDBModel to a file, set the 'source' field to
        this file name and mark atoms, xyz, and profiles  as unchanged.
        Normal pickling of the object will only dump those data that can not
        be reconstructed from the source of this model (if any).
        saveAs creates a 'new source' without further dependencies.

        @param path: target file name
        @type  path: str OR LocalPath instance
        """
        try:
            self.update()

            ## pickle all atoms, coordinates, profiles, even if unchanged
            self.forcePickle = 1

            self.setSource( path )

            self.xyzChanged = 0

            for p in self.rProfiles:
                self.rProfiles.setInfo( p, changed=0 )
            for p in self.aProfiles:
                self.aProfiles.setInfo( p, changed=0 )

            T.Dump( self, str(path) )

        except IOError, err:
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

        use:

          r = m.maskFrom( 'name', 'CA' ) * m.maskFrom('residue_name', 'ALA')

        @param atomFunction: function( dict_from_aProfiles.toDict() ),
                             true || false (Condition)
        @type  atomFunction: 1||0
        @param numpy: 1(default)||0, convert result to Numpy array of int
        @type  numpy: int

        @return: Numpy N.array( [0,1,1,0,0,0,1,0,..], N.Int) or list
        @rtype: array or list
        """
        try:
            result = map( atomFunction, self.aProfiles.toDicts() )
        except:

            ## fall-back solution: assign 0 to all entries that raise
            ## exception
            EHandler.warning("mask(): Error while mapping funtion "+
                       "to all atoms.")
            result = []

            for a in self.aProfiles.iterDicts():
                try:
                   result.append( atomFunction( a ) )
                ## put 0 if something goes wrong
                except :
                    EHandler.warning("mask(): Error while save-mapping ")
                    result.append(0)

        if numpy:
            return N.array( result )
        return result


    def maskFrom( self, key, cond ):
        """
        Create an atom mask from the values of a specific profile.
        Example, the following three statements are equivalent::

          >>> mask = m.maskFrom( 'name', 'CA' )
          >>> mask = m.maskFrom( 'name', lambda a: a == 'CA' )
          >>> mask = N.array( [ a == 'CA' for a in m.aProfiles['name'] ] )

        People having numpy installed can also simply use::
          >>> mask = numpy.array(m.aProfiles['name']) == 'CA'

        @param key: the name of the profile to use
        @type  key: str
        @param cond: either a function accepting a single value or a value or
                     an iterable of values (to allow several alternatives)
        @type  cond: function OR any OR [ any ]
        @return: array or list of indices where condition is met
        @rtype : list or N.array of int
        """

        if type( cond ) is types.FunctionType:
            return N.array( map( cond, self.aProfiles[ key ] ) )

        ## several allowed values given
        elif type( cond ) in [ list, tuple ]:
            return N.array( [ x in cond for x in self.aProfiles[key] ] )

        ## one allowed value given
        ## Numeric splits lists of str into 2-D char arrays, 'O' prevents that
        else:
            return N.array( self.aProfiles[key] ) == cond


    def maskCA( self, force=0 ):
        """
        Short cut for mask of all CA atoms.

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.__maskCA == None or force:
            self.__maskCA = self.maskFrom( 'name', 'CA' )

        return self.__maskCA


    def maskBB( self, force=0 ):
        """
        Short cut for mask of all backbone atoms.

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.__maskBB == None or force:
            self.__maskBB = self.maskFrom( 'name',
                                         ['CA', 'C', 'N', 'O', 'H','OXT'] )

        return self.__maskBB


    def maskHeavy( self, force=0 ):
        """
        Short cut for mask of all heavy atoms. ('element' <> H)

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.__maskHeavy == None or force:
            self.__maskHeavy = self.maskFrom( 'element', lambda a: a != 'H' )

        return self.__maskHeavy

    def maskH( self ):
        """
        Short cut for mask of hydrogens. ('element' == H)

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return N.logical_not( self.maskHeavy() )


    def _maskCB( self ):
        """
        Short cut for mask of all CB I{and} CA of GLY.

        @return: mask of all CB and CA of GLY
        @rtype: array
        """
        f = lambda a: a['name'] == 'CB' or\
            a['residue_name'] == 'GLY' and a['name'] == 'CA' 

        return self.maskF( f )


    def _maskCB( self ):
        """
        Short cut for mask of all CB I{and} CA of GLY.

        @return: mask of all CB plus CA of GLY
        @rtype: array
        """
        m_cb = self.maskFrom( 'name', 'CB' )
        m_g  = self.maskFrom( 'residue_name', 'GLY' )
        m_ca = self.maskCA()

        return m_cb + (m_g * m_ca)


    def maskH2O( self ):
        """
        Short cut for mask of all atoms in residues named TIP3, HOH and  WAT

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskFrom( 'residue_name', ['TIP3','HOH','WAT'] )

    def maskSolvent( self ):
        """
        Short cut for mask of all atoms in residues named
        TIP3, HOH, WAT, Na+, Cl-

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskFrom('residue_name', ['TIP3','HOH','WAT','Na+', 'Cl-'])
    
    def maskHetatm( self ):
        """
        Short cut for mask of all HETATM 

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskFrom( 'type', 'HETATM' )

    def maskProtein( self, standard=0 ):
        """
        Short cut for mask containing all atoms of amino acids.

        @param standard: only standard residue names (not CYX, NME,..)
                         (default 0)
        @type  standard: 0|1

        @return: N.array( 1 x N_atoms ) of 0||1,
                 mask of all protein atoms (based on residue name)
        @rtype: array
        """
        d = molUtils.aaDic
        if standard:
            d = molUtils.aaDicStandard

        names = map( string.upper, d.keys() )
        return N.array(
            [ n.upper() in names for n in self.aProfiles['residue_name'] ] )


    def indicesFrom( self, key, cond ):
        """
        Get atom indices conforming condition applied to an atom profile.
        Corresponds to::
          >>> Numeric.nonzero( m.maskFrom( key, cond) )

        @param key: the name of the profile to use
        @type  key: str
        @param cond: either a function accepting a single value or a value or
                     an iterable of values
        @type  cond: function OR any OR [any]
        @return: array of indices where condition is met
        @rtype : N.array of int
        """
        return N.nonzero( self.maskFrom( key, cond) )


    def indices( self, what ):
        """
        Get atom indices conforming condition. This is a convenience method
        to 'normalize' different kind of selections (masks, atom names,
        indices, functions) to indices as they are e.g. required by
        L{PDBModel.take}. 

        @param what: Selection::
             - function applied to each atom entry,
                e.g. lambda a: a['residue_name']=='GLY'
             - list of str, allowed atom names
             - list of int, allowed atom indices OR mask with only 1 and 0
             - int, single allowed atom index
        @type  what: function OR list of str or int OR int

        @return: N_atoms x 1 (0||1 )
        @rtype: Numeric array

        @raise PDBError: if what is neither of above
        """
        ## lambda funcion
        if type( what ) is types.FunctionType:
            return N.nonzero( self.maskF( what) )

        if type( what ) is list or type( what ) is N.arraytype:

            ## atom names
            if type( what[0] ) == str:
                return self.indicesFrom( 'name', what )

            if isinstance( what[0] , int) or \
               (isinstance(what, N.ndarray) and what.dtype in [int, bool]):
                ## mask
                if len( what ) == self.lenAtoms() and max( what ) < 2:
                    return N.nonzero( what )
                ## list of indices
                else:
                    return what

        ## single index
        if isinstance( what , int):
            return N.array( [what], N.Int )

        raise PDBError("PDBModel.indices(): Could not interpret condition ")


    def mask( self, what ):
        """
        Get atom mask. This is a convenience method to 'normalize'
        different kind of selections (masks, atom names, indices,
        functions) to a mask as it is e.g. required by L{PDBModel.compress}.

        @param what: Selection::
                     - function applied to each atom entry,
                        e.g. lambda a: a['residue_name']=='GLY'
                     - list of str, allowed atom names
                     - list of int, allowed atom indices OR mask with
                       only 1 and 0
                     - int, single allowed atom index
        @type  what: function OR list of str or int OR int

        @return: N_atoms x 1 (0||1 )
        @rtype: Numeric array

        @raise PDBError: if what is neither of above
        """
        ## lambda funcion
        if type( what ) == types.FunctionType:
            return self.maskF( what )

        if type( what ) == list or type( what ) is N.arraytype:

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
                    r = N.zeros( self.lenAtoms(),N.Int )
                    N.put( r, what, 1 )
                    return r

        ## single index
        if isinstance( what , int):
            return self.mask( [what] )

        raise PDBError, "PDBModel.mask(): Could not interpret condition "


    def index2map( self, index, len_i ):
        """
        Create a map of len_i length, giving the residue(/chain) numer of
        each atom, from list of residue(/chain) starting positions.

        @param index: list of starting positions, e.g. [0, 3, 8]
        @type  index: [ int ] or N.array of int
        @param len_i: length of target map, e.g. 10
        @type  len_i: int

        @return: list mapping atom positions to residue(/chain) number,
                 e.g. [0,0,0, 1,1,1,1,1, 2,2] from above example
        @rtype: N.array of int (and of len_i length)
        """
        index = N.concatenate( (index, [len_i]) )
        delta = index[1:] - index[:-1] 
        ## Numeric: delta = N.take( index, range(1, len(index) ) ) - index[:-1]
        return N.repeat( range(len(delta)), delta)


    def map2index( self, imap ):
        """
        Identify the starting positions of each residue(/chain) from a map
        giving the residue(/chain) number of each atom.

        @param imap: something like [0,0,0,1,1,1,1,1,2,2,2,...]
        @type  imap: [ int ]

        @return: list of starting positions, e.g. [0, 3, 8, ...] in above ex.
        @rtype: N.array of int
        """
        try:
            imap = N.concatenate( (imap, [imap[-1]] ) )
            delta = imap[1:] - imap[:-1] 
            # Numeric: delta = N.take( imap, range(1, len(imap) ) ) - imap[:-1]
            r = N.nonzero( delta ) + 1
            return N.concatenate( ( [0], r ) )

        except IndexError:
            ## handle empty imap parameter
            return N.zeros(0)


    def extendMask( self, mask, index, len_i ):
        """
        Translate a mask that is defined,e.g., on residues(/chains) to a mask
        that is defined on atoms.

        @param mask : mask marking positions in the list of residues or chains
        @type  mask : [ bool ] or N.array of bool or of 1||0
        @param index: starting positions of all residues or chains
        @type  index: [ int ] or N.array of int
        @param len_i: length of target mask
        @type  len_i: int

        @return: mask that blows up the residue / chain mask to an atom mask
        @rtype: N.array of bool
        """
        index = N.concatenate( (index, [len_i]) )
        delta = index[1:] - index[:-1] 
        # Numeric: delta = N.take( index, range(1, len(index) ) ) - index[:-1]
        return N.repeat( mask, delta )


    def atom2resMask( self, atomMask ):
        """
        Mask (set 0) residues for which all atoms are masked (0) in atomMask.

        @param atomMask: list/array of int, 1 x N_atoms
        @type  atomMask: list/array of int

        @return: 1 x N_residues (0||1 )
        @rtype: array of int
        """
        res_indices = self.atom2resIndices( N.nonzero( atomMask) )
        r = N.zeros( self.lenResidues() )
        N.put( r, res_indices, 1 )
        return r
        

    def atom2resIndices( self, indices ):
        """
        Get list of indices of residues for which any atom is in indices.

        @param indices: list of atom indices
        @type  indices: list of int

        @return: indices of residues
        @rtype: list of int
        """
        new_resmap = N.take( self.resMap(), indices )
        resIndex   = self.map2index( new_resmap )

        return N.take( new_resmap, resIndex )


    def res2atomMask( self, resMask ):
        """
        Convert residue mask to atom mask.

        @param resMask: list/array of int, 1 x N_residues
        @type  resMask: list/array of int

        @return: 1 x N_atoms
        @rtype: array of int
        """
        return self.extendMask( resMask, self.resIndex(), self.lenAtoms() )


    def res2atomIndices( self, indices ):
        """
        Convert residue indices to atom indices.

        @param indices: list/array of residue indices
        @type  indices: list/array of int

        @return: list of atom indices
        @rtype: list of int
        """
        if max( indices ) > self.lenResidues() or min( indices ) < 0:
            raise PDBError, "invalid residue indices"

        resMask = N.zeros( self.lenResidues() )
        N.put( resMask, indices, 1 )

        atomMask= self.res2atomMask( resMask )

        return N.nonzero( atomMask )


    def atom2chainIndices( self, indices, breaks=0 ):
        """
        Convert atom indices to chain indices. Each chain is only
        returned once.

        @param indices: list of atom indices
        @type  indices: list of int
        @param breaks: look for chain breaks in backbone coordinates (def. 0)
        @type  breaks: 0||1

        @return: chains any atom which is in indices
        @rtype: list of int
        """
        new_map = N.take( self.chainMap( breaks=breaks ), indices )
        index   = self.map2index( new_map )

        return N.take( new_map, index )


    def atom2chainMask( self, atomMask, breaks=0 ):
        """
        Mask (set to 0) chains for which all atoms are masked (0) in atomMask.
        Put another way: Mark all chains that contain any atom that is marked
        '1' in atomMask.

        @param atomMask: list/array of int, 1 x N_atoms
        @type  atomMask: list/array of int

        @return: 1 x N_residues (0||1 )
        @rtype: array of int
        """
        indices = self.atom2chainIndices( N.nonzero( atomMask), breaks=breaks )
        r = N.zeros( self.lenChains(breaks=breaks) )
        N.put( r, indices, 1 )
        return r


    def chain2atomMask( self, chainMask, breaks=0 ):
        """
        Convert chain mask to atom mask.

        @param chainMask: list/array of int, 1 x N_chains
        @type  chainMask: list/array of int
        @param breaks: look for chain breaks in backbone coordinates (def. 0)
        @type  breaks: 0||1

        @return: 1 x N_atoms
        @rtype: array of int
        """
        return self.extendMask( chainMask, self.chainIndex( breaks=breaks ),
                                self.lenAtoms() )


    def chain2atomIndices( self, indices, breaks=0 ):
        """
        Convert chain indices into atom indices.

        @param indices: list of chain indices
        @type  indices: list of int
        @param breaks: look for chain breaks in backbone coordinates (def. 0)
        @type  breaks: 0||1

        @return: all atoms belonging to the given chains
        @rtype: list of int
        """
        if max( N.absolute(indices) ) > self.lenChains( breaks=breaks ):
            raise PDBError, "invalid chain indices"

        chainMask = N.zeros( self.lenChains( breaks=breaks ) )
        N.put( chainMask, indices, 1 )

        atomMask= self.chain2atomMask( chainMask, breaks=breaks )

        return N.nonzero( atomMask )


    def res2atomProfile( self, p ):
        """
        Get an atom profile where each atom has the value its residue has
        in the residue profile.

        @param p: name of existing residue profile OR ...
                  [ any ], list of lenResidues() length
        @type  p: str

        @return: [ any ] OR array, atom profile
        @rtype: list or array
        """
        if type( p ) is str:
            p = self.rProfiles.get( p )

        isArray = isinstance( p, N.arraytype )

        resMap = self.resMap()

        r = [ p[ resMap[a] ] for a in range( len(resMap) ) ]

        if isArray:
            r = N.array( r )

        return r


    def atom2resProfile( self, p ):
        """
        Get a residue profile where each residue has the value that its first
        atom has in the atom profile.
        @param p: name of existing atom profile OR ...
                  [ any ], list of lenAtoms() length
        @type  p: str

        @return: [ any ] OR array, residue profile
        @rtype: list or array
        """
        if type( p ) is str:
            p = self.aProfiles.get( p )

        isArray = isinstance( p, N.arraytype )
        
        r = N.take( p, self.resIndex() )

        if not isArray:
            return r.tolist()
        return r
    

    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2mask( str_profname, [cutoff_min, cutoff_max=None])

        @param cutoff_min: low value cutoff
        @type  cutoff_min: float
        @param cutoff_max: high value cutoff
        @type  cutoff_max: float

        @return: mask len( profile(profName) ) x 1||0
        @rtype: array

        @raise ProfileError: if no profile is found with name profName
        """
        p = self.profile( profName )

        cutoff_min = cutoff_min or min( p ) - 1
        cutoff_max = cutoff_max or max( p ) + 1

        return N.greater( p, cutoff_min ) * N.less( p, cutoff_max )


    def profile2atomMask( self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2atomMask( str_profname, [cutoff_min, cutoff_max=None])
        Same as L{profile2mask}, but converts residue mask to atom mask.

        @param cutoff_min: low value cutoff
        @type  cutoff_min: float
        @param cutoff_max: high value cutoff
        @type  cutoff_max: float

        @return: mask N_atoms x 1|0
        @rtype: array

        @raise ProfileError: if no profile is found with name profName
        """
        r = self.profile2mask( profName, cutoff_min, cutoff_max )

        if len( r ) == self.lenResidues():
            r = self.res2atomMask( r )

        return r


    def concat( self, *models ):
        """
        Concatenate atoms, coordinates and profiles. source and fileName
        are lost, so are profiles that are not available in all models.
        model0.concat( model1 [, model2, ..]) -> single PDBModel.

        @param models: models to concatenate
        @type  models: model OR list of models

        @note: info records of given models are lost.
        """

        if len( models ) == 0:
            return self

        m = models[0]

        r = self.__class__()

        self.update()  ## trigger update if xyz or any profile is None
        
        r.setXyz( N.concatenate( ( self.getXyz(), m.getXyz() )  ) )

        r.setPdbCode( self.pdbCode )

        r.rProfiles = self.rProfiles.concat( m.rProfiles, )
        r.aProfiles = self.aProfiles.concat( m.aProfiles, )

        r.rProfiles.model = r
        r.aProfiles.model = r

        r._resIndex   = N.concatenate(
            (self.resIndex(), m.resIndex() + self.lenResidues())) 
        r._chainIndex =N.concatenate(
            (self.chainIndex(), m.chainIndex() +self.lenChains()))

        r.info = copy.deepcopy( self.info )

        return r.concat( *models[1:] )


    def take( self, i, *initArgs, **initKw ):
        """
        Extract a PDBModel with a subset of atoms::
          take( atomIndices ) -> PDBModel / sub-class.

        @param i: atomIndices, positions to take in the order to take
        @type  i: list/array of int

        @return: PDBModel / sub-class
        @rtype: PDBModel
        """
        r = self.__class__( *initArgs, **initKw )

        ## the easy part: extract coordinates and atoms
        r.xyz = N.take( self.getXyz(), i )
        r.xyzChanged = self.xyzChanged or not N.all(r.xyz == self.xyz)

        r.aProfiles = self.aProfiles.take( i, r )

        ## more tricky: rescue residue borders and extract residue profiles
        new_resmap   = N.take( self.resMap(), i )
        r._resIndex = self.map2index( new_resmap )

        i_res     = N.take( new_resmap, r._resIndex )
        r.rProfiles = self.rProfiles.take( i_res, r )

        ## now the same with chain borders (and later profiles)
        new_chainmap   = N.take( self.chainMap(), i )
        r._chainIndex = self.map2index( new_chainmap )

        ## copy non-sequential infos
        r.info = copy.deepcopy( self.info )
        r.pdbCode = self.pdbCode
        r.fileName = self.fileName
        r.source = self.source

        return r
        


    def keep( self, i ):
        """
        Replace atoms,coordinates,profiles of this(!) model with sub-set.
        (in-place version of N.take() )

        @param i: atom positions to be kept
        @type  i: list or array of int
        """
        if len(i)==self.lenAtoms() and max(i)<2:
            EHandler.warning('dont use PDBModel.keep() with mask.', trace=0) 

        r = self.take( i )

        self.xyz = r.xyz
        self.xyzChanged = r.xyzChanged

        self._resIndex  = r._resIndex
        self._chainIndex= r._chainIndex

        self.aProfiles = r.aProfiles
        self.rProfiles = r.rProfiles

        self.info = r.info

        self.__maskCA = self.__maskBB = self.__maskHeavy = None


    def clone( self ):
        """
        Clone PDBModel.

        @return: PDBModel / subclass, copy of this model,
                 see comments to Numeric.take()
        @rtype: PDBModel
        """
        return self.take( self.atomRange() )


    def compress( self, mask, *initArgs, **initKw ):
        """
        Compress PDBmodel using mask.
        compress( mask ) -> PDBModel

        @param mask: N.array( 1 x N_atoms of 1 or 0 )
                     1 .. keep this atom
        @type  mask: array

        @return: compressed PDBModel using mask
        @rtype: PDBModel
        """
        return self.take( N.nonzero( mask ), *initArgs, **initKw )


    def remove( self, what ):
        """        
        Convenience access to the 3 different remove methods.
        The mask used to remove atoms is returned. This mask can be used
        to apply the same change to another array of same dimension as
        the old(!) xyz and atoms.

        @param what: Decription of what to remove::
              - function( atom_dict ) -> 1 || 0    (1..remove) OR
              - list of int [4, 5, 6, 200, 201..], indices of atoms to remove
              - list of int [11111100001101011100..N_atoms], mask (1..remove)
              - int, remove atom with this index
        @type  what: list of int or int

        @return: N.array(1 x N_atoms_old) of 0||1, mask used to compress the
                 atoms and xyz arrays. 
        @rtype: array

        @raise PDBError: if what is neither of above
        """
        mask = N.logical_not( self.mask( what ) )
        self.keep( N.nonzero(mask) )
        return mask


    def takeChains( self, chains, breaks=0 ):
        """
        Get copy of this model with only the given chains.

        @param chains: list of chains, e.g. [0,2] for first and third
        @type  chains: list of int
        @param breaks: split chains at chain breaks (default 0)
        @type  breaks: 0|1
        @param maxDist: (if breaks=1) chain break threshold in Angstrom
        @type  maxDist: float

        @return: PDBModel consisting of the given chains in the given order
        @rtype : PDBModel / subclass
        """
        return self.take( self.chain2atomIndices( chains, breaks=breaks ) )
    

    def addChainFromSegid(self, verbose=1):
        """
        Takes the last letter of the segment ID and adds it as chain ID.
        """
        chain_ids   = self.aProfiles['chain_id']
        segment_ids = self.aProfiles['segment_id']

        for i in self.atomRange():

            try:
                chain_ids[i] = segment_ids[i][-1]
            except:
                if verbose:
                    EHandler.warning("addChainId(): Problem with atom "+str(a))


    def addChainId( self, first_id=None, keep_old=0, breaks=0 ):
        """
        Assign consecutive chain identifiers A - Z to all atoms.

        @param first_id: str (A - Z), first letter instead of 'A'
        @type  first_id: str 
        @param keep_old: don't override existing chain IDs (default 0)
        @type  keep_old: 1|0
        @param breaks: consider chain break as start of new chain (default 0)
        @type  breaks: 1|0
        """
        ids = self.aProfiles['chain_id']

        old_chains = []
        if keep_old:
            old_chains = N.take( ids, self.resIndex() )
            old_chains = mathUtils.nonredundant( old_chains )
            if '' in old_chains: old_chains.remove('')

        letters = string.uppercase
        if first_id:
            letters = letters[ letters.index( first_id ): ]
        letters = mathUtils.difference( letters, old_chains )

        chainMap = self.chainMap( breaks=breaks )

        try:
            for i in self.atomRange():

                if not (keep_old and ids[i] in old_chains):
                    ids[i] = letters[ chainMap[i] ]

        except IndexError:
            raise PDBError, 'Too many chains, running out of letters.'


    def renumberResidues( self, mask=None, start=1, addChainId=1 ):
        """
        Make all residue numbers consecutive and remove any insertion code
        letters. Note that a backward jump in residue numbering is interpreted
        as end of chain by chainMap() and chainIndex(). Chain borders might
        hence get lost if there is no change in chain label or segid.

        @param mask: [ 0||1 x N_atoms ] atom mask to apply BEFORE
        @type  mask: list of int
        @param start: starting number (default 1)
        @type  start: int
        @param addChainId: add chain IDs if they are missing
        @type  addChainId: 1|0
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
        @return: integer range for lenght of this model
        @rtype: [ int ]
        """
        return range( self.lenAtoms() )


    def lenAtoms( self ):
        """
        Number of atoms in model.

        @return: number of atoms
        @rtype: int
        """
        if not self.xyz is None:
            return len( self.xyz )

        if len( self.aProfiles ) > 0:
            r = self.aProfiles.profLength( default=-1 )
            if r != -1:
                return r

        if self.source is None:
            return 0

        return len( self.getXyz() )


    def lenResidues( self ):
        """
        Number of resudies in model.

        @return: total number of residues
        @rtype: int
        """
        if self._resIndex is None:
            return 0

        return len( self.resIndex() )


    def lenChains( self, breaks=0, maxDist=None ):
        """
        Number of chains in model.

        @param breaks: detect chain breaks from backbone atom distances (def 0)
        @type  breaks: 0||1
        @param maxDist: maximal distance between consequtive residues
                        [ None ] .. defaults to twice the average distance
        @type  maxDist: float

        @return: total number of chains
        @rtype: int
        """
        try:
            return len( self.chainIndex( breaks=breaks, maxDist=maxDist ) )
        except IndexError:  ## empty residue map
            return 0


    def resList( self, mask=None ):
        """
        Return list of lists of atom pseudo dictionaries per residue,
        which allows to iterate over residues and atoms of residues.

        @param mask: [ 0||1 x N_atoms ] atom mask to apply BEFORE
        @type  mask: 

        @return: A list of dictionaries::
        [ [ CrossView{'name':'N', ' residue_name':'LEU', ..},          
            CrossView{'name':'CA', 'residue_name':'LEU', ..} ],        
          [ CrossView{'name':'CA', 'residue_name':'GLY', ..}, .. ] ]      
        @rtype: list of list of CrossView    
        """
        ri = N.concatenate( (self.resIndex( mask=mask ), [self.lenAtoms()] ) )
        resLen = len( ri ) - 1
        atoms = self.getAtoms()
        if mask is not None:
            atoms = N.compress( mask, atoms ).tolist()

        return [ atoms[ ri[res] : ri[res+1] ] for res in range( resLen ) ] 


    def resModels( self ):
        """
        Creates a PDBModel per residue in PDBModel.

        @return: list of PDBModels, one for each residue
        @rtype: list of PDBModels
        """
        ri = self.resIndex()
        resLen = len( ri )
        xyz = self.getXyz()

        result = []
        for res in  range( resLen ):

            a = ri[res]
            if res == resLen - 1:
                e = self.lenAtoms()
            else:
                e = ri[res + 1]

            result += [ self.take( range(a,e) ) ]

        return result


    def resMapOriginal(self, mask=None):
        """
        Generate list to map from any atom to its ORIGINAL(!) PDB
        residue number.

        @param mask: [00111101011100111...] consider atom: yes or no
                     len(mask) == N_atoms
        @type  mask: list of int (1||0)

        @return: list all [000111111333344444..] with residue number
                 for each atom
        @rtype: list of int
        """
        ## by default take all atoms
        if mask is None: mask = N.ones( self.lenAtoms() , N.Int )

        ## apply mask to this list
        return N.compress( mask, self.aProfiles['residue_number'] )


    def __inferResIndex( self ):
        """
        Determine residue borders.

        @return: starting position of each residue
        @rtype:  list of int
        """
        result = []

        lastResNumber = -100
        lastResName   = ''
        index = -1
        lastAlt = 'x'

        res_nrs = self.aProfiles['residue_number']
        res_nam = self.aProfiles['residue_name']
        ins_cod = self.aProfiles['insertion_code']

        ## create residue numbering for selected atoms
        for i in range( self.lenAtoms() ):

            if res_nrs[i] != lastResNumber or \
               res_nam[i] != lastResName   or \
               ins_cod[i] != lastAlt:

                ## start of new residue
                lastResNumber = res_nrs[i]
                lastResName   = res_nam[i]
                lastAlt       = ins_cod[i]
                index += 1

                result.append( i )

        return N.array(result, N.Int)        


    def resIndex( self, mask=None, force=0, cache=1 ):
        """
        Get the position of the each residue's first atom.

        @param force: re-calculate even if cached result is available (def 0)
        @type  force: 1||0
        @param cache: cache the result if new (def 1)
        @type  cache: 1||0
        @param mask: atom mask to apply before (i.e. result indices refer to
                     compressed model)
        @type  mask: list of int (1||0)

        @return: index of the first atom of each residue
        @rtype: list of int
        """

        if self._resIndex is not None and not force and mask is None:
            return self._resIndex

        r = self.__inferResIndex()

        if mask is not None:
            m = self.index2map( r, len( mask ) )
            m = N.compress( mask, m )
            r = self.map2index( m )

        if mask is None and cache:
            self._resIndex = r

        return r


    def resMap(  self, force=0, cache=1 ):
        """
        Get list to map from any atom to a continuous residue numbering
        (starting with 0). A new residue is assumed to start whenever the
        'residue_number' or the 'residue_name' record changes between 2
        atoms. The mask is applied BEFORE looking for residue borders,
        i.e. it can change the residue numbering.

        See L{resList()} for an example of how to use the residue map.

        @param mask: [0000011111001111...] include atom: yes or no
                     len(atom_mask) == number of atoms in self.xyz
        @type  mask: list of int (1||0)
        @param force: recalculate map even if cached one is available (def 0)
        @type  force: 0||1
        @param cache: cache new map (def 1)
        @type  cache: 0||1

        @return: array [00011111122223333..], residue index for each
                 unmasked atom
        @rtype:  list of int
        """
        return self.index2map( self.resIndex( force=force,cache=cache ),
                               self.lenAtoms() )

    def resEndIndex( self ):
        """
        Get the position of the each residue's last atom.
        @return: index of the last atom of each residue
        @rtype: list of int
        """
        r = self.resIndex()
        N.concatenate( (r[1:], [self.lenAtoms()]) ) - 1 


    def __inferChainIndex( self ):
        result = []

        lastResidue = -100
        lastChainID = None
        lastSegID = None

        chn_ids = self.aProfiles['chain_id']
        seg_ids = self.aProfiles['segment_id']
        res_nrs = self.aProfiles['residue_number']
        ter_atm = self.aProfiles['after_ter']

        for i in range( self.lenAtoms() ):

            if chn_ids[i] != lastChainID or \
               seg_ids[i] != lastSegID   or \
               res_nrs[i] <  lastResidue or \
               ter_atm[i]:

                result.append( i )

            lastResidue = res_nrs[i]
            lastChainID = chn_ids[i]
            lastSegID   = seg_ids[i]

        return N.array( result, N.Int )


    def chainIndex( self, breaks=0, maxDist=None, force=0, cache=1 ):
        """
        Get indices of first atom of each chain.

        @param breaks: split chains at chain breaks (def 0)
        @type  breaks: 1||0
        @param maxDist: (if breaks=1) chain break threshold in Angstrom
        @type  maxDist: float
        @param force  : re-analyze residue numbering, chain and segids to
                        find chain boundaries, use with care! (def 0)
        @type  force  : 1||0
        @param cache  : cache result if new (breaks are not cached) (def 1)
        @type  cache  : 1||0

        @return: array (1 x N_chains) of int
        @rtype: list of int
        """
        ## fast track
        if not breaks and not force and self._chainIndex is not None:
            return self._chainIndex
        
        r = self._chainIndex
        
        if self._chainIndex is None or force:
            r = self.__inferChainIndex()

            if cache:
                self._chainIndex = r

        if breaks:
            break_pos = self.chainBreaks( breaks_only=1, maxDist=None )
            break_pos = break_pos + 1  ## chainBreaks reports last atom of each chain
            r = mathUtils.union( break_pos, r )
            r.sort()

        return N.array( r, N.Int )


    def chainMap( self, breaks=0, maxDist=None ):
        """
        Get chain index of each atom. A new chain is started between 2 atoms if
        the chain_id or segment_id changes, the residue numbering jumps back or
        a TER record was found.

        @param breaks: split chains at chain breaks (def 0)
        @type  breaks: 1||0
        @param maxDist: (if breaks=1) chain break threshold in Angstrom
        @type  maxDist: float

        @return: array 1 x N_atoms of int, e.g. [000000011111111111122222...]
        @rtype: list of int
        """
        return self.index2map( self.chainIndex( breaks=breaks, maxDist=maxDist ),
                               self.lenAtoms() )
            

    def chainBreaks( self, breaks_only=1, maxDist=None, force=0 ):
        """
        Identify discontinuities in the molecule's backbone.

        @param breaks_only: don't report ends of regular chains (def 1)
        @type  breaks_only: 1|0
        @param maxDist: maximal distance between consequtive residues
                        [ None ] .. defaults to twice the average distance
        @type  maxDist: float

        @return: atom indices of last atom **before** a probable chain break
        @rtype: list of int
        """
        if self.__chainBreaks is not None and not force and \
               maxDist is None and breaks_only:
            r = self.__chainBreaks

        else:
        
            i_bb = N.nonzero( self.maskBB() )
            bb   = self.take( i_bb )
            bb_ri= bb.resIndex()
            xyz = [ bb.xyz[ bb_ri[i] : bb_ri[i+1] ] for i in range(len(bb_ri)-1) ]
            xyz +=[ bb.xyz[ bb_ri[-1]: len(bb) ] ]

            ## get distance between last and first backbone atom of each residue
            dist = lambda a,b: N.sqrt( N.sum( N.power( a - b ,2 ) ) )

            d = [ dist( xyz[i][-1], xyz[i+1][0] ) for i in range(len(bb_ri)-1) ]

            ## get distances above mean
            cutoff = maxDist or MLab.mean(d)*2
            breaks = N.nonzero( N.greater( d, cutoff ) )

            if len(breaks) == 0:
                return N.array( [] )

            ri = self.resIndex()
            ri_to_e = {}
            for i in range( len(ri)-1 ):
                ri_to_e[ ri[i] ] = ri[ i+1 ]-1

            ## map back to the original atom indices
            r = [ ri_to_e[ i_bb[ bb_ri[i] ] ] for i in breaks ]

            if breaks_only:
                ri = self.chainIndex( breaks=0 )
                r = [ x for x in r if not x+1 in ri ]

                if maxDist is None:
                    self.__chainBreaks = r

        return N.array( r )


    def removeRes( self, resname ):
        """
        Remove all atoms with a certain residue name.

        @param resname: name of residue to be removed
        @type  resname: str OR list of str
        """
        resname = T.toList( resname )
        self.remove( self.maskFrom( 'residue_name', resname) )


    def rms( self, other, mask=None, mask_fit=None, fit=1, n_it=1 ):
        """
        Rmsd between two PDBModels.

        @param other: other model to compare this one with
        @type  other: PDBModel
        @param mask: atom mask for rmsd calculation
        @type  mask: list of int
        @param mask_fit: atom mask for superposition (default: same as mask)
        @type  mask_fit: list of int
        @param fit: superimpose first (default 1)
        @type  fit: 1||0
        @param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence, kicking out outliers on the way
        @type  n_it: int

        @return: rms in Angstrom
        @rtype: float
        """
        x, y = self.getXyz(), other.getXyz()

        if mask_fit is None: mask_fit = mask

        if fit:

            fx, fy = x, y

            if mask_fit is not None:
                fx = N.compress( mask_fit, x, 0 )
                fy = N.compress( mask_fit, y, 0 )
            ## find transformation for best match
            r, t = rmsFit.match( fx, fy, n_iterations=n_it )[0]

            ## transform coordinates
            y = N.dot(y, N.transpose(r)) + t

        if mask is not None:
            x = N.compress( mask, x, 0 )
            y = N.compress( mask, y, 0 )

        ## calculate row distances
        d = N.sqrt(N.sum(N.power(x - y, 2), 1))

        return N.sqrt( N.average(d**2) )


    def transformation( self, refModel, mask=None, n_it=1,
                        z=2, eps_rmsd=0.5, eps_stdv=0.05,
                        profname='rms_outlier'):
        """
        Get the transformation matrix which least-square fits this model
        onto the other model.

        @param refModel: reference PDBModel
        @type  refModel: PDBModel
        @param mask: atom mask for superposition
        @type  mask: list of int
        @param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence
        @type  n_it: int        
        @param z: number of standard deviations for outlier definition
                  (default 2)
        @type  z: float
        @param eps_rmsd: tolerance in rmsd (default 0.5)
        @type  eps_rmsd: float
        @param eps_stdv: tolerance in standard deviations (default 0.05)
        @type  eps_stdv: float
        @param profname: name of new atom profile getting outlier flag
        @type  profname: str
        @return: array(3 x 3), array(3 x 1) - rotation and translation matrices
        @rtype: array, array
        """
        x, y = refModel.getXyz(), self.getXyz()

        if mask is not None:
            x = N.compress( mask, x, 0 )
            y = N.compress( mask, y, 0 )
            outlier_mask = N.zeros( N.sum(mask)  )

        else:
            outlier_mask = N.zeros( len(self) )

        r, iter_trace = rmsFit.match( x, y, n_iterations=n_it, z=z,
                                      eps_rmsd=eps_rmsd, eps_stdv=eps_stdv)

        N.put( outlier_mask, iter_trace[-1][-1], 1 )

        if n_it != 1:
            self.aProfiles.set( profname, outlier_mask, mask=mask,
                                 default=1,
                                 comment='outliers in last iterative fitting',
                                 n_iterations=len( iter_trace ) )
        return r


    def transform( self, *rt ):
        """
        Transform coordinates of PDBModel.

        @param rt: rotational and translation array:
                   array( 4 x 4 ) OR array(3 x 3), array(3 x 1)
        @type  rt: array OR array, array

        @return: PDBModel with transformed coordinates
        @rtype: PDBModel
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
        result.setXyz( N.dot( self.getXyz(), N.transpose(r) ) + t  )

        return result


    def fit( self, refModel, mask=None, n_it=1,
             z=2, eps_rmsd=0.5, eps_stdv=0.05,
             profname='rms_outlier'):
        """
        Least-square fit this model onto refMode

        @param refModel: reference PDBModel
        @type  refModel: PDBModel
        @param mask: atom mask for superposition
        @type  mask: list of int (1||0)
        @param n_it: number of fit iterations::
                       1 - classic single fit (default)
                       0 - until convergence
        @type  n_it: int        
        @param z: number of standard deviations for outlier definition
                  (default 2)
        @type  z: float
        @param eps_rmsd: tolerance in rmsd (default 0.5)
        @type  eps_rmsd: float
        @param eps_stdv: tolerance in standard deviations (default 0.05)
        @type  eps_stdv: float
        @param profname: name of new atom profile containing outlier flag
        @type  profname: str

        @return: PDBModel with transformed coordinates
        @rtype: PDBModel
        """
        return self.transform(
            self.transformation( refModel, mask, n_it, eps_rmsd=eps_rmsd,
                                 eps_stdv=eps_stdv, profname=profname ) )


    def magicFit( self, refModel, mask=None ):
        """
        Superimpose this model onto a ref. model with similar atom content.
        magicFit( refModel [, mask ] ) -> PDBModel (or subclass )

        @param refModel: reference PDBModel
        @type  refModel: PDBModel
        @param mask: atom mask to use for the fit
        @type  mask: list of int (1||0)

        @return: fitted PDBModel or sub-class
        @rtype: PDBModel
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


    def centered( self, mask=None ):
        """
        Get model with centered coordinates.

        @param mask: atom mask applied before calculating the center
        @type  mask: list of int (1||0)

        @return: model with centered coordinates
        @rtype: PDBModel
        """
        r = self.clone()
        if mask is None: mask = N.ones( len(self) )

        avg = N.average( N.compress( mask, r.getXyz(), 0 ) )

        r.setXyz( r.getXyz() - avg )

        return r


    def center( self, mask=None ):
        """
        Geometric centar of model.

        @param mask: atom mask applied before calculating the center
        @type  mask: list of int (1||0)

        @return: xyz coordinates of center
        @rtype: (float, float, float)
        """
        if mask is None:
            return N.average( self.getXyz() )

        return N.average( N.compress( mask, self.getXyz(), axis=0 ) )


    def centerOfMass( self ):
        """
        Center of mass of PDBModel.

        @return: array(N.Float32)
        @rtype: (float, float, float)
        """
        M = self.masses()
        return mathUtils.wMean( self.getXyz(), M )


    def masses( self ):
        """
        Collect the molecular weight of all atoms in PDBModel.

        @return: 1-D array with mass of every atom in 1/12 of C12 mass.
        @rtype: array of floats

        @raise PDBError: if the model contains elements of unknown mass
        """
        try:
            M = [ molUtils.atomMasses[e] for e in self.aProfiles['element'] ]

        except KeyError, why:
            raise PDBError('Cannot find mass for '+str(why))

        return N.array( M )


    def mass( self ):
        """
        Molecular weight of PDBModel.

        @return: total mass in 1/12 of C12 mass
        @rtype: float

        @raise PDBError: if the model contains elements of unknown mass
        """
        return N.sum( self.masses() )


    def residusMaximus( self, atomValues, mask=None ):
        """
        Take list of value per atom, return list where all atoms of any
        residue are set to the highest value of any atom in that residue.
        (after applying mask)

        @param atomValues: values per atom
        @type  atomValues: list 
        @param mask: atom mask
        @type  mask: list of int (1||0)

        @return: array with values set to the maximal intra-residue value
        @rtype: array of float
        """
        if mask is None:
            mask = N.ones( len(atomValues) )

        ## eliminate all values that do not belong to the selected atoms
        masked = atomValues * mask

        result = []

        ## set all atoms of each residue to uniform value
        for res in range( 0, self.lenResidues() ):

            ## get atom entries for this residue
            resAtoms = N.compress( N.equal( self.resMap(), res ), masked )

            ## get maximum value
            masterValue = max( resAtoms )

            result += resAtoms * 0.0 + masterValue

        return N.array( result, N.Float32 )


    def argsort( self, cmpfunc=None ):
        """
        Prepare sorting atoms within residues according to comparison function.

        @param cmpfunc: function( m.aProfiles[i], m.aProfiles[j] ) -> -1, 0, +1
        @type  cmpfunc: function

        @return: suggested position of each atom in re-sorted model 
                 ( e.g. [2,1,4,6,5,0,..] )
        @rtype: list of int
        """
        ## by default sort alphabetically by atom name
        cmpfunc = cmpfunc or ( lambda a1, a2: cmp(a1['name'],a2['name']) )

        result = []
        pos = 0

        ## get sort list for each residue
        for resAtoms in self.resList():

            resIndex = range(0, len( resAtoms ) )

            ## convert atom-based function into index-based function
            f_cmp = lambda i,j,atoms=resAtoms : cmpfunc( atoms[i], atoms[j] )

            ## get sortMap for this residue (numbering starts with 0)
            resIndex.sort( f_cmp )

            ## add first residue atom's position in self.atoms to each entry 
            resIndex = map( lambda i, delta=pos: i + delta, resIndex )

            ## concatenate to result list
            result += resIndex

            pos += len( resAtoms )

        return result


    def sort( self, sortArg=None ):
        """
        Apply a given sort list to the atoms of this model.

        @param sortArg: comparison function
        @type  sortArg: function

        @return: copy of this model with re-sorted atoms (see Numeric.take() )
        @rtype: PDBModel
        """
        sortArg = sortArg or self.argsort()

        r = self.take( sortArg )

        return r


    def unsort( self, sortList ):
        """
        Undo a previous sorting on the model itself (no copy).

        @param sortList: sort list used for previous sorting.
        @type  sortList: list of int

        @return: the (back)sort list used ( to undo the undo...)
        @rtype: list of int

        @raise PDBError: if sorting changed atom number
        """
        ## get new sort list that reverts the old one
        backSortList = range(0, self.lenAtoms() )
        backSortList.sort( lambda i,j, l=sortList: cmp(l[i],l[j]) )

        ## get re-sorted PDBModel copy
        self.keep( backSortList )

        return backSortList


    def atomNames(self, start=None, stop=None):
        """
        Return a list of atom names from start to stop RESIDUE index

        @param start: index of first residue 
        @type  start: int
        @param stop: index of last residue
        @type  stop: int

        @return: ['C','CA','CB' .... ]
        @rtype: list of str
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

        return self.aProfiles['name'][i:j]


    def __testDict_and( self, dic, condition ):
        """
        Test if B{all} key-value pairs of condition are matched in dic

        @param condition: {..}, key-value pairs to be matched
        @type  condition: dictionary
        @param dic: {..}, dictionary to be tested
        @type  dic: dictionary

        @return: 1|0, 1 if all key-value pairs of condition are matched in dic
        @rtype: 1|0
        """
        for k,v in condition.items():
            if dic.get( k, None ) not in v:
                return 0
        return 1


    def __testDict_or( self, dic, condition ):
        """
        Test if B{any} key-value pairs of condition are matched in dic

        @param condition: {..}, key-value pairs to be matched
        @type  condition: dictionary
        @param dic: {..}, dictionary to be tested
        @type  dic: dictionary

        @return: 1|0, 1 if any key-value pairs of condition are matched in dic
        @rtype: 1|0
        """
        for k,v in condition.items():
            if dic.get( k, None ) in v:
                return 1
        return 0


    def filterIndex( self, mode=0, **kw ):
        """
        Get atom positions that match a combination of key=values.
        E.g. filter( chain_id='A', name=['CA','CB'] ) -> index

        @param mode: 0 combine with AND (default), 1 combine with OR
        @type  mode: 0||1
        @param kw: combination of atom dictionary keys and
                   values/list of values that will be used to filter
        @type  kw: filter options, see example

        @return: sort list
        @rtype: list of int
        """
        ## cache to minimize function lookup
        atoms = self.aProfiles.toDicts()
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

        @param mode: 0 combine with AND (default), 1 combine with OR
        @type  mode: 0||1
        @param kw: combination of atom dictionary keys and
                   values/list of values that will be used to filter
        @type  kw: filter options, see example

        @return: filterd PDBModel
        @rtype: PDBModel
        """
        return self.take( self.filterIndex( mode=mode, **kw ) )


    def equals(self, ref, start=None, stop=None):
        """
        Compares the residue and atom sequence in the given range.
        Coordinates are not checked, other profiles are not checked.

        @param start: index of first residue 
        @type  start: int
        @param stop: index of last residue
        @type  stop: int 

        @return: [ 1||0, 1||0 ],
                 first position sequence identity 0|1,
                 second positio atom identity 0|1
        @rtype: list if int
        """
        ## By default compare all residues
        start = start or 0
        if stop == None:    ## don't use stop = stop or .. stop might be 0!
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


##     def equalAtoms( self, ref ):
##         """
##         Apply to SORTED models without HETATOMS. Coordinates are not checked.

##         @note: in some rare cases m1.equalAtoms( m2 ) gives a different result
##                than m2.equalAtoms( m1 ). This is due to the used
##                SequenceMatcher class.      
##         @todo: option to make sure atoms are also in same order

##         @param ref: reference PDBModel
##         @type  ref: PDBModel

##         @return: (mask, mask_ref), two atom masks for all equal (1) atoms
##                   in models
##         @rtype: (array, array)
##         """
##         ## compare sequences
##         seqMask, seqMask_ref = match2seq.compareModels(self, ref)

##         ## get residue mask on atom level
##         mask = self.res2atomMask(seqMask)
##         mask_ref = ref.res2atomMask(seqMask_ref)

##         ## get list of matching RESIDUES
##         equal = N.nonzero(seqMask)
##         equal_ref = N.nonzero(seqMask_ref)

##         ## check that all atoms are equal in matching residues
##         for i in range(0, len(equal)):

##             ## atom name lists for current residue
##             aa = self.atomNames(equal[i],equal[i])
##             aa_ref = ref.atomNames(equal_ref[i],equal_ref[i])

##             ## starting atom of current residue
##             ind = self.resIndex()[ equal[i] ]
##             ind_ref = ref.resIndex()[ equal_ref[i] ]

##             ## check that the atoms are the same
##             for j in range( len(aa) ):

##                 try:
##                     mask[ind] = aa[j] == aa_ref[j]
##                     ind += 1

##                 except IndexError:
##                     mask[ind] = 0
##                     ind += 1

##             for j in range( len(aa_ref) ):

##                 try:
##                     mask_ref[ind_ref] = aa_ref[j] == aa[j]
##                     ind_ref += 1
##                 except IndexError:
##                     mask_ref[ind_ref] = 0
##                     ind_ref += 1


##         return mask, mask_ref


    def compareAtoms( self, ref ):
        """
        Get list of atom indices for this and reference model that converts
        both into 2 models with identical residue and atom content.

        E.g.
         >>> m2 = m1.sort()    ## m2 has now different atom order
         >>> i2, i1 = m2.compareAtoms( m1 )
         >>> m1 = m1.take( i1 ); m2 = m2.take( i2 )
         >>> m1.atomNames() == m2.atomNames()  ## m2 has again same atom order

        @return: indices, indices_ref
        @rtype: ([int], [int])
        """
        ## compare sequences
        seqMask, seqMask_ref = match2seq.compareModels(self, ref)

        ## get list of matching RESIDUES
        equal = N.nonzero(seqMask)
        equal_ref = N.nonzero(seqMask_ref)

        result, result_ref = [], []

        ## check that all atoms are equal in matching residues
        for i in range(0, len(equal)):

            ## atom name lists for current residue
            aa = self.atomNames( equal[i],equal[i] )
            aa_ref = ref.atomNames( equal_ref[i],equal_ref[i] )

            ## starting atom of current residue
            ind = self.resIndex()[ equal[i] ]
            ind_ref = ref.resIndex()[ equal_ref[i] ]

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


    def __chainFraction( self, chain, ref ):
        """
        Look how well a given chain matches a continuous stretch of residues
        in ref.

        @param chain: chain index
        @type  chain: int
        @param ref: reference PDBModel
        @type  ref: PDBModel

        @return: average relative length of matching chain fragments
        @rtype: float
        """
        m = self.takeChains([chain])

        m_cast = m.take( m.compareAtoms( ref )[0] )

        f = 1. * len( m_cast ) / m_cast.lenChains( breaks=1, maxDist=5.)
        f = f / len( m )

        return f


    def compareChains( self, ref, breaks=0, fractLimit=0.2 ):
        """
        Get list of corresponding chain indices for this and reference model.
        Use takeChains() to create two models with identical chain content and
        order from the result of this function.

        @param ref: reference PDBModel
        @type  ref: PDBModel
        @param breaks: look for chain breaks in backbone coordinates
        @type  breaks: 1||0
        @param fractLimit: 
        @type  fractLimit: float

        @return: chainIndices, chainIndices_ref
        @rtype: ([int], [int])
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



#############
##  TESTING        
#############

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    #: load test PDB once into class rather than 3 times into every instance
    MODEL = None

    def prepare(self):
        import tempfile

        ## loading output file from X-plor
        t = time.clock()
        self.MODEL = self.MODEL or B.PDBModel( T.testRoot()+'/com/1BGS.pdb')
        self.m = self.MODEL
        self.fout_pdb = tempfile.mktemp( '_test1.pdb' )
        self.fout1 = tempfile.mktemp( '_test1.model' )
        self.fout2 = tempfile.mktemp( '_test2.model' )
        print "prepare: ", time.clock() - t

    def cleanUp( self ):
        T.tryRemove( self.fout1 )
        T.tryRemove( self.fout2 )
        T.tryRemove( self.fout_pdb )

    def test_removeRes(self):
        """PDBModel.removeRes test"""
        t = time.clock()
        self._m = self.m.clone()
        self._m.removeRes(['TIP3', 'HOH'])
        self.assertEqual( len(self._m), 1968)
        self.assertAlmostEqual( self._m.mass(), 21325.90004, 3 )
        print "removeRes: ", time.clock() - t

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
            print 'The molecule consists of %i chains'% self.m.lenChains()
            print '\tChainId \tFirst atom'
            for i in chainIdx:
                print '\t%s \t\t%i'%(self._m.aProfiles['chain_id'][i], int(i))

        ## iterate over all chains
        for c in range( 0, len( chainIdx ) ):

            if self.local:
                print "chain ", c, " starts with ", 
                print self._m.aProfiles['residue_name'][ chainIdx[c] ],

                print " and has sequence: "

        
            ## mask out atoms of all other chains
            chainMask  = N.equal( self._m.chainMap( breaks=1 ), c )
            if self.local:
                print self._m.sequence( chainMask )

        self.assertEqual( self._m.lenChains(), 4)


    def test_sorting(self):
        """PDBModel sorting test"""
        if self.local:
            print "sorting atoms alphabetically..."
        m2 = self.m.compress( self.m.maskProtein() )
        sort = m2.argsort()
        m2 = m2.sort( sort )

        self.assertAlmostEqual( N.sum( m2.centerOfMass() ),  23.1032009125)


    def test_chainHandling(self):
        """PDBModel chain handling and writing test"""
        self.m4 = B.PDBModel( T.testRoot()+'/com/1BGS_original.pdb')
        self.assertEqual( self.m4.lenChains(), 9 )
        self.assertEqual( self.m4.lenChains( breaks=1 ), 138 )
        self.m4.writePdb( self.fout_pdb, ter=2 )

    def test_rename(self):
        """PDBModel renameAmberRes tests"""
        self.m3 = B.PDBModel( T.testRoot()+'/amber/1HPT_0dry.pdb')

        n_cyx = self.m3.aProfiles['residue_name'].count('CYX')
        n_hid = self.m3.aProfiles['residue_name'].count('HID')
        n_hip = self.m3.aProfiles['residue_name'].count('HIP')
        n_hie = self.m3.aProfiles['residue_name'].count('HIE')
        n_hix = n_hid + n_hie + n_hip
        
        self.m3.renameAmberRes()

        self.assertEqual(n_cyx, self.m3.aProfiles['residue_name'].count('CYS'))
        self.assertEqual(n_hix, self.m3.aProfiles['residue_name'].count('HIS'))


    def test_sourceHandling(self):
        """PDBModel source / disconnection tests"""
        self._m = self.m.clone()
        anames = self.m.aProfiles['name']
        xyz0   = self.m.getXyz()[0]

        self._m.slim()
        
        ## _m2 uses _m1 as source
        self._m2 = B.PDBModel( self._m )
        l1 = self._m2.aProfiles['name']
        self.assertEqual( l1, anames )
        
        ## remove unchanged profiles and coordinates
        self._m2.slim()
        
        ## fetch them again from source (of source)
        l2 = self._m2.aProfiles['name']
        self.assertEqual( l2, anames )
        
        ## disconnect _m from PDB file source
        self._m.saveAs( self.fout2 )

        self._m2.slim()
        self._m.slim()
        
        ## this should now trigger the reloading of fout2
        self.assertEqual( self._m2.aProfiles['name'], anames )
        self.assert_( N.all( self._m2.getXyz()[0] == xyz0) )
        
        ## after disconnection, slim() should not have any effect
        self._m2.disconnect()
        self._m2.slim()
        self.assert_( self._m2.aProfiles.profiles['name'] is not None )
        
        
        
        
def clock( s, ns=globals() ):
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

##     import Biskit as B

##     print "Loading"
##     m = B.PDBModel( T.testRoot()+'/com/1BGS.pdb')
##     print "Loaded"
##     m.slim()
    
##     m_ = B.PDBModel( m )
##     l1 = m_.aProfiles['name']
##     m_.slim()
##     l2 = m_.aProfiles['name']
    
##     print "Done"
        
    BT.localTest()
