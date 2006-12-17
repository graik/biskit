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
Store and manipulate coordinates and atom information.
"""

import tools as t
import molUtils
import mathUtils
import match2seq
import rmsFit
from LocalPath import LocalPath
from Errors import BiskitError
from Biskit import EHandler
from ProfileCollection import ProfileCollection, ProfileError
from PDBParserFactory import PDBParserFactory

import Numeric as N
import MLab
import os, sys
import copy
import time
import string
import types
import Scientific.IO.PDB as IO

from multiarray import arraytype


class PDBProfiles( ProfileCollection ):
    def version( self ):
        return ProfileCollection.version(self) + '; PDBModel $Revision$'

class PDBError(BiskitError):
    pass

class PDBModel:
    """
    Store and manipulate coordinates and atom infos stemming from a
    PDB file. Coordinates are stored in the Numeric array 'xyz'; the
    additional atom infos from the PDB are stored in the list of
    dictionaries 'atoms'.
    Methods are provided to remove items from  both atoms and xyz
    simultaniously, to restore the PDB file, to get masks for certain
    atom types/names/residues, to iterate over residues, sort atoms etc.

    Atom- or residue-related values can be put into 'profile' arrays.
    See setAtomProfile() and setResProfile(). Any reordering or removal of 
    atoms is also applied to the profiles so that they should always match
    the current atom/residue order.

    The object remembers its source (a PDB file or a PDBModel in memory
    or a pickled PDBModel on disc) and keeps track of whether atoms, xyz, or
    (which) profiles have been changed with respect to the source.
    The slim() method is called before pickling a PDBModel. It sets to None
    the atoms and/or xyz array or any profiles if they have not been changed
    since beeing read from a source on disc. The change status of xyz and
    atoms is reported by isChanged() (with respect to the direct source) and
    isChangedFromDisc() (with respect to the source of the source... just in
    case).
    You can trick this mechanism by setting atomsChanged or xyzChanged back
    to 0 if you want to make only temporary changes that are lost after a
    call to slim().

    Additional infos about the model can be put into a dictionary 'info'.

    @todo: clean up, rename returnPdb()
    @todo: clean up the __terAtoms mess
    @todo: use Profiles instead of space consuming atom dictionaries
    """

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
        self.atoms = None

        ## save atom-/residue-based values
        self.rProfiles = PDBProfiles()
        self.aProfiles = PDBProfiles()

        ## pre-defined atom masks, calculated when first needed
        self.caMask = None
        self.bbMask = None
        self.heavyMask = None

        ## prepare caching
        self.__resMap = None
        self.__resIndex = None

        ## positions of atoms followed by TER in PDB file
        self.__terAtoms = []

        if noxyz:
            ## trick update to leave xyz untouched
            self.xyz = 0

        ## monitor changes of coordinates and atoms
        self.xyzChanged = 0
        self.atomsChanged = 0

        self.forcePickle = 0

        ## version as of creation of this object
        self.initVersion = self.version()

        ## to collect further informations
        self.info = { 'date':t.dateSortString() }

        if source <> None:

            self.update( skipRes=skipRes, lookHarder=1 )

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

    def __defaults(self ):
        """
        backwards compatibility to earlier pickled models
        """

        self.__resMap = getattr( self, '_PDBModel__resMap', None)
        self.__resIndex = getattr( self, '_PDBModel__resIndex', None)

        if type( self.source ) == str:
            self.source = LocalPath( self.source )

        self.__validSource = getattr( self, '_PDBModel__validSource', 0)

        self.initVersion = getattr( self, 'initVersion', 'old PDBModel')

        ## convert old profile dictionaries into new ProfileCollections
        if 'resProfiles' in self.__dict__:
            self.rProfiles=PDBProfiles(
                                   profiles=getattr(self,'resProfiles',{}),
                                   infos=getattr(self,'resProfiles_info',{}) )
            try:
                del self.resProfiles; del self.resProfiles_info
            except: pass

        if 'atomProfiles' in self.__dict__:
            self.aProfiles=PDBProfiles(
                                   profiles=getattr(self,'atomProfiles',{}),
                                   infos=getattr(self,'atomProfiles_info',{}) )
            try:
                del self.atomProfiles; del self.atomProfiles_info
            except: pass

        ## if there are not even old profiles...
        if not getattr( self, 'aProfiles', None):
            self.aProfiles = PDBProfiles()
        if not getattr( self, 'rProfiles', None):
            self.rProfiles = PDBProfiles()

        self.forcePickle = getattr( self, 'forcePickle', 0 )

        self.info = getattr( self, 'info', { 'date':t.dateSortString() } )

        ## not simplified to 1 line because __pdbTer would otherwise always
        ## been executed
        if getattr( self, '_PDBModel__terAtoms', -1) == -1:
            self.__terAtoms = self.__pdbTer(1)

        ## fix previous bug: slim() was creating a self.xyx instead of xyz
        if getattr( self, 'xyx', 0 ):
            del self.xyx


    def update( self, skipRes=None, lookHarder=0, force=0 ):
        """
        Read coordinates, atoms, fileName, etc. from PDB or
        pickled PDBModel - but only if they are currently empty.
        The atomsChanged and xyzChanged flags are not changed.
        
        @param skipRes: names of residues to skip if updating from PDB
        @type  skipRes: list of str
        @param lookHarder: 0(default): update only existing profiles
        @type  lookHarder: 0|1
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
        parser.update(self, source, skipRes=skipRes, lookHarder=lookHarder )
        
        
    def __pdbTer( self, rmOld=0 ):
        """
        @param rmOld: 1, remove after_ter=0 flags from all atoms
        @type  rmOld: 1||0

        @return: list of atom indices that are followed by a TER record
                 (marked with 'after_ter' flag of the next atom
                 by __collectAll).
        @rtype: list of int
        """
        atoms = self.getAtoms()

        ## all atoms preceeding a TER record
        ater = [ i-1 for i in range( self.lenAtoms() )
                 if atoms[i].get('after_ter',0)]

        ## remove after_ter=0 flags from all atoms
        if rmOld:
            for a in self.atoms:
                if a.has_key('after_ter') and a['after_ter']==0:
                    del a['after_ter']

        ## remove negative indices
        return N.compress( N.greater( ater, -1 ), ater )


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


    def setAtoms( self, atoms ):
        """
        Replace atoms list of dictionaries.
        self.__terAtoms is re-created from the 'after_ter' records in atoms

        @param atoms: list of dictionaries as returned by PDBFile.readLine()
                      (nearly, for differences see L{__collectAll()} )
        @type  atoms: list of dictionaries

        @return: [ dict ], old atom dictionaries
        @rtype: list of dict
        """
        old = self.atoms
        self.atoms = atoms

        changedNow = not ( self.atoms == old )

        if self.atomsChanged:
            self.__resMap = self.__resIndex = None
            self.caMask = None
            self.bbMask = None
            self.heavyMask = None
            self.__terAtoms = self.__pdbTer()

        self.atomsChanged = self.atomsChanged or changedNow

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

        @return: xyz-coordinates, N.array( 3 x N_atoms, 'f' )
        @rtype: array 
        """
        if self.xyz is None:
            self.update( force=1 )

        if self.xyz is None:
            return N.array( [], 'f' )

        if mask is None:
            return self.xyz

        return N.compress( mask, self.xyz, 0 )


    def getAtoms( self, mask=None ):
        """
        Get atom dictionaries, fetch from source PDB or pickled PDBModel,
        if necessary. See also L{__collectAll()}.

        @param mask: atom mask
        @type  mask: list of int OR array of 1||0

        @return: atom dictionary, list of dictionaries from PDBFile.readline()
        @rtype: list of dic
        """
        if self.atoms is None:
            self.update( force=1 )

        if self.atoms is None:
            return []

        if mask is None:
            return self.atoms

        return [ self.atoms[i] for i in N.nonzero( mask ) ]


    def setResProfile( self, name, prof, mask=None, default=None, asarray=1,
                       comment=None, **moreInfo):
        """
        Add/override residue-based profile.

        @param name: name to access profile
        @type  name: str
        @param prof: list/array of values
        @type  prof: list OR array
        @param mask: list/array 1 x N_residues of 0|1, N.sum(mask)==len(prof)
        @type  mask: list OR array
        @param default: value for masked residues 
                        default: None for lists, 0 for arrays
        @type  default: any
        @param asarray: store as list (0), as array (2) or store numbers as
                        array but everything else as list (1) default: 1
        @type  asarray: 0|1|2
        @param comment: goes into aProfiles_info[name]['comment']
        @type  comment: str
        @param moreInfo: additional key-value pairs for aProfiles_info[name]
        @type  moreInfo: (key, value)

        @raise ProfileError: if length of prof != N_atoms
        """
        self.rProfiles.set( name, prof, mask=mask, default=default,
                            asarray=asarray, comment=comment, **moreInfo )

    def setAtomProfile( self, name, prof, mask=None, default=None, asarray=1,
                        comment=None, **moreInfo):
        """
        Add/override atom-based profile.

        @param name: name to access profile
        @type  name: str
        @param prof: list/array of values
        @type  prof: list OR array
        @param mask: list/array 1 x N_residues of 0 || 1, profile items
        @type  mask: list OR array
        @param default: value for masked residues
                        default: None for lists, 0 for arrays
        @type  default: any
        @param asarray: store as list (0), as array (2) or store numbers as
                        array but everything else as list (1) default: 1
        @type  asarray: 0|1|2
        @param comment: goes into aProfiles_info[name]['comment']
        @type  comment: str
        @param moreInfo: additional key-value pairs for aProfiles_info[name]
        @type  moreInfo: (key, value)

        @raise ProfileError: if length of prof
        """
        self.aProfiles.set( name, prof, mask=mask, default=default,
                            asarray=asarray, comment=comment, **moreInfo )


    def resProfile( self, name, default=None ):
        """
        resProfile( profile_name ) -> array 1 x N_res with residue values

        @param name: name to access profile
        @type  name: str

        @return: residue profile array 1 x N_res with residue values
        @rtype:  array 

        @raise ProfileError: if rProfiles contains no entry for |name|
        """
        return self.rProfiles.get( name, default=default )


    def atomProfile(self, name, default=None ):
        """
        atomProfile( profile_name ) -> array 1 x N_atoms with atom values

        @param name: name to access profile
        @type  name: str

        @return: atom profile array 1 x N_atoms with atom values
        @rtype:  array

        @raise ProfileError: if aProfiles contains no entry for |name|
        """
        return self.aProfiles.get( name, default=default )


    def profile( self, name, default=None, lookHarder=0 ):
        """
        profile( name, lookHarder=0) -> atom or residue profile

        @param name: name to access profile
        @type  name: str        
        @param default: default result if no profile is found,
                        if None, raise exception
        @type  default: any
        @param lookHarder: update from source before reporting missing profile
        @type  lookHarder: 0||1

        @raise ProfileError: if neither atom- nor rProfiles
                                   contains |name|
        """
        if lookHarder and not name in self.aProfiles and \
               not name in self.rProfiles:
            self.update()

        if name in self.aProfiles:
            return self.aProfiles[ name ]

        return self.rProfiles.get( name, default=default )


    def profileInfo( self, name, lookHarder=0 ):
        """
        profileInfo( name ) -> dict with infos about profile

        @param name: name to access profile
        @type  name: str       
        @param lookHarder:update from source before reporting missing profile::
                           Guaranteed infos: 'version'->str,
                                             'comment'->str,
                                             'changed'->1||0
        @type  lookHarder: 0|1

        @raise ProfileError: if neither atom - nor rProfiles
                                   contains |name|
        """
        if lookHarder and not name in self.Profiles and \
               not name in self.rProfiles:
            self.update()

        if name in self.aProfiles:
            return self.aProfiles.getInfo( name )

        if name in self.rProfiles:
            return self.rProfiles.getInfo( name )

        raise ProfileError( 'No profile info found with name '+str(name))


    def setProfileInfo( self, name, **args ):
        """
        Add/Override infos about a given profile
        E.g. setProfileInfo('relASA', comment='new', params={'bin':'whatif'})

        @param name: name to access profile
        @type  name: str

        @raise ProfileError: if neither atom
        """
        d = self.profileInfo( name )
        for key, value in args.items():
            d[key] = value


    def removeProfile( self, *names ):
        """
        Remove residue or atom profile(s)
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


    def isChanged(self):
        """
        Tell if xyz or atoms have been changed compared to source file or
        source object (which can be still in memory).

        @return: (1,1)..both xyz and atoms field have been changed
        @rtype: (1||0, 1||0)
        """
        if (self.xyzChanged==0 or self.atomsChanged==0) \
           and ( self.source is not None and self.validSource() is None ):
            raise PDBError('Invalid source: '+str(self.source) )

        if self.validSource() is None:
            return (1,1)

        return (self.xyzChanged, self.atomsChanged)


    def isChangedFromDisc(self):
        """
        Tell whether xyz and atoms can currently be reconstructed from a
        source on disc. Same as isChanged() unless source is another not yet
        saved PDBModel instance that made changes relative to its own source
        ...

        @return: (1,1)..both xyz and atoms field have been changed
        @rtype: (1||0, 1||0)
        """
        if isinstance( self.source, PDBModel ):
            xChanged = self.isChanged()[0] or \
                       self.source.isChangedFromDisc()[0]
            aChanged = self.isChanged()[1] or \
                       self.source.isChangedFromDisc()[1]

            return xChanged, aChanged

        return self.isChanged()


    def profileChangedFromDisc(self, pname):
        """
        Check if profile has changed compared to source.

        @return: 1, if profile |pname| can currently not be
                 reconstructed from a source on disc.
        @rtype: int

        @raise ProfileError: if there is no atom or res profile with pname
        """
        if self.validSource() is None:
            return 1

        if isinstance( self.source, PDBModel ):
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
                self.rProfiles[ key ] = None

            if type( self.rProfiles[ key ] ) == arraytype:
                self.rProfiles[ key ] = self.rProfiles[key].tolist()

        for key in self.aProfiles:

            if not self.profileChangedFromDisc( key ):
                self.aProfiles[ key ] = None

            if type( self.aProfiles[ key ] ) == arraytype:
                self.aProfiles[ key ] = self.aProfiles[key].tolist()

    def slim( self ):
        """
        Remove xyz array and list of atoms if they haven't been changed and
        could hence be loaded from the source file (only if there is a source
        file...). Remove any unchanged profiles.
        B{AUTOMATICALLY CALLED BEFORE PICKLING}
        """
        ## remove atoms/coordinates if they are unchanged from an existing
        ## source
        ## override this behaviour with forcePickle
        if not self.forcePickle:

            xChanged, aChanged = self.isChangedFromDisc()

            if not xChanged:
                self.xyz = None

            if type( self.xyz ) == arraytype and self.xyz.typecode() != 'f':
                self.xyz = self.xyz.astype('f')

            if not aChanged:
                self.atoms = None

            if type( self.atoms ) == arraytype:
                self.atoms = self.atoms.tolist()

            self.__slimProfiles()

        self.__resMap = self.__resIndex = None
        self.caMask = self.bbMask = self.heavyMask = None
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
                if isinstance( self.source, PDBModel ):
                    self.__validSource = self.source
                else:
                    self.__validSource = None

        return self.__validSource


    def sourceFile( self ):
        """
        Name of pickled source or PDB file.

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

        self.xyzChanged, self.atomsChanged = 1, 1

        for p in self.rProfiles:
            self.setProfileInfo( p, changed=1 )
        for p in self.aProfiles:
            self.setProfileInfo( p, changed=1 )


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
        @type  xtable: dict

        @return: 1-letter-code AA sequence (based on first atom of each res).
        @rtype: str
        """
        mask = mask or N.ones( self.lenAtoms(), 'i' )

        firstAtm = N.zeros( self.lenAtoms() )
        N.put( firstAtm, self.resIndex(), 1 )

##        mask = mask * firstAtm * self.maskProtein()
        mask = mask * firstAtm

        l = [ a['residue_name'] for a in N.compress( mask, self.getAtoms()) ]

        return ''.join( molUtils.singleAA( l ) )


    def xplor2amber( self, change=1, aatm=1 ):
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

        residues = self.resList()

        if not change:
            residues = copy.deepcopy( residues )

        atoms = []

        for r in residues:

            for a in r:

                if aatm and len(a['name'])>2 and a['name'][0] in numbers:
                    a['name'] = a['name'][1:] + a['name'][0]

            resname = r[0]['residue_name']
            if resname == 'HIS':
                anames = [ a['name'] for a in r ]

                if 'HE2' in anames:  resname = 'HIE'
                if 'HD1' in anames:  resname = 'HID'
                if 'HE2' in anames and 'HD1' in anames:
                    resname = 'HIP'

            for a in r:
                a['residue_name'] = resname

            atoms.extend( r )

        return atoms


    def writePdb( self, fname, ter=1, amber=0, original=0, left=0, wrap=0,
                  headlines=None, taillines=None):
        """
        Save model as PDB file.

        @param fname: name of new file
        @type  fname: str
        @param ter: Option of how to treat the retminal record::
                    0, don't write any TER statements
                    1, restore original TER statements (doesn't work,
                         if preceeding atom has been deleted) [default]
                    2, put TER between all detected chains
                    3, as 2 but also detect and split discontinuous chains
        @type  ter: 0, 1, 2 or 3
        @param amber: amber formatted atom names
                      (ter=3, left=1, wrap=0) (default 0)
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

            atoms = self.getAtoms()
            if amber:
                atoms = self.xplor2amber( change=0 )
                ter = 3
                wrap = 0
                left = 1

            if ter == 2 or ter == 3:
                terIndex = N.array( self.chainIndex( breaks=(ter==3) )[1:] ) -1
            if ter == 1:
                terIndex = self.__terAtoms

            if headlines:
                for l in headlines:
                    f.writeLine( l[0], l[1] )

            if taillines:
                for l in taillines:
                    f.writeLine( l[0], l[1] )

            for i in range( len( atoms ) ):

                a = atoms[i]

                ## fetch coordinates Vector
                a['position'] = self.xyz[ i ]

                ## change atom name if requested...
                real_name = a['name']
                aname = real_name

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

                ## restore atom record
                a['name'] = real_name
                del( a['position'] )

                ## write TER line with details from previous atom
                if (ter>0 and  i in terIndex):
                    f.writeLine('TER', a )

            f.close()

        except:
            EHandler.error( "Error writing "+fname )


    def returnPdb( self, out=None, ter=1, headlines=None, taillines=None):
        """
        Restore PDB file from (possibly transformed) coordinates and pdb line
        dictionaries in self.atoms. This is an older version of writePdb that
        returns a list of PDB lines instead of writing to a file.

        @param out: stdout or None, if None a list is returned
        @type  out: stdout or None
        @param ter: Option of how to treat the retminal record::
                    0, don't write any TER statements
                    1, restore original TER statements (doesn't work,
                         if preceeding atom has been deleted)
                    2, put TER between all detected chains
        @type  ter: 0, 1 or 2
        @param headlines: [( str, dict or str)], list of record / data tuples::
                          e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
        @type  headlines: list of tuples
        @param taillines: same as headlines but appended at the end of file
        @type  taillines: list of tuples

        @return: [ str ], lines of a PDB file
        @rtype: list of strings
        """
        try:
            pdb_lst=[]
            i=0

            if ter == 2:
                terIndex = N.array( self.chainIndex()[1:] ) - 1
            if ter == 1:
                terIndex = self.__terAtoms

            if headlines:
                for l in headlines:
                    pdb_lst += [ '%6s%s\n'%(l[0],l[1]) ]

            if taillines:
                for l in taillines:
                    pdb_lst += [ '%6s%s\n'%(l[0],l[1]) ]

            for a in self.getAtoms():

                ## fetch coordinates Vector
                a['positionX'], a['positionY'], a['positionZ'] = self.xyz[ i ]

                ## switch original and short atom name
                a['name'], a['name_original'] = a['name_original'], a['name']

                ## write line
                atom_line = 'ATOM  %(serial_number)5i %(name)-4s %(residue_name)3s %(chain_id)1s%(residue_number)4i%(insertion_code)1s   %(positionX)8.3f%(positionY)8.3f%(positionZ)8.3f%(occupancy)6.2f%(temperature_factor)6.2f      %(segment_id)-4s%(element)2s\n'
                pdb_lst += [ atom_line%a ]

                ## switch back
                a['name'], a['name_original'] = a['name_original'], a['name']
                del( a['positionX'], a['positionY'], a['positionZ'] )

                ## write TER line with details from previous atom
                if (ter and  i in terIndex):
                    ter_line = 'TER   %(serial_number)5i      %(residue_name)3s %(chain_id)1s%(residue_number)4i%(insertion_code)1s\n'
                    pdb_lst += [ ter_line % self.getAtoms()[i] ]

                i += 1
            ## write to stdout or return list 
            if out == 'stdout':
                sys.stdout.writelines( pdb_lst )
            else:
                return pdb_lst

        except:
            EHandler.error( "PDBModel.returnPdb(): " )


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

            self.xyzChanged, self.atomsChanged = 0, 0

            for p in self.rProfiles:
                self.setProfileInfo( p, changed=0 )
            for p in self.aProfiles:
                self.setProfileInfo( p, changed=0 )

            t.Dump( self, str(path) )

        except IOError, err:
            raise PDBError("Can't open %s for writing." % t.absfile(str(path)))


    def maskF(self, atomFunction, numpy=1 ):
        """
        Create list whith result of atomFunction( atom ) for each atom.

        @param atomFunction: function( dict_from_PDBFile.readline() ),
                             true || false (Condition)
        @type  atomFunction: 1||0
        @param numpy: 1(default)||0, convert result to Numpy array of int
        @type  numpy: int

        @return: Numpy N.array( [0,1,1,0,0,0,1,0,..], 'i') or list
        @rtype: array or list
        """
        try:
            result = map( atomFunction, self.getAtoms() )
        except:

            ## fall-back solution: assign 0 to all entries that raise
            ## exception
            EHandler.warning("mask(): Error while mapping funtion "+
                       "to all atoms.")
            result = []

            for a in self.getAtoms():
                try:
                   result.append( atomFunction( a ) )
                ## put 0 if something goes wrong
                except :
                    EHandler.warning("mask(): Error while save-mapping ")
                    result.append(0)

        if numpy:
            return N.array( result )
        return result


    def maskCA( self, force=0 ):
        """
        Short cut for mask of all CA atoms.

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.caMask == None or force:
            self.caMask = self.maskF( lambda a: a['name'] == 'CA' )

        return self.caMask


    def maskBB( self, force=0 ):
        """
        Short cut for mask of all backbone atoms.

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.bbMask == None or force:
            self.bbMask = self.maskF( 
                lambda a: a['name'] in ['CA', 'C', 'N', 'O', 'H','OXT'] )

        return self.bbMask


    def maskHeavy( self, force=0 ):
        """
        Short cut for mask of all heavy atoms. ('element' <> H)

        @param force: force calculation even if cached mask is available
        @type  force: 0||1

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        if self.heavyMask == None or force:
            self.heavyMask = self.maskF( lambda a: a.get('element','') <> 'H' )

        return self.heavyMask

    def maskH( self ):
        """
        Short cut for mask of hydrogens. ('element' == H)

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return N.logical_not( self.maskHeavy() )


    def maskCB( self ):
        """
        Short cut for mask of all CB I{and} CA of GLY.

        @return: mask of all CB and CA of GLY
        @rtype: array
        """
        f = lambda a: a['name'] == 'CB' or\
            a['residue_name'] == 'GLY' and a['name'] == 'CA' 

        return self.maskF( f )


    def maskH2O( self ):
        """
        Short cut for mask of all atoms in residues named TIP3, HOH and  WAT

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskF( lambda a: a['residue_name'] in ['TIP3','HOH','WAT'])

    def maskSolvent( self ):
        """
        Short cut for mask of all atoms in residues named
        TIP3, HOH, WAT, Na+, Cl-

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskF( lambda a: a['residue_name'] in ['TIP3','HOH','WAT',
                                                          'Na+', 'Cl-'] )
    def maskHetatm( self ):
        """
        Short cut for mask of all HETATM 

        @return: N.array( 1 x N_atoms ) of 0||1
        @rtype: array
        """
        return self.maskF( lambda a: a['type'] == 'HETATM' )

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
        return self.maskF( lambda a, n=names: a['residue_name'].upper() in n )


    def indices( self, what ):
        """
        Get atom indices conforming condition

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
        if type( what ) == types.FunctionType:
            return N.nonzero( self.maskF( what) )

        if type( what ) == list or type( what ) == type( N.zeros(1) ):

            ## atom names
            if type( what[0] ) == str:
                return self.indices(
                    lambda a, allowed=what: a['name'] in allowed )

            if type( what[0] ) == int:
                ## mask
                if len( what ) == self.lenAtoms() and max( what ) < 2:
                    return N.nonzero( what )
                ## list of indices
                else:
                    return what

        ## single index
        if type( what ) == int:
            return N.array( [what], 'i' )

        raise PDBError("PDBModel.indices(): Could not interpret condition ")


    def mask( self, what, numpy=1 ):
        """
        Get atom mask.

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
            return self.maskF( what, numpy )

        if type( what ) == list or type( what ) == type( N.zeros(1) ):

            ## atom names
            if type( what[0] ) == str:
                return self.mask(
                    lambda a, allowed=what: a['name'] in allowed,   numpy)

            if type( what[0] ) == int:
                ## mask
                if len( what ) == self.lenAtoms() and max( what ) < 2:
                    return what
                ## list of indices
                else:
                    r = N.zeros( self.lenAtoms(),'i' )
                    N.put( r, what, 1 )
                    return r

        ## single index
        if type( what ) == int:
            return self.mask( [what] )

        raise PDBError, "PDBModel.mask(): Could not interpret condition "


    def atom2resMask( self, atomMask ):
        """
        Mask (0) residues for which all atoms are masked (0) in atomMask.

        @param atomMask: list/array of int, 1 x N_atoms
        @type  atomMask: list/array of int

        @return: 1 x N_residues (0||1 )
        @rtype: array of int
        """
        result = N.zeros( self.lenResidues(), 'i' )

        resMap = self.resMap()

        if len( resMap ) == 0:
            return result

        for res in range( 0, resMap[-1]+1 ):

            subMask = N.take( atomMask, N.nonzero( N.equal( resMap, res ) ) )
            result[res] = N.sum( subMask ) > 0

        return result


    def atom2resIndices( self, indices ):
        """
        Get list of indices of residue for which any atom is in indices.

        @param indices: list of atom indices
        @type  indices: list of int

        @return: indices of residues
        @rtype: list of int
        """
        result = []
        rI = self.resIndex()

        for i in range( len( rI ) ):

            a = rI[i]

            if i < len( rI )-1:
                e = rI[i+1]
            else:
                e = self.lenAtoms()

            if N.sum( N.less( indices, e ) * N.greater( indices, a-1 ) ):
                result += [i]

        return result


    def res2atomMask( self, resMask ):
        """
        convert residue mask to atom mask.

        @param resMask: list/array of int, 1 x N_residues
        @type  resMask: list/array of int

        @return: 1 x N_atoms
        @rtype: array of int
        """
        result = N.zeros( self.lenAtoms(), 'i')

        resMap = self.resMap()

        if len( resMap ) == 0:
            return result

        if len( resMask ) != resMap[-1]+1:
            raise PDBError, "invalid residue mask"

        N.put( result, self.res2atomIndices( N.nonzero( resMask ) ), 1 )

        return result


    def res2atomIndices( self, indices ):
        """
        Convert residue indices to atom indices.

        @param indices: list/array of residue indices
        @type  indices: list/array of int

        @return: list of atom indices
        @rtype: list of int
        """
        result = []
        resMap = self.resMap()

        if max( indices ) > resMap[-1] or min( indices ) < 0:
            raise PDBError, "invalid residue indices"

        for res in indices:
            result += N.nonzero( N.equal( resMap, res ) ).tolist()

        return result


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
            p = self.resProfile( p )

        isArray = isinstance( p, N.arraytype )

        resMap = self.resMap()

        r = [ p[ resMap[a] ] for a in range( len(resMap) ) ]

        if isArray:
            r = N.array( r )

        return r


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
        cm = self.chainMap( breaks=breaks )
        r = []

        for i in indices:
            chain = cm[i]
            if not chain in r:
                r += [ chain ]

        return r


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
        cm = self.chainMap( breaks=breaks )

        r = []
        for i in indices:
            r += N.nonzero( N.equal( cm, i ) ).tolist()

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


    def __takeAtomIndices( self, oldI, takeI ):
        """
        Translate atom positions so that they point to the same atoms after
        a call to N.take() (if they are kept at all).

        @param oldI: indices to translate
        @type  oldI: list of int
        @param takeI: indices for N
        @type  takeI: list of int

        @return: indices of current atoms after take, resorted
        @rtype: list of int
        """
        ind = N.take( range( self.lenAtoms() ), takeI )

        return [ i for i in range( len(ind) ) if ind[i] in oldI ]


    def __takeResIndices( self, oldI, takeI ):

        ind = N.take( range( self.lenResidues() ), takeI )
        return [ i for i in range( len(ind) ) if ind[i] in oldI ]


##     def __takeTerIndices( self, oldI, takeI, new_model ):
##         """
##         Pretty much a hack. Better would be perhaps having a profile
##         with chain index.
##         """
##         i_map = []

##         previous = 0
##         for i in oldI:
##             i_map +=  (i+1 - previous) * [i]
##             previous = i

##         i_map = N.take( i_map, takeI )

##         i_ter = []
##         for i in oldI:
##             eq = N.nonzero( N.equal( i_map, i ) )
##             if N.sum( eq ) > 0:
##                 i_ter += [ eq[-1] ]

##         ## make sure we got the last atom of a residue
##         atoms = new_model.getAtoms()

##         result = []
##         for i in i_ter:
##             res = atoms[i]['residue_name']
##             rn  = atoms[i]['residue_number']

##             while rn == atoms[i+1]['residue_number'] and \
##                   atoms[i+1]['residue_name'] == res:
##                 i += 1

##             r += [i]

##         return result


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

        r.setXyz( N.concatenate( ( self.getXyz(), m.getXyz() )  ) )
        r.setAtoms( self.getAtoms() + m.getAtoms() )
        r.setPdbCode( self.pdbCode )

        if None in self.rProfiles.values() or \
           None in self.aProfiles.values():
            self.update()

        r.rProfiles = self.rProfiles.concat( m.rProfiles )

        r.aProfiles = self.aProfiles.concat( m.aProfiles )

        r.__terAtoms = self.__terAtoms +\
                       (N.array( m.__terAtoms )+self.lenAtoms()).tolist()

        return r.concat( *models[1:] )


    def take( self, i, deepcopy=0 ):
        """
        All fields of the result model are shallow copies of this model's
        fields. I.e. removing or reordering of atoms does not affect the
        original model but changes to entries in the atoms dictionaries
        would also change the atom dictionaries of this model.

        take( atomIndices, [deepcopy=0] ) -> PDBModel / sub-class.

        @note: the position of TER records is translated to the new atoms.
        Chain boundaries can hence be lost (if the terminal atoms are not
        taken) or move into the middle of a  residue (if atoms are resorted).

        @param i: atomIndices, positions to take in the order to take
        @type  i: list/array of int
        @param deepcopy: deepcopy atom dictionaries (default 0)
        @type  deepcopy: 0||1

        @return: PDBModel / sub-class
        @rtype: PDBModel
        """
        ## needs current (old) atoms to create correct residue mask
        rIndices = self.atom2resIndices( i )

        r = self.__class__()

        r.source = self.source

        r.xyz = N.take( self.getXyz(), i, 0 )
        r.xyzChanged = self.xyzChanged or not (r.xyz == self.xyz)

        if deepcopy:
            r.atoms = [ copy.copy( a ) for a in self.getAtoms() ]
            r.atoms = N.take( r.atoms, i ).tolist()
        else:
            r.atoms = N.take( self.getAtoms(), i ).tolist()
        r.atomsChanged = self.atomsChanged or not (r.atoms == self.atoms)

        r.aProfiles = self.aProfiles.take( i )
        r.rProfiles = self.rProfiles.take( rIndices )

        r.info = copy.deepcopy( self.info )

        if self.__terAtoms:
            l = len( self.__terAtoms)
            r.__terAtoms = self.__takeAtomIndices( self.__terAtoms, i )

##             if len( r.__terAtoms ) != l:
##                 EHandler.warning('PDBModel lost position of %i TER records.'\
##                                  % (l - len( r.__terAtoms )) )

##             r.__terAtoms = self.__takeTerIndices( self.__terAtoms, i, r )

##           ## mark last atom of residues containing TER
##           terRes = self.atom2resIndices( self.__terAtoms )
##           terRes = self.__takeResIndices( terRes, self.atom2resIndices(i) )

##           rI = r.resIndex() + [ r.lenAtoms() ]
##           r.__terAtoms = [ rI[res+1] - 1 for res in terRes ]

        r.pdbCode = self.pdbCode
        r.fileName = self.fileName

        return r


    def keep( self, i ):
        """
        Replace atoms,coordinates,profiles of this(!) model with sub-set.
        (in-place version of N.take() )

        @param i: lst/array of int
        @type  i: 
        """
        if len(i)==self.lenAtoms() and max(i)<2:
            EHandler.warning('dont use PDBModel.keep() with mask.', trace=0) 

        r = self.take( i )

        self.atoms = r.atoms
        self.xyz = r.xyz
        self.atomsChanged = r.atomsChanged
        self.xyzChanged = r.xyzChanged

        self.aProfiles = r.aProfiles
        self.rProfiles = r.rProfiles

        self.info = r.info

        self.caMask = self.bbMask = self.heavyMask = None
        self.__resMap = self.__resIndex = None

        self.__terAtoms = r.__terAtoms


    def clone( self, deepcopy=0 ):
        """
        Clone PDBModel.

        @return: PDBModel / subclass, copy of this model,
                 see comments to Numeric.take()
        @rtype: PDBModel
        """
        return self.take( range( self.lenAtoms() ), deepcopy )


    def compress( self, mask, deepcopy=0 ):
        """
        Compress PDBmodel using mask.
        compress( mask, [deepcopy=0] ) -> PDBModel

        @param mask: N.array( 1 x N_atoms of 1 or 0 )
                     1 .. keep this atom
        @type  mask: array
        @param deepcopy: deepcopy atom dictionaries of this model and result
                         (default 0 )
        @type  deepcopy: 1||0

        @return: compressed PDBModel using mask
        @rtype: PDBModel
        """
        return self.take( N.nonzero( mask ), deepcopy )


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


    def takeChains( self, chainLst, deepcopy=0, breaks=0, maxDist=None ):
        """
        Get copy of this model with only the given chains.

        @param chainLst: list of chains
        @type  chainLst: list of int
        @param deepcopy: deepcopy atom dictionaries (default 0)
        @type  deepcopy: 0||1
        @param breaks: split chains at chain breaks (default 0)
        @type  breaks: 0|1
        @param maxDist: (if breaks=1) chain break threshold in Angstrom
        @type  maxDist: float

        @return: PDBModel (by default, all fields are shallow copies,
                           see Numeric.take() )
        @rtype: PDBModel
        """
        chainMap = self.chainMap( breaks=breaks, maxDist=None )

        atomMask = map( lambda i, allowed=chainLst: i in allowed, chainMap)

        return self.compress( atomMask, deepcopy )


    def addChainFromSegid(self, verbose=1):
        """
        Takes the last letter of the segment ID and adds it as chain ID.
        """
        for a in self.getAtoms():

            try:
                a['chain_id'] = a['segment_id'][-1]
            except:
                if verbose:
                    EHandler.warning("addChainId(): Problem with atom "+str(a))

        self.atomsChanged = 1


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
        atoms = self.getAtoms()

        old_chains = []
        if keep_old:
            old_chains = [ atoms[i]['chain_id'] for i in self.resIndex(breaks)]
            old_chains = mathUtils.nonredundant( old_chains )
            if '' in old_chains: old_chains.remove('')

        letters = string.uppercase
        if first_id:
            letters = letters[ letters.index( first_id ): ]
        letters = mathUtils.difference( letters, old_chains )

        chainMap = self.chainMap( breaks=breaks )

        for i in range( len( atoms ) ):

            if not (keep_old and atoms[i]['chain_id'] in old_chains):
                atoms[i]['chain_id'] = letters[ chainMap[i] ]

        self.atomsChanged = 1


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


    def lenAtoms( self ):
        """
        Number of atoms in model.

        @return: number of atoms
        @rtype: int
        """
        if not self.xyz is None:
            return len( self.xyz )

        if not self.atoms is None:
            return len( self.atoms )

        if self.source is None:
            return 0

        return len( self.getXyz() )


    def lenResidues( self ):
        """
        Number of resudies in model.

        @return: total number of residues
        @rtype: int
        """
        try:
            return self.resMap()[-1] + 1
        except IndexError:  ## empty residue map
            return 0


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
            return self.chainMap( breaks=breaks, maxDist=maxDist )[-1] + 1
        except IndexError:  ## empty residue map
            return 0


    def resList( self, mask=None ):
        """
        Return list of lists of atom dictionaries per residue,
        which allows to iterate over residues and atoms of residues.

        @param mask: [ 0||1 x N_atoms ] atom mask to apply BEFORE
        @type  mask: 

        @return: A list of dictionaries::
        [ [ {'name':'N', 'residue_name':'LEU', ..},          
            {'name':'CA','residue_name':'LEU', ..} ],        
          [ {'name':'CA', 'residue_name':'GLY', ..}, .. ] ]      
        @rtype: list of dictionaries    
        """
        ri = N.concatenate( (self.resIndex( mask=mask ), [self.lenAtoms()] ) )
        resLen = len( ri ) - 1
        atoms = self.getAtoms()
        if mask:
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
        atoms = self.getAtoms()
        xyz = self.getXyz()

        result = []
        for res in  range( resLen ):
            m = PDBModel()
            a = ri[res]
            if res == resLen - 1:
                e = self.lenAtoms()
            else:
                e = ri[res + 1]
            m.setAtoms( atoms[a:e] )
            m.setXyz(     xyz[a:e] )
            result += [m]
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
        result = []
        for a in self.getAtoms():
            result.append( a['residue_number'] )

        ## by default take all atoms
        mask = mask or N.ones( len(self.getAtoms() ) , 'i' )

        ## apply mask to this list
        return N.compress( mask, N.array(result, 'i') )


    def __calcResMap( self, mask=None ):
        """
        Create a map of residue residue for atoms in model.

        @param mask: atom mask
        @type  mask: list of int (1||0)

        @return: array [00011111122223333..], residue index for each
                 unmasked atom
        @rtype:  list of int        
        """
        ## By default consider all atoms
        mask = mask or N.ones( self.lenAtoms() )

        result = []

        lastResNumber = -100
        lastResName   = ''
        index = -1
        lastAlt = 'x'

        ## create residue numbering for selected atoms
        for a in N.compress( mask, self.getAtoms() ):

            if lastResNumber <> a['residue_number'] or \
               lastResName   <> a['residue_name'] or \
               lastAlt <> a.get('insertion_code', None ):

                ## start of new residue
                lastResNumber = a['residue_number']
                lastResName = a['residue_name']
                lastAlt = a.get( 'insertion_code', None )
                index += 1

            result.append( index )

        return N.array(result, 'i')        


    def resMap(  self, mask=None, force=0, cache=1 ):
        """
        Get list to map from any atom to a continuous residue numbering
        (starting with 0). A new residue is assumed to start whenever the
        'residue_number' or the 'residue_name' record changes between 2
        atoms. The mask is applied BEFORE looking for residue borders,
        i.e. it can change the residue numbering.

        See L{resList()} for an example of how to use the residue map.

        @param mask: [0000011111001111...] include atom: yes or no
                     len(atom_mask) == number of atoms in self.atoms/self.xyz
        @type  mask: list of int (1||0)
        @param force: recalculate map even if cached one is available (def 0)
        @type  force: 0||1
        @param cache: cache new map (def 1)
        @type  cache: 0||1

        @return: array [00011111122223333..], residue index for each
                 unmasked atom
        @rtype:  list of int
        """
        if self.__resMap and not force and mask == None:
            return self.__resMap

        result = self.__calcResMap( mask )

        if cache and mask == None:
            self.__resMap = result

        return result


    def resIndex( self, mask=None, force=0, cache=1 ):
        """
        Get the position of the each residue's first atom. The result is by
        default cached. That's not really necessary - The calculation is fast.

        @param force: re-calculate even if cached result is available (def 0)
        @type  force: 1||0
        @param cache: cache the result (def 1)
        @type  cache: 1||0
        @param mask: atom mask to apply before (i.e. result indices refer to
                     compressed model)
        @type  mask: list of int (1||0)

        @return: index of the first atom of each residue
        @rtype: list of int
        """

        if self.__resIndex and not force and mask == None:
            return self.__resIndex

        cache = cache and mask == None

        m = self.resMap()
        if mask:
            m = N.compress( mask, m )

        r = [ i for i in range(len(m)) if i==0 or m[i] != m[i-1] ]

        if cache:
            self__resIndex = r

        return r


    def resEndIndex( self, mask=None ):
        """
        Get the position of the each residue's last atom.

        @param mask: atom mask to apply before (i.e. result indices refer to
                     compressed model)
        @type  mask: list of int (1||0)

        @return: index of the last atom of each residue
        @rtype: list of int
        """

        m = self.resMap()
        if mask:
            m = N.compress( mask, m )

        l = len(m)
        return [ i-1 for i in range(1, l) if i==l or m[i] != m[i-1] ]


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
        result = []

        i = 0
        lastChain = -1
        lastResidue = -100
        lastChainID = None
        lastSegID = None

        ter_atoms = self.__terAtoms
        if breaks:
            ter_atoms = N.concatenate( (ter_atoms,
                                        self.chainBreaks( maxDist=maxDist )) )

        for a in self.getAtoms() :

            if a.get('chain_id',None) <> lastChainID or \
               a.get('segment_id', None) <> lastSegID or \
               a['residue_number'] < lastResidue or \
               i-1 in ter_atoms:

                lastChain += 1

            result.append( lastChain )

            lastResidue = a['residue_number']
            lastChainID = a.get('chain_id',None)
            lastSegID = a.get( 'segment_id', None )
            i += 1

        return N.array( result, 'i' )


    def chainIndex( self, breaks=0, maxDist=None ):
        """
        Get indices of first atom of each chain.

        @param breaks: split chains at chain breaks (def 0)
        @type  breaks: 1||0
        @param maxDist: (if breaks=1) chain break threshold in Angstrom
        @type  maxDist: float

        @return: array (1 x N_chains) of int
        @rtype: list of int
        """
        ci = self.chainMap( breaks=breaks, maxDist=maxDist )

        result = []

        try:
            for chain in range( 0, ci[-1]+1 ):
                result.append( N.nonzero( N.equal( ci, chain ))[0] )

        except IndexError:   ## empty chainMap -> empty model?
            pass

        return N.array( result, 'i' )


    def chainBreaks( self, breaks_only=1, maxDist=None ):
        """
        Identify discontinuities in the molecule's backbone.

        @param breaks_only: don't report ends of regular chains (def 1)
        @type  breaks_only: 1|0
        @param maxDist: maximal distance between consequtive residues
                        [ None ] .. defaults to twice the average distance
        @type  maxDist: float

        @return: atom indices of last atom before a probable chain break
        @rtype: list of int
        """
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

        if not breaks:
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

        return N.array( r )


    def removeRes( self, resname ):
        """
        Remove all atoms with a certain residue name.

        @param resname: name of residue to be removed
        @type  resname: str OR list of str
        """
        resname = t.toList( resname )
        self.remove( lambda a, res=resname: a['residue_name'] in res )


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

        mask_fit = mask_fit or mask

        if fit:

            fx, fy = x, y

            if mask_fit:
                fx = N.compress( mask_fit, x, 0 )
                fy = N.compress( mask_fit, y, 0 )
            ## find transformation for best match
            r, t = rmsFit.match( fx, fy, n_iterations=n_it )[0]

            ## transform coordinates
            y = N.dot(y, N.transpose(r)) + t

        if mask:
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

        if mask:
            x = N.compress( mask, x, 0 )
            y = N.compress( mask, y, 0 )
            outlier_mask = N.zeros( N.sum(mask)  )

        else:
            outlier_mask = N.zeros( len(self) )

        r, iter_trace = rmsFit.match( x, y, n_iterations=n_it, z=z,
                                      eps_rmsd=eps_rmsd, eps_stdv=eps_stdv)

        N.put( outlier_mask, iter_trace[-1][-1], 1 )

        if n_it != 1:
            self.setAtomProfile( profname, outlier_mask, mask=mask,
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
        if mask != None:
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
        mask = mask or N.ones( len(self) )

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
        if mask == None:
            return N.average( self.getXyz() )

        return N.average( N.compress( mask, self.getXyz(), axis=0 ) )


    def centerOfMass( self ):
        """
        Center of mass of PDBModel.

        @return: array('f')
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
            M = [ molUtils.atomMasses[a['element']]for a in self.getAtoms() ]

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
        if mask == None:
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

        return N.array( result, 'f' )


    def argsort( self, cmpfunc=None ):
        """
        Prepare sorting atoms within residues according to comparison function.

        @param cmpfunc: function( self.atoms[i], self.atoms[j] ) -> -1, 0, +1
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


    def sort( self, sortArg=None, deepcopy=0  ):
        """
        Apply a given sort list to the atoms of this model.

        @param sortArg: comparison function
        @type  sortArg: function
        @param deepcopy: deepcopy atom dictionaries (default 0)
        @type  deepcopy: 0||1

        @return: copy of this model with re-sorted atoms (see Numeric.take() )
        @rtype: PDBModel
        """
        sortArg = sortArg or self.argsort()

        terA = self.__terAtoms

        r = self.take( sortArg, deepcopy )

        r.__terAtoms = terA

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

        return [ a['name'] for a in self.getAtoms()[i:j] ]


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
        atoms = self.getAtoms()
        if mode == 0:
            f_test = self.__testDict_and
        else:
            f_test = self.__testDict_or

        for k in kw:
            kw[ k ] = t.toList( kw[ k ] )

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
        Coordinates are not checked, profiles are not checked.

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


    def equalAtoms( self, ref ):
        """
        Apply to SORTED models without HETATOMS. Coordinates are not checked.

        @note: in some rare cases m1.equalAtoms( m2 ) gives a different result
               than m2.equalAtoms( m1 ). This is due to the used
               SequenceMatcher class.      
        @todo: option to make sure atoms are also in same order

        @param ref: reference PDBModel
        @type  ref: PDBModel

        @return: (mask, mask_ref), two atom masks for all equal (1) atoms
                  in models
        @rtype: (array, array)
        """
        ## compare sequences
        seqMask, seqMask_ref = match2seq.compareModels(self, ref)

        ## get residue mask on atom level
        mask = self.res2atomMask(seqMask)
        mask_ref = ref.res2atomMask(seqMask_ref)

        ## get list of matching RESIDUES
        equal = N.nonzero(seqMask)
        equal_ref = N.nonzero(seqMask_ref)

        ## check that all atoms are equal in matching residues
        for i in range(0, len(equal)):

            ## atom name lists for current residue
            aa = self.atomNames(equal[i],equal[i])
            aa_ref = ref.atomNames(equal_ref[i],equal_ref[i])

            ## starting atom of current residue
            ind = self.resIndex()[ equal[i] ]
            ind_ref = ref.resIndex()[ equal_ref[i] ]

            ## check that the atoms are the same
            for j in range( len(aa) ):

                try:
                    mask[ind] = aa[j] == aa_ref[j]
                    ind += 1

                except IndexError:
                    mask[ind] = 0
                    ind += 1

            for j in range( len(aa_ref) ):

                try:
                    mask_ref[ind_ref] = aa_ref[j] == aa[j]
                    ind_ref += 1
                except IndexError:
                    mask_ref[ind_ref] = 0
                    ind_ref += 1


        return mask, mask_ref


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
        l_0 = len( m )
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
        
        @return: coordinates of center of mass
        @rtype:  array
        """
        ## loading output file from X-plor
        if local: print 'Loading pdb file ..'
        m = PDBModel( t.testRoot()+'/rec/1A2P.pdb')#"/com_wet/1BGS.pdb")

        ## remove solvent
        if local: print "removing waters.."
        m.removeRes(['TIP3', 'HOH'])

        ## X-plor doesn't write chainIds, so during the simulation
        ## we store them in the last letter of the segId. Here we
        ## restore the chainId.
        m.addChainFromSegid()

        ## start positions of all chains
        chainIdx = m.chainIndex().tolist()

        ## print some chain info
        if local:
            print 'The molecule consists of %i chains'% m.lenChains()
            print '\tChainId \tFirst atom'
            for i in chainIdx:
                print '\t%s \t\t%i'%(m.atoms[i]['chain_id'], int(i))

        ## iterate over all chains
        for c in range( 0, len( chainIdx ) ):

            if local:
                print "chain ", c, " starts with ", 
                print m.atoms[ chainIdx[c] ]['residue_name'],

                print " and has sequence: "

            ## mask out atoms of all other chains
            chainMask  = N.equal( m.chainMap( breaks=1 ), c )
            if local: print m.sequence( chainMask )

        ## test sorting
        if local: print "sorting atoms alphabetically..."
        sort = m.argsort()
        m2 = m.sort( sort )

        if local:
            globals().update( locals() )

        return N.sum( m2.centerOfMass() )


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: coordinates of center of mass
        @rtype:  array
        """
        return N.sum( N.array([ 29.13901386,  46.70021977,  38.34113311]) )
    
        

if __name__ == '__main__':

    test = Test()

    assert abs( test.run( local=1 ) - test.expected_result() ) < 1e-8

## last line count:: 3229
