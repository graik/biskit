##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

import tools as t
import molUtils
import mathUtils
import match2seq
import rmsFit
from LocalPath import LocalPath
from Errors import BiskitError
from Biskit import EHandler
from ProfileCollection import ProfileCollection, ProfileError

import Numeric as N
import MLab
from Scientific.IO.PDB import PDBFile
import os
import copy
import time
import string
import types

from multiarray import arraytype

## strange bugfix #2 (unpickling of old PDBModels)
## from Numeric import array_constructor


class PDBProfiles( ProfileCollection ):
    def version( self ):
        return ProfileCollection.version(self) + '; PDBModel $Revision$'

class PDBError(BiskitError):
    pass

class PDBModel:
    """
    Store and manipulate coordinates and atom infos stemming from a
    PDB file. Coordinates are stored in the Numeric array 'xyz'; the
    additional atom infos from the PDB are stored in the list of dictionaries
    'atoms'.
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

    @todo
    - clean up, re-organize parsing
    - clean up, rename returnPdb()
    - use Profiles instead of space consuming atom dictionaries ?
    """

    def __init__( self, source=None, pdbCode=None, noxyz=0, skipRes=None ):
        """
        * PDBModel() creates an empty Model to which coordinates (field xyz)
          and PDB infos (field atoms) have still to be added.
        * PDBModel( file_name ) creates a complete model with coordinates
          and PDB infos taken from file_name (pdb, pdb.gz, pickled PDBModel)
        * PDBModel( PDBModel ) creates a copy of the given model
        * PDBModel( PDBModel, noxyz=1 ) creates a copy without coordinates
        
        source  - str, file name of pdb/pdb.gz file OR pickled PDBModel OR
                - PDBModel, template structure to copy atoms/xyz field from
        pdbCode - str, PDB code, is extracted from file name otherwise
        noxyz   - 0 (default) || 1, create without coordinates
        
        !! raise PDBError, if file can't be read
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

        if source <> None:

            self.update( skipRes=skipRes, lookHarder=1 )

        if noxyz:
            ## discard coordinates, even when read from PDB file
            self.xyz = None

        ## monitor changes of coordinates and atoms
        self.xyzChanged = 0
        self.atomsChanged = 0

        self.forcePickle = 0

        ## version as of creation of this object
        self.initVersion = self.version()

        ## to collect further informations
        self.info = { 'date':t.dateSortString() }


    def version( self ):
        return 'PDBModel $Revision$'


    def __getstate__(self):
        """called before pickling the object."""
        self.slim()
        self.forcePickle = 0
        return self.__dict__

    def __setstate__(self, state ):
        """called for unpickling the object."""
        self.__dict__ = state
        ## backwards compability
        self.__defaults() 

    def __len__(self):
        return self.lenAtoms()

    def __defaults(self ):
        """backwards compatibility to earlier pickled models"""
        
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
        skipRes    - lst of str, names of residues to skip if updating from PDB
        lookHarder - 0|1, 0(default): update only existing profiles
        force      - 0|1, ignore invalid source (0) or report error (1)
        
        !! raise PDBError if file can't be unpickled or read
        """
        source = self.validSource()

        if source is None and force:
            raise PDBError( str(self.source) + ' is not a valid source.')

        if source is None:
            return

        isPdb = type( source ) == str and\
                (source[-4:].upper() == '.PDB' or
                 source[-7:].upper() == '.PDB.GZ')

        try:
            ## updating from PDBModel (pickled or not)
            if not isPdb:

                if type( source ) == str:
                    s = t.Load( source )
                else:
                    s = source

                self.fileName = self.fileName or s.fileName

                self.pdbCode = self.pdbCode or s.pdbCode

                self.atoms = self.atoms or s.getAtoms()

                self.xyz = self.xyz or s.getXyz()

                self.__terAtoms = getattr(self, '_PDBModel__terAtoms',[]) or \
                                  getattr(s,'_PDBModel__terAtoms',[])

                self.rProfiles.updateMissing( s.rProfiles,
                                              copyMissing=lookHarder)
                self.aProfiles.updateMissing( s.aProfiles,
                                              copyMissing=lookHarder)

                if skipRes:
                    self.removeRes( skipRes )
                
                return

        except ProfileError, why:
            EHandler.warning("Cannot read/update profiles from source.")
        except Exception, why:
            EHandler.warning("Cannot unpickle source model from %s, "\
                   % str(source) +
                   "Reason:\n" + str(why) + "\n" +
                   "I try opening file as ordinary PDB -- fingers crossed.",
                    error=0 ) 
            ## try again as PDB
            isPdb = 1

        try:
            ## atoms and/or coordinates need to be updated from PDB
            if (self.atoms==None or self.xyz==None) and isPdb:

                atoms, xyz = self.__collectAll( source, skipRes )

                self.atoms = self.atoms or atoms

                self.xyz = self.xyz or xyz

                self.__terAtoms = self.__pdbTer()
                
                self.fileName = self.fileName or source

                self.pdbCode = self.pdbCode or \
                               os.path.basename( self.fileName )[:4]
        except:
            raise PDBError('Cannot read ' + str(source) + ' as PDB\n'\
                           '\ERROR: ' + t.lastError() )


    def __pdbTer( self, rmOld=0 ):
        """
        -> list of atom indices that are followed by a TER record (marked
        with 'after_ter' flag of the next atom by __collectAll).
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
        xyz - Numpy array ( 3 x N_atoms ) of float
        -> N.array( 3 x N_atoms ) or None, old coordinates
        """
        old = self.xyz
        self.xyz = xyz

        self.xyzChanged = self.xyzChanged or \
                          not mathUtils.arrayEqual(self.xyz,old )
        return old


    def setAtoms( self, atoms ):
        """
        Replace atoms list of dictionaries.
        atoms - list of dictionaries as returned by PDBFile.readLine()
                (nearly, for differences see __collectAll() )
        self.__terAtoms is re-created from the 'after_ter' records in atoms
        -> [ dict ], old atom dictionaries
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
        source - LocalPath OR PDBModel OR str
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
        mask - [int] OR array of 1||0, atom mask
        -> N.array( 3 x N_atoms, 'f' ) 
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
        if necessary. See also __collectAll().
        mask - [int] OR array of 1||0, atom mask
        -> list of  dictionaries from PDBFile.readline()
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
        name     - str, name to access profile
        prof     - list/array of values
        mask     - list/array 1 x N_residues of 0|1, N.sum(mask)==len(prof)
        default  - any, value for residues masked.
        asarray  - 0|1|2, store as list (0), as array (2) or store
                   numbers as array but everything else as list (1) [1]
        comment  - str, goes into aProfiles_info[name]['comment']
        moreInfo - additional key-value pairs for aProfiles_info[name]
        !! raise ProfileError, if length of prof != N_residues
        """
        self.rProfiles.set( name, prof, mask=mask, default=default,
                            asarray=asarray, comment=comment, **moreInfo )

    def setAtomProfile( self, name, prof, mask=None, default=None, asarray=1,
                        comment=None, **moreInfo):
        """
        Add/override atom-based profile.
        name     - str, name to access profile
        prof     - list/array of values
        mask     - list/array 1 x N_residues of 0 || 1, profile items
        default  - any, value masked residues [None for lists, 0 for arrays]
        asarray  - 0|1|2, store as list (0), as array (2) or store
                   numbers as array but everything else as list (1) [1]
        comment  - str, goes into aProfiles_info[name]['comment']
        moreInfo - additional key-value pairs for aProfiles_info[name]
        !! raise ProfileError, if length of prof != N_atoms
        """
        self.aProfiles.set( name, prof, mask=mask, default=default,
                            asarray=asarray, comment=comment, **moreInfo )


    def resProfile( self, name, default=None ):
        """resProfile( profile_name ) -> array 1 x N_res with residue values
        !! raise ProfileError, if rProfiles contains no entry for |name|
        """
        return self.rProfiles.get( name, default=default )


    def atomProfile(self, name, default=None ):
        """atomProfile( profile_name ) -> array 1 x N_atoms with atom values
        !! raise ProfileError, if aProfiles contains no entry for |name|
        """
        return self.aProfiles.get( name, default=default )


    def profile( self, name, default=None, lookHarder=0 ):
        """
        profile( name, lookHarder=0) -> atom or residue profile
        default - default result if no profile is found,
                  if None, raise exception
        lookHarder - 0|1, update from source before reporting missing profile
        !! raise ProfileError, if neither atom- nor rProfiles contains |name|
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
        lookHarder - 0|1, update from source before reporting missing profile
        Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1||0
        !! raise ProfileError, if neither atom- nor rProfiles contains |name|
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
        """Add/Override infos about a given profile
        e.g. setProfileInfo('relASA', comment='new', params={'bin':'whatif'})
        !! raise ProfileError, if neither atom- nor rProfiles contains |name|
        """
        d = self.profileInfo( name )
        for key, value in args.items():
            d[key] = value


    def removeProfile( self, *names ):
        """
        removeProfile( str_name [,name2, name3] ) -> 1|0,
        remove res or atom profile(s)
        -> int, 1 if at least 1 profile has been deleted, 0,
           if none has been found
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
        -> (1||0, 1||0), (1,1)..both xyz and atoms field have been changed
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
        -> (1||0, 1||0), (1,1)..both xyz and atoms field have been changed
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
        -> 1, if profile |pname| can currently not be reconstructed from a
        source on disc.
        !! raise ProfileError if there is no atom or res profile with pname
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
        AUTOMATICALLY CALLED BEFORE PICKLING
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
        AUTOMATICALLY CALLED BEFORE PICKLING
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
        """ -> str or PDBModel, None if this model has no valid source."""
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
        -> str, file name of pickled source or PDB file
        !! PDBError if there is no valid source
        """
        s = self.validSource()

        if s is None:
            raise PDBError('no valid source')

        if type( s ) == str:
            return s

        return self.source.sourceFile()


    def disconnect( self ):
        """
        Disconnect this model from its source (if any). Note: If this model
        has an (in-memory) PDBModel instance as source, the entries of
        'atoms' could still reference the same dictionaries.
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
        """-> string"""
        return self.pdbCode

    def setPdbCode(self, code ):
        """ code - str, new pdb code """
        self.pdbCode = code
        

    def __firstLetter( self, aName ):
        """ -> first letter (i.e. not a number) from a string. """
        try:
            i = int( aName[0] )
            return self.__firstLetter( aName[1:] )
        except:
            return  aName[0]


    def __collectAll( self, fname, skipRes=None ):
        """
        Parse ATOM/HETATM lines from PDB. Collect coordinates plus
        dictionaries with the other pdb records of each atom.
        REMARK, HEADER, etc. lines are ignored.
        fname - name of pdb file
        
        Some changes are made to the dictionary from PDBFile.readline():
            - the 'position' entry (with the coordinates) is removed
            - leading and trailing spaces are removed from 'name' ..
            - .. but a 'name_original' entry keeps the old name with spaces
            - a 'type' entry is added. Its value is 'ATOM' or 'HETATM'
            - a 'after_ter' entry is added. Its value is 1, if atom is
              preceeded by a 'TER' line, otherwise 0
            - empty 'element' entries are filled with the first non-number
              letter from the atom 'name'

        -> (list of  dictionaries from PDBFile.readline() , xyz array N x 3)
        """
        items = []
        xyz   = []

        f = PDBFile( fname )
            
        try:
            line, i = ('',''), 0
            
            while line[0] <> 'END' and line[0] <> 'ENDMDL':

                i += 1
                try:
                    line = f.readLine()
                except ValueError, what:
                    t.errWriteln('Warning: Error parsing line %i of %s' %
                               (i, t.stripFilename( fname )) )
                    t.errWriteln('\tError: '+str(what) )
                    continue

                ## preserve position of TER records
                newChain = line[0] == 'TER'
                if newChain:
                    line = f.readLine()
                    
                if (line[0] in ['ATOM','HETATM'] ):

                    a = line[1]

                    if skipRes and a['residue_name'] in skipRes:
                        continue

                    a['name_original'] = a['name']
                    a['name'] = a['name'].strip()

                    a['type'] = line[0]
                    if newChain: a['after_ter'] = 1

                    if a['element'] == '':
                        a['element'] = self.__firstLetter( a['name'] )

                    if a['position'].is_vector:
                        lst = [ a['position'][0],
                                a['position'][1],
                                a['position'][2]]
                        xyz.append( lst )
                    else:
                        xyz.append( a['position'] )

                    del a['position']

                    items += [ a ]
            
        except:
            raise PDBError("Error parsing file "+fname+": " + t.lastError() )

        try:
            f.close()
        except:
            pass

        if len( xyz ) == 0:
            raise PDBError("Error parsing file "+fname+": "+
                            "Couldn't find any atoms.")
            
        return items, N.array( xyz, 'f' )
    

    def sequence(self, mask=None, xtable=molUtils.xxDic ):
        """
        mask   - atom mask, to apply before (default None)
        xtable - dict {str:str}, additional residue:single_letter mapping
                 for non-standard residues (default molUtils.xxDic)
        -> str, 1-letter-code AA sequence (based on first atom of each res).
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
        change - 1|0, change this model's atoms directly      [1]
        aatm   - 1|0, use, for example, HG23 instead of 3HG2  [1]
        -> [ {..} ], list of atom dictionaries
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
        fname - name of new file
        ter   - 0, don't write any TER statements
              - 1, restore original TER statements (doesn't work, if preceeding
                atom has been deleted) [default]
              - 2, put TER between all detected chains
              - 3, as 2 but also detect and split discontinuous chains
        amber - 1|0, amber formatted atom names (ter=3, left=1, wrap=0)   [0]
        original  - 1|0, revert atom names to the ones parsed in from PDB [0]
        left      - 1|0, left-align atom names (as in amber pdbs)  [0]
        wrap      - 1|0, write e.g. 'NH12' as '2NH1'               [0]
        headlines - [ ( str, dict or str) ], list of record / data tuples
                    e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
        taillines - same as headlines but appended at the end of file
        """
        try:
            f = PDBFile( fname, mode='w' )
            
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
        dictionaries in self.atoms. Older version of writePdb that returns
        a list of PDB lines instead of writing to a file.
        out   - stdout or None, if None a list is returned
        ter   - 0, don't write any TER statements
              - 1, restore original TER statements (doesn't work, if preceeding
                atom has been deleted)
              - 2, put TER between all detected chains
        headlines - [ ( str, dict or str) ], list of record / data tuples
                    e.g. [ ('SEQRES', '  1 A 22  ALA GLY ALA'), ]
        taillines - same as headlines but appended at the end of file
        -> [ str ], lines of a PDB file
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
        path - str OR LocalPath instance, target file name
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
        atomF - function( dict_from_PDBFile.readline()  ) -> true || false
                (Condition)
        numpy - 1(default)||0, convert result to Numpy array of int
        -> Numpy N.array( [0,1,1,0,0,0,1,0,..], 'i') or list
        """
        try:
            result = map( atomFunction, self.getAtoms() )
        except:

            ## fall-back solution: assign 0 to all entries that raise
            ## exception
            t.errWriteln("Warning mask(): Error while mapping funtion "+
                       "to all atoms: " + t.lastError() )
            result = []

            for a in self.getAtoms():
                try:
                   result.append( atomFunction( a ) )
                ## put 0 if something goes wrong
                except :
                    t.errWriteln("Warning mask(): Error while save-mapping "+
                               t.lastError() )
                    result.append(0)

        if numpy:
            return N.array( result, 'i' )
        return result


    def maskCA( self, force=0 ):
        """
        Short cut for mask of all CA atoms.
        force - 0 || 1, force calculation even if cached mask is available
        -> N.array( 1 x N_atoms ) of 0 || 1
        """
        if self.caMask == None or force:
            self.caMask = self.mask( lambda a: a['name'] == 'CA' )

        return self.caMask


    def maskBB( self, force=0 ):
        """
        Short cut for mask of all backbone atoms.
        force - 0 || 1, force calculation even if cached mask is available
        -> N.array( 1 x N_atoms ) of 0 || 1
        """
        if self.bbMask == None or force:
            self.bbMask = self.mask( 
                lambda a: a['name'] in ['CA', 'C', 'N', 'O', 'H','OXT'] )

        return self.bbMask


    def maskHeavy( self, force=0 ):
        """
        Short cut for mask of all heavy atoms. ('element' <> H)
        force - 0 || 1, force calculation even if cached mask is available
        -> N.array( 1 x N_atoms ) of 0 || 1
        """
        if self.heavyMask == None or force:
            self.heavyMask = self.mask( lambda a: a.get('element','') <> 'H' )

        return self.heavyMask

    def maskH( self ):
        return N.logical_not( self.maskHeavy() )


    def maskCB( self ):
        """-> mask of all CB and CA of GLY """
        f = lambda a: a['name'] == 'CB' or\
            a['residue_name'] == 'GLY' and a['name'] == 'CA' 

        return self.mask( f )
    

    def maskH2O( self ):
        """-> mask of all atoms in residues named TIP3, HOH, WAT """
        return self.mask( lambda a: a['residue_name'] in ['TIP3','HOH','WAT'] )

    def maskSolvent( self ):
        """-> all atoms in residues named TIP3, HOH, WAT, Na+, Cl-"""
        return self.mask( lambda a: a['residue_name'] in ['TIP3','HOH','WAT',
                                                          'Na+', 'Cl-'] )
    def maskHetatm( self ):
        """-> mask of all HETATMs """
        return self.maskF( lambda a: a['type'] == 'HETATM' )

    def maskProtein( self, standard=0 ):
        """
        standard - 0|1, only standard residue names (not CYX, NME,..) [0]
        -> mask of all protein atoms (based on residue name)
        """
        d = molUtils.aaDic
        if standard:
            d = molUtils.aaDicStandard
            
        names = map( string.upper, d.keys() )
        return self.maskF( lambda a, n=names: a['residue_name'].upper() in n )
    

    def indices( self, what ):
        """
        Get atom indices conforming condition
        what - function applied to each atom entry,
                e.g. lambda a: a['residue_name']=='GLY'
             - list of str, allowed atom names
             - list of int, allowed atom indices OR mask with only 1 and 0
             - int, single allowed atom index
        -> Numeric array, N_atoms x 1 (0 || 1 )
        !! raise PDBError, if what is neither of above
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
        what - function applied to each atom entry,
                e.g. lambda a: a['residue_name']=='GLY'
             - list of str, allowed atom names
             - list of int, allowed atom indices OR mask with only 1 and 0
             - int, single allowed atom index
        -> Numeric array, N_atoms x 1 (0 || 1 )
        !! raise PDBError, if what is neither of above
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
        atomMask - list/array of int, 1 x N_atoms
        -> array of int, 1 x N_residues
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
        -> [int], indices of residues for which any atom is in indices
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
        resMask - list/array of int, 1 x N_residues
        -> array of int, 1 x N_atoms
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
        indices - list/array of int
        -> [ int ]
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
        p - str, name of existing residue profile OR...
            [ any ], list of lenResidues() length
        -> [ any ] OR array, atom profile 
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
        Convert atom indices to chain indices. Each chain is only returned
        once.
        indices - [ int ]
        breaks  - 0|1, look for chain breaks in backbone coordinates [0]
        -> [ int ], chains any atom which is in indices
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
        indices - [ int ]
        breaks  - 0|1, look for chain breaks in backbone coordinates [0]
        -> [ int ], all atoms belonging to the given chains
        """
        cm = self.chainMap( breaks=breaks )

        r = []
        for i in indices:
            r += N.nonzero( N.equal( cm, i ) ).tolist()

        return r


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2mask( str_profname, [cutoff_min, cutoff_max=None])
        -> mask len( profile(profName) ) x 1|0
        !! ProfileError if no profile is found with name profName
        """
        p = self.profile( profName )

        cutoff_min = cutoff_min or min( p ) - 1
        cutoff_max = cutoff_max or max( p ) + 1

        return N.greater( p, cutoff_min ) * N.less( p, cutoff_max )


    def profile2atomMask( self, profName, cutoff_min=None, cutoff_max=None ):
        """
        profile2atomMask( str_profname, [cutoff_min, cutoff_max=None])
        Same as profile2mask, but converts residue mask to atom mask.
        -> mask N_atoms x 1|0
       !! ProfileError if no profile is found with name profName
        """
        r = self.profile2mask( profName, cutoff_min, cutoff_max )

        if len( r ) == self.lenResidues():
            r = self.res2atomMask( r )

        return r


    def __takeAtomIndices( self, oldI, takeI ):
        """
        Translate atom positions so that they point to the same atoms after
        a call to N.take() (if they are kept at all).
        oldI - list of int, indices to translate
        takeI- list of int, indices for N.take()
        -> list of int, indices of current atoms after take, resorted
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
        model0.concat( model1 [, model2, ..]) -> single PDBModel.
        Concatenate atoms, coordinates and profiles. source and fileName
        are lost, so are profiles that are not available in all models.
        Note: info records of given models are lost.
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
        take( atomIndices, [deepcopy=0] ) -> PDBModel / sub-class.
        All fields of the result model are shallow copies of this model's
        fields. I.e. removing or reordering of atoms does not affect the
        original model but changes to entries in the atoms dictionaries
        would also change the atom dictionaries of this model.
        Note: the position of TER records is translated to the new atoms.
        Chain boundaries can hence be lost (if the terminal atoms are not
        taken) or move into the middle of a  residue (if atoms are resorted).

        atomIndices - list/array of int, positions to take in the order to take
        deepcopy    - 0||1, deepcopy atom dictionaries (0)
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
        i - lst/array of int
        """
        if len(i)==self.lenAtoms() and max(i)<2:
            t.errWriteln('WARNING: dont use PDBModel.keep() with mask.') 
        
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
        -> PDBModel / subclass, copy of this model, see comments to N.take()
        """
        return self.take( range( self.lenAtoms() ), deepcopy )


    def compress( self, mask, deepcopy=0 ):
        """
        compress( mask, [deepcopy=0] ) -> PDBModel
        mask - N.array( 1 x N_atoms of 1 or 0 )
               1 .. keep this atom
        deepcopy - deepcopy atom dictionaries of this model and result
        -> PDBModel
        """
        return self.take( N.nonzero( mask ), deepcopy )


    def remove( self, what ):
        """
        Convenience access to the 3 different remove methods. 
        what - function( atom_dict ) -> 1 || 0    (1..remove) OR
             - list of int [4, 5, 6, 200, 201..], indices of atoms to remove
             - list of int [11111100001101011100..N_atoms], mask (1..remove)
             - int, remove atom with this index

        -> N.array(1 x N_atoms_old) of 0||1,  mask used to compress the atoms
           and xyz arrays. This mask can be used to apply the same change
           to another array of same dimension as the old(!) xyz and atoms.
        !! raise PDBError, if what is neither of above
        """
        mask = N.logical_not( self.mask( what ) )
        self.keep( N.nonzero(mask) )
        return mask
    
    
    def takeChains( self, chainLst, deepcopy=0, breaks=0, maxDist=None ):
        """
        Get copy of this model with only the given chains.
        chainLst - list of int
        deepcopy - 0||1, deepcopy atom dictionaries (0)
        breaks   - 0|1, split chains at chain breaks
        maxDist  - float, (if breaks=1) chain break threshold in A
        -> PDBModel (by default, all fields are shallow copies, see N.take() )
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
                    t.errWriteln("Warning addChainId(): Problem with atom "\
                                 +str(a)+"Error: " + t.lastError() )

        self.atomsChanged = 1

    def addChainId( self, first_id=None, keep_old=0, breaks=0 ):
        """
        Assign chain identifiers A - Z to all atoms.
        first_id  - str (A - Z), first letter instead of 'A'
        keep_old  - 1|0, don't override existing chain IDs [0]
        breaks    - 1|0, consider chain break as start of new chain [0]
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
        mask  - [ 0||1 x N_atoms ] atom mask to apply BEFORE
        start - int, starting number [1]
        addChainId - 1|0, add chain IDs if they are missing
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
        """-> int, number of atoms"""
        if not self.xyz is None:
            return len( self.xyz )

        if not self.atoms is None:
            return len( self.atoms )

        if self.source is None:
            return 0
        
        return len( self.getXyz() )

    def lenResidues( self ):
        """-> int, total number of residues"""
        try:
            return self.resMap()[-1] + 1
        except IndexError:  ## empty residue map
            return 0

    def lenChains( self, breaks=0, maxDist=None ):
        """
        breaks  - 0|1, detect chain breaks from backbone atom distances [0]
        maxDist - float, maximal distance between consequtive residues
                  [ None ] .. defaults to twice the average distance

        -> int, total number of chains
        """
        try:
            return self.chainMap( breaks=breaks, maxDist=maxDist )[-1] + 1
        except IndexError:  ## empty residue map
            return 0


    def resList( self, mask=None ):
        """
        Return list of lists of atom dictionaries per residue,
        which allows to iterate over residues and atoms of residues.
        mask - [ 0||1 x N_atoms ] atom mask to apply BEFORE
        
        -> [ [ {'name':'N', 'residue_name':'LEU', ..},
               {'name':'CA','residue_name':'LEU', ..} ],
             [ { 'name':'CA', 'residue_name':'GLY', .. }, .. ] ]
        """
        ri = N.concatenate( (self.resIndex( mask=mask ), [self.lenAtoms()] ) )
        resLen = len( ri ) - 1
        atoms = self.getAtoms()
        if mask:
            atoms = N.compress( mask, atoms ).tolist()

        return [ atoms[ ri[res] : ri[res+1] ] for res in range( resLen ) ] 


    def resModels( self ):
        """
        -> list of PDBModels, one for each residue
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
        mask - [00111101011100111...] consider atom: yes or no
               len(mask) == N_atoms
        -> list ala [000111111333344444..] with residue number for each atom 
        """
        result = []
        for a in self.getAtoms():
            result.append( a['residue_number'] )
        
        ## by default take all atoms
        mask = mask or N.ones( len(self.getAtoms() ) , 'i' )

        ## apply mask to this list
        return N.compress( mask, N.array(result, 'i') )


    def __calcResMap( self, mask=None ):
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

        See resList() for an example of how to use the residue map.

        mask - [0000011111001111...] include atom: yes or no
               len(atom_mask) == number of atoms in self.atoms/self.xyz
        force- 0|1, recalculate map even if cached one is available [0]
        cache- 0|1, cache new map [1]

        -> array [00011111122223333..], residue index for each unmasked atom
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
        force - 1|0, re-calculate even if cached result is available [0]
        cache - 1|0, cache the result [1]
        mask  - [0|1], atom mask to apply before (i.e. result indices refer to
                compressed model)
        -> [int], index of the last atom of each residue
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
        mask - [0|1], atom mask to apply before (i.e. result indices refer to
               compressed model)
        -> [int], index of the last atom of each residue
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
        breaks  - 1|0, split chains at chain breaks [0]
        maxDist - float, (if breaks=1) chain break threshold in A
        -> array 1 x N_atoms of int, e.g. [000000011111111111122222...]
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
        breaks  - 1|0, split chains at chain breaks [0]
        maxDist - float, (if breaks=1) chain break threshold in A
        -> array (1 x N_chains) of int
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
        breaks_only - 1|0, don't report ends of regular chains [1]
        maxDist     - float, maximal distance between consequtive residues
                      [ None ] .. defaults to twice the average distance
        -> [ int ], atom indices of last atom before a probable chain break
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
        resname - str or list of str, name of residue to be removed.
        """
        resname = t.toList( resname )
        self.remove( lambda a, res=resname: a['residue_name'] in res )
        

    def rms( self, other, mask=None, mask_fit=None, fit=1, n_it=1 ):
        """
        other    - PDBModel, other model to compare this one with
        mask     - [ int ], atom mask for rmsd calculation
        mask_fit - [ int ], atom mask for superposition (default: same as mask)
        fit      - 1||0, superimpose first (default: 1)
        n_it - int, number of fit iterations, kicking out outliers on the way
               [1] -> classic single fit, 0 -> until convergence
        -> float, rms in A
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
        refModel - PDBModel
        mask  - [ int ], atom mask for superposition
        n_it  - int, number of fit iterations, kicking out outliers on the way
                [1] -> classic single fit, 0 -> until convergence
        z - float, number of standard deviations for outlier definition [2]
        eps_rmsd - float, tolerance in rmsd [0.5]
        eps_stdv - float, tolerance in standard deviations [0.05]
        profname - str, name of new atom profile getting outlier flag
        -> array(3 x 3), array(3 x 1) - rotation and translation matrices
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
        rt - array( 4 x 4 ) OR array(3 x 3), array(3 x 1)
        -> PDBModel, with transformed coordinates
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
        Least-square fit this model onto refModel.
        refModel - PDBModel
        mask  - [ 1|0 ], atom mask for superposition
        n_it  - int, number of fit iterations, kicking out outliers on the way
                [1] -> classic single fit, 0 -> until convergence
        z - float, number of standard deviations for outlier definition [2]
        eps_rmsd - float, tolerance in rmsd [0.5]
        eps_stdv - float, tolerance in standard deviations [0.05]
        profname - str, name of new atom profile containing outlier flag
        -> PDBModel, with transformed coordinates
        """
        return self.transform(
            self.transformation( refModel, mask, n_it, eps_rmsd=eps_rmsd,
                                 eps_stdv=eps_stdv, profname=profname ) )
    

    def magicFit( self, refModel, mask=None ):
        """
        magicFit( refModel [, mask ] ) -> PDBModel (or subclass )
        Superimpose this model onto a ref. model with similar atom content.
        refModel - PDBModel
        mask     - [ 0||1 ] or array of 0||1, atom mask to use for the fit
        -> PDBModel or sub-class
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
        mask - [ 1|0 ], atom mask applied before calculating the center
        -> PDBModel, model with centered coordinates
        """
        r = self.clone()
        mask = mask or N.ones( len(self) )

        avg = N.average( N.compress( mask, r.getXyz(), 0 ) )

        r.setXyz( r.getXyz() - avg )

        return r


    def center( self, mask=None ):
        """
        mask - [ 1|0 ], atom mask applied before calculating the center
        -> (float, float, float), xyz coordinates of center
        """
        if mask == None:
            return N.average( self.getXyz() )

        return N.average( N.compress( mask, self.getXyz(), axis=0 ) )

    
    def centerOfMass( self ):
        """
        -> array('f')
        """
        M = self.masses()
        return mathUtils.wMean( self.getXyz(), M )


    def masses( self ):
        """
        -> array('f'), 1-D array with mass of every atom in 1/12 of C12 mass.
        !! PDBError, if the model contains elements of unknown mass
        """
        try:
            M = [ molUtils.atomMasses[a['element']]for a in self.getAtoms() ]

        except KeyError, why:
            raise PDBError('Cannot find mass for '+str(why))

        return N.array( M )


    def mass( self ):
        """
        -> float, total mass in 1/12 of C12 mass
        !! PDBError, if the model contains elements of unknown mass
        """
        return N.sum( self.masses )
    

    def residusMaximus( self, atomValues, mask=None ):
        """
        Take list of value per atom, return list where all atoms of any
        residue are set to the highest value of any atom in that residue.
        (after applying mask)
        atomValues - list 1 x N, values per atom
        mask - list 1 x N of 1||0, 'master' atoms of each residue
        -> array 1 x N of float
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
        cmpfunc  - function( self.atoms[i], self.atoms[j] ) -> -1, 0, +1
        
        -> [2,1,4,6,5,0,..] suggested position of each atom in re-sorted model
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
        sortArg  - comparison function
        deepcopy - 0||1, deepcopy atom dictionaries (0)
        -> copy of this model with re-sorted atoms (see N.take() )
        """
        sortArg = sortArg or self.argsort()

        terA = self.__terAtoms

        r = self.take( sortArg, deepcopy )
        
        r.__terAtoms = terA

        return r
    

    def unsort( self, sortList ):
        """
        Undo a previous sorting on the model itself (no copy).
        sortList - sort list used for previous sorting.

        -> the (back)sort list used ( to undo the undo...)
        !! raise PDBError, if sorting changed atom number
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
        -> ['C','CA','CB' .... ] 
        """
        ## By default return list of all atoms
        start = start or 0

        if stop == None:    ## don't use "stop = stop or ..", stop might be 0!
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
        condition_dic - {..}, key-value pairs to be matched
        dic           - {..}, dictionary to be tested
        -> 1|0, 1 if all key-value pairs of condition are matched in dic
        """
        for k,v in condition.items():
            if dic.get( k, None ) not in v:
                return 0
        return 1

    def __testDict_or( self, dic, condition ):
        """
        condition_dic - {..}, key-value pairs to be matched
        dic           - {..}, dictionary to be tested
        -> 1|0, 1 if any key-value pairs of condition are matched in dic
        """
        for k,v in condition.items():
            if dic.get( k, None ) in v:
                return 1
        return 0

    def filterIndex( self, mode=0, **kw ):
        """
        Get atom positions that match a combination of key=values.
        E.g. filter( chain_id='A', name=['CA','CB'] ) -> PDBModel
        mode - 0|1, 0 combine with AND, 1 combine with OR
        kw   - combination of atom dictionary keys and values/list of values
        -> [ int ]
        """
        ## cache to minimize function lookup
        atoms = self.getAtoms()
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
        mode - 0|1, 0 combine with AND, 1 combine with OR
        kw   - combination of atom dictionary keys and values/list of values
        -> PDBModel
        """
        return self.take( self.filterIndex( mode=mode, **kw ) )


    def equals(self, ref, start=None, stop=None):
        """
        Compares the residue and atom sequence in the given range. Coordinates
        are not checked, profiles are not checked.
        start, stop - res indices
        -> [1,0] - sequence identity 0|1, atom identity 0|1
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
        @todo option to make sure atoms are also in same order
        Apply to SORTED models without HETATOMS. Coordinates are not checked.
        Note: in some rare cases m1.equalAtoms( m2 ) gives a different result
              than m2.equalAtoms( m1 ). This is due to the used SequenceMatcher
              class.
        -> (mask, mask_ref), two atom masks for all equal (1) atoms in models
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
        E.g. m2 = m1.sort()    ## m2 has now different atom order
             i2, i1 = m2.compareAtoms( m1 )
             m1 = m1.take( i1 ); m2 = m2.take( i2 )
             m1.atomNames() == m2.atomNames()  ## m2 has again same atom order 
        -> ([int],[int]), indices, indices_ref 
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
        -> float, average relative length of matching chain fragments 
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
        ref      - PDBModel
        breaks   - 1|0, look for chain breaks in backbone coordinates [0]
        fractLimit - 
        -> ([int], [int]), chainIndices, chainIndices_ref
        """
        i, i_ref = self.compareAtoms( ref )

        c0  = self.atom2chainIndices(i, breaks=breaks)
        c_r = ref.atom2chainIndices(i_ref,breaks=breaks)

        ## dirty hack to throw out single matching residues
        c0 = [ c for c in c0 \
               if self.__chainFraction( c, ref ) > fractLimit  ]
        c_r= [ c for c in c_r \
               if ref.__chainFraction( c, self ) > fractLimit  ]

        return c0, c_r
        
###############
## Testing
if __name__ == '__main__':

    t0 = time.time()

    ## loading output file from X-plor
    print 'Loading pdb file ..'
    m = PDBModel( t.testRoot()+'/rec/1A2P.pdb')#"/com_wet/1BGS.pdb")

    ## remove solvent
    print "removing waters.."
    m.removeRes(['TIP3', 'HOH'])

    ## X-plor doesn't write chainIds, so during the simulation
    ## we store them in the last letter of the segId. Here we
    ## restore the chainId.
    m.addChainFromSegid()
    
    ## start positions of all chains
    chainIdx = m.chainIndex()

    ## print some chain info
    print 'The molecule consists of %i chains'% m.lenChains()
    print '\tChainId \tFirst atom'
    for i in chainIdx:
        print '\t%s \t\t%i'%(m.atoms[i]['chain_id'], i)

    print "\ndone in ", time.time() - t0, "s"

    ## add accessibility and curvature
    from PDBDope import *
    m = m.compress( N.logical_not(m.maskSolvent() ) )
    d = PDBDope( m )
    d.addASA()
    d.addSurfaceRacer()   

    ## iterate over all chains
    for c in range( 0, len( chainIdx ) ):

        print "chain ", c, " starts with ", 
        print m.atoms[ chainIdx[c] ]['residue_name'],

        print " and has sequence: "

        ## mask out atoms of all other chains
        chainMask  = N.equal( m.chainMap( breaks=1 ), c )
        print m.sequence( chainMask )

    ## test sorting
    print "sorting atoms alphabetically..."
    sort = m.argsort()
    m2 = m.sort( sort )

    print "done in ", time.time() - t0, "s"
