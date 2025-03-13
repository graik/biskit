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
Trajectory - Collection of coordinate frames of a molecule 
"""

import biskit
import biskit.core.oldnumeric as N0
from biskit.core.trajparserFactory import TrajParserFactory

## superposition module from M. Habeck
from biskit import rmsFit

from biskit import tools as T
from biskit import mathUtils as MU
from biskit.errors import BiskitError
from biskit import EHandler, PDBModel, PDBError, ProfileCollection

import string
import re
import copy
import tempfile, os, types

## PCA
import numpy as N
import numpy.linalg as LA

## read Amber CRD and NC (netcdf) trajectories
import parmed.amber

class TrajError( BiskitError ):
    pass

class TrajProfiles( ProfileCollection ):
    pass

class Trajectory:
    """
    Manage many conformations of one molecule.
    Read coordinates and PDB infos from many PDB files,
    superimpose all coordinates to reference structure,
    calculate pairwise RMSD of all conformations ..
    """

    ## regExpression needed for sorting frames by their names
    ## used by __cmpFileNames()
    ex_numbers = re.compile(r'\D*([0-9]+)\D*')

    def __init__( self, source=None, refpdb=None, hasbox=True,
                  rmwat=1, castAll=0,  verbose=True ):
        """
        Collect coordinates into Numpy array. By default, the atom content
        of the first PDB is compared to the reference PDB to look for atoms
        that have to be removed or re-ordered. It is assumed that the other
        PDB files in the list have an identical atom content/layout and the
        same re-ordering / removing is applied to all of them. Set castAll
        to 1, in order to check each PDB seperately.

        :param source: file name of Amber MD or list of PDBs (files or models)
        :type  pdbs: str OR [ str ] OR [ PDBModel ]
        :param refpdb: file name of reference pdb
        :type  refpdb: str
        :param hasbox: for ASCII-formatted input trajectory - assume periodic 
                       box info in the input file (default True)
        :type hasbox: boolean
        :param rmwat: skip all TIP3, HOH, Cl-, Na+ from all files (default: 1)
        :type  rmwat: 0|1
        :param castAll: re-analyze atom content of each frame (default: 0)
        :type  castAll: 0|1
        :param verbose: verbosity level (default: 1)
        :type  verbose: 1|0
        """
        self.ref = None
        self.frames = None
        self.resIndex = None
        self.frameNames = None
        self.pc = None
        self.profiles = TrajProfiles()
        self.verbose= verbose

        if source is not None:
            p = TrajParserFactory.getParser(source, 
                                            verbose=verbose,
                                            rmwat=rmwat,
                                            hasbox=hasbox,
                                            analyzeEach=castAll)
            p.parse2new(source, refpdb, self)

        ## version as of creation of this object
        self.initVersion = T.dateString() + ';' + biskit.__version__


    def __getitem__( self, i ):
        """
        Get a single frame and support slicing of tajectories. 

        :param i: index OR SliceTyp
        :type  i: int OR [int]

        :return: model OR trajectory
        :rtype: PDBModel OR Trajectory  

        :raise TrajError: if out of memory OR invalid index
        """
        try:
            if isinstance( i , (int, N0.integer)):
                return self.getPDBModel( i )

            if type( i ) is slice:
                start = i.start or 0
                stop  = i.stop  or self.lenFrames()
                if stop > self.lenFrames():
                    stop = self.lenFrames()
                step  = i.step  or 1
                return self.takeFrames( list(range( start, stop, step)) )

        except MemoryError:
            raise TrajError("out of memory, cannot extract frames %s" % str(i))

        raise TrajError("unsupported index or slice: %s" % str( i ))


    def __len__( self ):
        """
        Number of frames in the trajectory.

        :return: length
        :rtype: int
        """
        return self.lenFrames()


    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        self.__dict__ = state
        ## backwards compability
        self.__defaults() 


    def __defaults(self ):
        """
        backwards compatibility to earlier pickled trajectories
        """
        self.pc = getattr( self, 'pc', None )
        self.frameNames = getattr( self, 'frameNames', None)
        self.profiles = getattr( self, 'profiles', TrajProfiles() )

        if type( self.frames ) is not N0.ndarray:
            self.frames = N0.array( self.frames )
        if type( self.resIndex ) is not N0.ndarray:
            self.resIndex = N0.array( self.resIndex )


    def __getstate__(self):
        """
        Called before pickling the object.
        """
        try:
            if type( self.frames ) == list or self.frames.dtype.char == 'd':
                EHandler.warning("Converting coordinates to float array.")
                self.frames = N0.array( self.frames ).astype(N0.Float32)
        except:
            EHandler.warning('Could not convert frames to float array.', 1)

        return self.__dict__


    def avgModel( self ):
        """
        Return a PDBModel with coordinates that are the average of
        all frames.

        :return: PDBModel with average structure of trajectory (no fitting!) 
                 this trajectory's ref is the source of result model
        :rtype: PDBModel
        """
        result = PDBModel( self.getRef(), noxyz=1 )
        result.setXyz( N0.average( self.frames ) )

        return result


    def getRef( self ):
        """
        :return: reference PDBModel
        :rtype: PDBModel
        """
        return self.ref


    def setRef( self, refModel ):
        """
        Assign new reference model.

        :param refModel: PDBModel with same number of atoms as in the frames.
        :type  refModel: PDBModel

        :return: old reference PDBModel
        :rtype: PDBModel
        """
        if (self.ref and refModel.equals( self.ref ) != [1,1] ) or \
           (self.frames is not None and refModel.lenAtoms() != len( self.frames[0] ) ):

            raise TrajError(\
                "Atoms of reference model don't match trajectory frames.")

        old, self.ref = self.ref, refModel

        return old


    def concat( self, *traj ):
        """
        Concatenate this with other trajectories. The ref model of the
        new Trajectory is a 'semi-deep' copy of this trajectorie's model.
        (see :class:`PDBModel.take()` )::
           concat( traj [, traj2, traj3, ..] ) -> Trajectory 

        :param traj: one or more Trajectory with identical atoms as this one
        :type  traj: Trajectories

        :return: concatenated trajecties
        :rtype: Trajectory
        """
        if len( traj ) == 0:
            return self

        r = self.__class__()

        r.frames = N0.concatenate( (self.frames, traj[0].frames), 0 )

        r.setRef( self.ref.clone())

        if self.frameNames and traj[0].frameNames:
            r.frameNames = self.frameNames + traj[0].frameNames

        try:
            if self.pc is not None and traj[0].pc is not None:
                r.pc['p'] = N0.concatenate( (self.pc['p'], traj[0].pc['p']),0)
                r.pc['u'] = N0.concatenate( (self.pc['u'], traj[0].pc['u']),0)
        except TypeError as why:
            EHandler.error('cannot concat PC '+str(why) )

        r.profiles = self.profiles.concat( traj[0].profiles )

        ## recursively add other trajectories
        return r.concat( *traj[1:] )


    def concatAtoms( self, *traj ):
        """
        Concatenate 2 trajectories of same (frame) length 'horizontally', i.e.
        for each frame the coordinates of one are appended to the coordinates
        of the other. The ref model of the new trajectory is a 'semi-deep' copy
        of this trajectory's model (see :class:`PDBModel.take()` )::
          concatAtoms( traj1 [traj2, traj3..]) -> Trajectory

        :param traj: one or more Trajectory of the same number of frames
        :type  traj: Trajectories

        :return: trajectory with concatenated atoms
        :rtype: Trajectory        
        """
        if len( traj ) == 0:
            return self

        r = self.__class__()

        r.frames = N0.concatenate( (self.frames, traj[0].frames), 1 )
        r.setRef( self.ref.concat( traj[0].getRef() ) )

        r.profiles = self.profiles.clone()

        r.frameNames = self.frameNames

        return r.concatAtoms( *traj[1:] )


    def lenFrames( self ):
        """
        :return: number of frames in trajectory
        :rtype: int
        """
        return len( self.frames )


    def lenAtoms( self ):
        """
        :return: number of atoms in frames
        :rtype: int
        """
        return N0.shape( self.frames )[1]


    def atomMask( self, what ):
        """
        Get atom mask.

        :param what: Create the mask using::
                      - funct( self.ref.atoms[i] ) -> int
                      - list int ( indices )
                      - mask
                      - list of str (atom names)
        :type  what: any

        :return: list of 1|0  1 x N_atoms
        :rtype: list[1|0]
        """
        return self.getRef().mask( what )


    def resMap( self, force=0 ):
        """
        A list with residue numbers mapped to each position. i.e.
        C{ [00001111222222222333..] }

        :param force: calculate even if it already been calculated
        :type  force: 1|0

        :return: list of int
        :rtype:  [int]
        """
        if not self.resIndex or not force:
            self.resIndex = self.getRef().resMap()

        return self.resIndex


    def takeFrames( self, indices ):
        """
        Return a copy of the trajectory containing only the specified frames.

        :param indices: positions to take
        :type  indices: [int]

        :return: copy of this Trajectory (fewer frames, semi-deep copy of ref)
        :rtype: Trajectory
        """
        ## remove out-of-bound indices
        indices = N0.compress( N0.less( indices, len( self.frames) ), indices )

        r = self.__class__()

        ## this step takes some time for large frames !
        r.frames = N0.take( self.frames, indices, 0 )

        ## semi-deep copy of reference model
        r.setRef( self.ref.take( list(range( self.ref.lenAtoms()))) )

        if self.frameNames is not None:
            r.frameNames = N0.take( self.frameNames, indices, 0 )
            r.frameNames = list(map( ''.join, r.frameNames.tolist() ))

        r.pc = self.__takePca( indices )

        r.profiles = self.profiles.take( indices )

        r.resIndex = self.resIndex

        return r


    def clone( self ):
        """
        Copy trajectory.

        :return: Trajectory (or sub-class), copy of this trajectory
        :rtype: Trajectory
        """
        return self.takeFrames( list(range( self.lenFrames())) )


    def compressFrames( self, mask ):
        """
        Compress trajectory with a frame mask. 

        :param mask: frame mask, 1 x N_frames
        :type  mask: [1|0]

        :return: copy of this Trajectory (fewer frames, semi-deep copy of ref)
        :rtype: Trajectory
        """
        return self.takeFrames( N0.nonzero( mask ) )


    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.

        :param traj: trajectory
        :type  traj: trajectory
        """
        self.frames = traj.frames
        self.ref = traj.ref
        self.frameNames = traj.frameNames
        self.pc = traj.pc
        self.profiles = traj.profiles


    def keepFrames( self, indices ):
        """
        in-place version of :class:`takeFrames`. keepFrames( indices ) -> None

        :param indices: frame numbers
        :type  indices: [int]
        """
        r = self.takeFrames( indices )
        self.replaceContent( r )


    def removeFrames( self, indices ):
        """
        Remove given frames from this trajectory object.

        :param indices: frame numbers
        :type  indices: [int]
        """
        i = N0.arange( self.lenFrames())
        i = N0.delete( i, indices )
        self.keepFrames( i )

    def takeAtoms( self, indices, returnClass=None ):
        """
        Take atoms from frames::
          takeAtoms( indices, type=None ) -> copy of Trajectory

        :param indices: list of atom indices
        :type  indices: [int]
        :param returnClass: default: None, same class as this object
        :type  returnClass: class OR None

        :return: copy of this Trajectory (with fewer atoms)
        :rtype: Trajectory        
        """

        returnClass = returnClass or self.__class__
        r = returnClass()

        ## copy over everything, so that child classes can preserve own fields
        r.__dict__.update( self.__dict__ )
        r.frames = r.ref = r.frameNames = r.profiles = None

        r.frames = N0.take( self.frames, indices, 1 )

        r.setRef( self.ref.take( indices ) )

        r.frameNames = copy.copy( self.frameNames )
        r.resIndex = None

        r.profiles = self.profiles.clone()

        r.pc = self.pc   ## this is not really clean

        return r


    def compressAtoms( self, aMask, returnClass=None ):
        """
        Get copy of this trajectory with only atoms marked 1 in aMask.

        :param aMask: atom mask [10011100101111...],
                      lst 1 x N_atoms of 1(keep) or 0
        :type  aMask: [1|0]
        :param returnClass: default: None, same class as this object
        :type  returnClass: class

        :return: copy of Trajectory with fewer atoms
        :rtype: Trajectory
        """
        return self.takeAtoms( N0.nonzero( aMask ), returnClass )


    def keepAtoms( self, indices ):
        """
        in-place version of :class:`takeAtoms`. keepAtoms( indices ) -> None

        :param indices: atom numbers
        :type  indices: [int]
        """
        r = self.takeAtoms( indices )
        self.replaceContent( r )


    def removeAtoms( self, what ):
        """
        Remove atoms from all frames of trajectory and from reference
        structure.

        :param what: Specify what atoms to remove::
                      - function( atom_dict ) -> 1 || 0    or (1..remove)
                      - list of int [4, 5, 6, 200, 201..], indices of atoms
                          to remove
                      - list of int [11111100001101011100..N_atoms], mask
                          (1..remove)
                      - int, remove atom with this index
        :type  what: any


        :return: N0.array(1 x N_atoms_old) of 0||1, mask used to compress the
                 atoms and xyz arrays. This mask can be used to apply the
                 same change to another array of same dimension as the
                 old(!) xyz and atoms.
        :rtype: array
        """
        ## pass what on to PDBModel, collect resulting mask
        mask = N0.logical_not( self.atomMask( what ) )

        self.keepAtoms( N0.nonzero( mask ) )

        return mask


    def takeChains(self, chainLst, returnClass=None ):
        """
        Extract some chains from complete Trajectory.

        :param chainLst: chains to isolate
        :type  chainLst: [int]
        :param returnClass: default: None, same class as this object
        :type  returnClass: class

        :return: Trajectory  with only those chains
        :rtype: Trajectory
        """
        chainMap = self.ref.chainMap()

        atomMask = list(map( lambda i, allowed=chainLst: i in allowed, chainMap))

        return self.compressAtoms( atomMask, returnClass )


    def fit( self, mask=None, ref=None, n_it=1,
             prof='rms', verbose=1, fit=1, **profInfos ):
        """
        Superimpose all coordinate frames on reference coordinates. Put rms
        values in a profile. If n_it > 1, the fraction of atoms considered
        for the fit is put into a profile called |prof|_considered
        (i.e. by default 'rms_considered').

        :param mask: atom mask, atoms to consider default: [all]
        :type  mask: [1|0]
        :param ref: use as reference, default: None, average Structure
        :type  ref: PDBModel
        :param n_it: number of fit iterations, kicking out outliers on the way
                     1 -> classic single fit, 0 -> until convergence
                     (default: 1)
        :type  n_it: int
        :param prof: save rms per frame in profile of this name, ['rms']
        :type  prof: str
        :param verbose: print progress info to STDERR (default: 1)
        :type  verbose: 1|0
        :param fit: transform frames after match, otherwise just calc rms
                    (default: 1)          
        :type  fit: 1|0
        :param profInfos: additional key=value pairs for rms profile info []
        :type profInfos: key=value
        """
        if ref is None:
            refxyz = N0.average( self.frames, 0 )
        else:
            refxyz = ref.getXyz()

        if mask is None:
            mask = N0.ones( len( refxyz ), N0.Int32 )

        refxyz = N0.compress( mask, refxyz, 0 )

        if verbose: T.errWrite( "rmsd fitting..." )

        rms = []          ## rms value of each frame
        non_outliers = [] ## fraction of atoms considered for rms and fit
        iterations = []   ## number of iterations performed on each frame

        for i in range(0, len( self.frames) ):

            xyz = self.frames[i]

            if n_it != 1:
                (r, t), rmsdList = rmsFit.match( refxyz,
                                                 N0.compress( mask, xyz, 0), n_it)
                iterations.append( len( rmsdList ) )
                non_outliers.append( rmsdList[-1][0] )

                xyz_transformed = N0.dot( xyz, N0.transpose(r)) + t

                rms += [ rmsdList[-1][1] ]

            else:
                r, t = rmsFit.findTransformation( refxyz,
                                                  N0.compress( mask, xyz, 0))

                xyz_transformed = N0.dot( xyz, N0.transpose(r)) + t

                d = N0.sqrt(N0.sum(N0.power( N0.compress(mask, xyz_transformed,0)\
                                             - refxyz, 2), 1))


                rms += [ N0.sqrt( N0.average(d**2) ) ]


            if fit:
                self.frames[i] = xyz_transformed.astype(N0.Float32)

            if verbose and i%100 == 0:
                T.errWrite( '#' )

        self.setProfile( prof, rms, n_iterations=n_it, **profInfos )

        if non_outliers:
            self.setProfile( prof+'_considered', non_outliers,
                             n_iterations=n_it,
                             comment='fraction of atoms considered for iterative fit' )

        if verbose: T.errWrite( 'done\n' )


    def transform( self, *rt ):
        """
        Apply given transformation to all frames (in place).

        :param rt: rotation translation matrix
        :type  rt: array( 4 x 4 ) OR array(3 x 3), array(3 x 1)
        """
        if len(rt) == 2:
            r, t = rt[0], rt[1]
        else:
            rt = rt[0]
            r, t = (rt[0:3,0:3], rt[0:3, 3])

        r = N0.transpose( r )
        r = r.astype(N0.Float32)
        t = t.astype(N0.Float32)

        for i in range( len( self.frames ) ):
            self.frames[ i ] = N0.array( N0.dot( self.frames[i], r ) ) + t 

##         self.frames = N0.array( [ N0.dot( f, r )  for f in self.frames ] )
##         self.frames += t


    def blockFit2ref( self, refModel=None, mask=None, conv=1e-6 ):
        """
        Fit trajectory until convergence onto it's own average and then
        transform the average of all frames onto the reference. To be used
        with parallell trajectories.

        :param refModel: Reference model (default: None)
        :type  refModel: PDBModel       
        :param mask: atom mask to apply before fitting
                     (default: None, all atoms)
        :type  mask: [1|0]
        :param conv: convergence cutoff creterion (default: 1e-6)
        :type  conv: float
        """
        self.fit( ref=self.ref )

        m_avg = self.avgModel()

        ## fit on average until it's not getting better
        d = 1.
        dd= 1.
        while dd >= conv:

            self.fit( ref=m_avg, mask=mask )
            m_new_avg = self.avgModel()

            oldD, d    = d, m_avg.rms( m_new_avg, mask=mask )
            T.flushPrint( "rms difference: %f" % d )

            dd = oldD - d
            m_avg = m_new_avg

        ## transform trajectory en block onto reference
        if refModel:
            T.flushPrint('fitting trajectory en-block onto reference...')

            if refModel.atomNames() != self.ref.atomNames():
                ref_i, i = refModel.compareAtoms( m_avg )

                refModel = refModel.take( ref_i )
                m_avg    = m_avg.take( i )

            r, t = m_avg.transformation( refModel, mask )
            self.transform( r, t )


    def writePdb( self, index, fname):
        """
        Write (possibly transformed) coordinates back to pdb.

        :param index: frame index in trajectory
        :type  index: int
        :param fname: name of new file
        :type  fname: str 
        """
        try:
            self.getPDBModel( index ).writePdb( fname )
        except:
            EHandler.error('Error writing %s.' % fname)


    def writePdbs( self, fname, frames=None):
        """
        Write coordinates to an NMR-style MODEL/ENDMDL pdb file.

        :param fname: name of new file
        :type  fname: str
        :param frames: frame indices (default: None, all)
        :type  frames: [int]
        """

        if frames is None:
            frames = list(range( self.lenFrames())) 

        ## open new file
        out = open(fname, 'w')
        ## fetch name for temporary file
        tmpFile = tempfile.mktemp()  

        for n in frames:

            ## write temporary pdb and open it
            self.writePdb( n, tmpFile )

            pdb_temp = open( tmpFile )

            out.write("MODEL%6i\n" % n)
            for line in pdb_temp:
                if line[:3] != 'END':       # filter out END
                    out.write(line)
            out.write("ENDMDL\n")

            ## remove temporary file
            pdb_temp.close()
            os.remove( tmpFile )

        out.write("END")
        out.close()


    def writeCrd( self, fname, frames=None ):
        """
        Write frames to Amber crd file (w/o box info).

        :param fname: output file name
        :type  fname: str
        :param frames: frame indices (default: all)
        :type  frames: [int]
        """
        if frames is None:
            frames = list(range( self.lenFrames()))

        template = " %7.3f" * 10 + '\n'

        ## open new file
        out = open( T.absfile(fname), 'wt')
        __write = out.write  ## cache function address for speed

        __write('\n')

        n_lines = None

        for fi in frames:
            f = N0.ravel( self.frames[ fi ] )

            if n_lines is None:
                n_lines = n_lines or len( f ) // 10
                overhang= ( 0 != len( f ) % 10 )
                i_lines = range( n_lines )

            for i in i_lines:
                __write( template % tuple( f[10*i:10*i+10] ) )

            if overhang:
                for x in f[i*10+10:]:
                    __write(" %7.3f" % x )
                __write('\n')

        out.close()


    def getPDBModel( self, index ):
        """
        Get PDBModel object for a particular frame of the trajectory.

        :param index: frame index
        :type  index: int

        :return: model
        :rtype: PDBModel
        """
        s = PDBModel( self.ref, noxyz=1 )
        s.setXyz( self.frames[ index ] )
        return s


    def setProfile( self, name, prof, mask=None, default=None, asarray=1,
                    comment=None, **moreInfo ):
        """
        Add/override profile.

        :param name: profile name
        :type  name: str
        :param prof: list of values
        :type  prof: [any]

        :param mask: list 1 x N_items of 0|1, if there are less values
                     than items, provide mask for missing values,
                     N0.sum(mask)==N_items
        :type  mask: [0|1]
        :param default: value for items masked.
        :type  default: any
        :param asarray: store as list (0), as array (2) or store numbers as
                        array but everything else as list (1) (default: 1)
        :type  asarray: 0|1|2
        :param comment: goes into info[name]['comment']
        :type  comment: str
        :param moreInfo: additional key-value pairs for info[name]
        :type  moreInfo:

        :raise ProfileError: if length of prof != N_residues
        """
        self.profiles.set( name, prof, mask, default, asarray=asarray,
                           comment=comment, **moreInfo )


    def profile( self, name, default=None ):
        """
        Get the values of a profile::
          get( name ) -> list of values

        :param name: profile name
        :type  name: str
        :param default: default result if no profile is found
        :type  default: any

        :raise ProfileError: if no profile is found with |name|
        """
        return self.profiles.get( name, default )


    def profileInfo( self, name ):
        """
        GEt information associated with a profile::
          profileInfo( name ) -> dict with infos about profile
        Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1|0

        :param name: profile name
        :type  name: str

        :raise ProfileError: if no profile is found with |name|
        """
        return self.profiles.getInfo( name )


    def setProfileInfo( self, name, **args ):
        """
        Add/Override infos about a given profile::
          e.g. setInfo('relASA', comment='new', params={'bin':'whatif'})

        :param name: profile name
        :type  name: str

        :raise ProfileError: if no profile is found with |name|
        """
        self.profiles.setInfo( name, **args )


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """
        Create a mask from a profile with optional cutoff values.

        :param profName: profile name
        :type  profName: str
        :param cutoff_min: lower cutoff value (default: None)
        :type  cutoff_min: float
        :param cutoff_max: upper cutoff value (default: None)
        :type  cutoff_max: float

        :return: mask lenFrames x 1|0
        :rtype: [1|0]
        """
        return self.profiles.profile2mask( profName, cutoff_min, cutoff_max )


    def plotProfile( self, *name, **args ):
        """
        Create a biggles.FramedPlot objetc from a profile. Display the plot
        with biggles.FramedPlot.show()::
          plotProfile(name1, [name2, ..],[arg1=x, arg2=y])->biggles.FramedPlot

        :param name: profile name
        :type  name: str

        :return: plot object
        :rtype: biggles.FramedPlot
        """
        return self.profiles.plot( *name, **args )


    def pairwiseRmsd( self, aMask=None, noFit=0 ):
        """
        Calculate rmsd between each 2 coordinate frames.

        :param aMask: atom mask
        :type  aMask: [1|0]
        :return: frames x frames array of float
        :rtype: array
        """
        frames = self.frames

        if aMask is not None:
            frames = N0.compress( aMask, frames, 1 )

        result = N0.zeros( (len( frames ), len( frames )), N0.Float32 )

        for i in range(0, len( frames ) ):

            for j in range( i+1, len( frames ) ):
                if noFit:
                    d = N0.sqrt(N0.sum(N0.power(frames[i]-frames[j], 2), 1))
                    result[i,j] = result[j,i] = N0.sqrt( N0.average(d**2) )

                else:
                    rt, rmsdLst = rmsFit.match( frames[i], frames[j], 1 )
                    result[i,j] = result[j,i] = rmsdLst[0][1]

        return result


    def getFluct_global( self, mask=None ):
        """
        Get RMS of each atom from it's average position in trajectory.
        The frames should be superimposed (fit() ) to a reference.

        :param mask: N x 1 list/Numpy array of 0|1, (N=atoms),
                     atoms to be considered.
        :type  mask: [1|0]

        :return: Numpy array ( N_unmasked x 1 ) of float.
        :rtype: array
        """
        frames = self.frames
        if mask is not None:
            frames = N0.compress( mask, frames, 1 )

        ## mean position of each atom in all frames
        avg = N0.average( frames )

        return N0.average(N0.sqrt(N0.sum(N0.power(frames - avg, 2), 2) ))


    def __resWindow( self, res, n_neighbores, rchainMap=None,
                     left_allowed=None, right_allowed=None ):
        """
        Get indices of all atoms of a residue and some atoms of its
        neighboring residues (if they belong to the same chain).

        :param res: residue index
        :type  res: int
        :param n_neighbores: number of residues to include right and left
        :type  n_neighbores: int
        :param right_allowed: array 1 x N_atoms of 1|0, possible neighbore
                              atoms
        :type  right_allowed: array
        :param left_allowed: array 1 x N_atoms of 1|0, possible neighbore atoms
        :type  left_allowed: array 
        :param rchainMap: array 1 x N_residues of int, chain id of each res
        :type  rchainMap: array

        :return: atoms of res, atoms of neighbores
        :rtype: [ int ], [ int ]
        """
        ## some defaults.. time-consuming..
        if rchainMap is None:
            rchainMap = N0.take( self.chainMap(), self.resIndex() )

        if left_allowed  is None: left_allowed = N0.nonzero( self.ref.maskBB() )
        if right_allowed is None: right_allowed= N0.nonzero( self.ref.maskBB() )

        ## atom indices of center residue
        result = self.ref.res2atomIndices( [ res ] ).tolist()

        ## get indices of neighbore residues that still belong to same chain
        l = self.ref.lenResidues()
        chain = rchainMap[res]

        outer_left = range( res-n_neighbores, res )
        outer_right= range( res+1, res+n_neighbores+1 )

        outer_left = [ i for i in outer_left  if i > 0 and rchainMap[i]==chain]
        outer_right= [ i for i in outer_right if i < l and rchainMap[i]==chain]

        ## convert to atom indices, filter them against allowed neighbore atoms
        if outer_left:
            outer_left = self.ref.res2atomIndices( outer_left )
            outer_left = MU.intersection( left_allowed,  outer_left )

        if outer_right:
            outer_right= self.ref.res2atomIndices( outer_right)
            outer_right= MU.intersection( right_allowed, outer_right)

        return result, outer_left + outer_right


    def getFluct_local( self, mask=None, border_res=1,
                        left_atoms=['C'], right_atoms=['N'], verbose=1 ):
        """
        Get mean displacement of each atom from it's average position after
        fitting of each residue to the reference backbone coordinates of itself
        and selected atoms of neighboring residues to the right and left.

        :param mask: N_atoms x 1 array of 0||1, atoms for which fluctuation
                     should be calculated
        :type  mask: array
        :param border_res: number of neighboring residues to use for fitting
        :type  border_res: int
        :param left_atoms: atoms (names) to use from these neighbore residues
        :type  left_atoms: [str]
        :param right_atoms: atoms (names) to use from these neighbore residues
        :type  right_atoms: [str]

        :return: Numpy array ( N_unmasked x 1 ) of float
        :rtype: array
        """
        if mask is None:
            mask = N0.ones( len( self.frames[0] ), N0.Int32 )

        if verbose: T.errWrite( "rmsd fitting per residue..." )

        residues = N0.nonzero( self.ref.atom2resMask( mask ) )

        ## backbone atoms used for fit
        fit_atoms_right = N0.nonzero( self.ref.mask( right_atoms ) )
        fit_atoms_left  = N0.nonzero( self.ref.mask( left_atoms ) )
        ## chain index of each residue
        rchainMap = N0.take( self.ref.chainMap(), self.ref.resIndex() )

        result = []

        for res in residues:

            i_res, i_border = self.__resWindow(res, border_res, rchainMap,
                                               fit_atoms_left, fit_atoms_right)

            try:
                if not len( i_res ): raise PDBError('empty residue')

                t_res = self.takeAtoms( i_res + i_border )

                i_center = range( len( i_res ) )

                mask_BB = t_res.ref.maskBB() * t_res.ref.maskHeavy()

                ## fit with border atoms ..
                t_res.fit( ref=t_res.ref, mask=mask_BB, verbose=0 )
                ## .. but calculate only with center residue atoms
                frames = N0.take( t_res.frames, i_center, 1 )

                avg = N0.average( frames )

                rmsd = N0.average(N0.sqrt(N0.sum(N0.power(frames - avg, 2), 2) ))

                result.extend( rmsd )

                if verbose: T.errWrite('#')

            except ZeroDivisionError:
                result.extend( N0.zeros( len(i_res), N0.Float32 ) )
                T.errWrite('?' + str( res ))

        if verbose: T.errWriteln( "done" )

        return result


    def residusMaximus( self, atomValues, mask=None ):
        """
        Take list of value per atom, return list where all atoms of any
        residue are set to the highest value of any atom in that residue.
        (after applying mask)

        :param atomValues: list 1 x N, values per atom
        :type  atomValues: [ float ]
        :param mask: list 1 x N, 0|1, 'master' atoms of each residue
        :type  mask: [1|0]

        :return: Numpy array 1 x N of float
        :rtype: array
        """
        if mask is None:
            mask = N0.ones( len( self.frames[0] ), N0.Int32 )

        ## eliminate all values that do not belong to the selected atoms
        masked = atomValues * mask

        result = []

        ## set all atoms of each residue to uniform value
        for res in range( 0, self.resMap()[-1]+1 ):

            ## get atom entries for this residue
            resAtoms = N0.compress( N0.equal( self.resMap(), res ), masked )

            ## get maximum value
            masterValue = max( resAtoms )

            result += resAtoms * 0.0 + masterValue

        return N0.array( result )


    def getGammaFluct( self, fluctList=None ):
        """
        Set value of all atoms of each residue to fluctuation of
        its gamma atom ( CG, SG, OG ).

        :param fluctList: 1x N, precalculated list of values or
                          None (will be calculated new)
        :type  fluctList: [float]

        :return: Numpy array 1 x N of float
        :rtype: [float]
        """
        if fluctList is None:
            fluctList = self.getFluct_local()

        ## define mask for gamma atoms in all Amino acids
        select = lambda a: a['name'] in ['CG', 'CG1', 'CG2', 'OG',
                                         'OG1', 'OG2','SG']
        CG_mask = self.mask( select )

        return self.residusMaximus( fluctList, CG_mask )


    def getResFluct( self, atomFluctList=None ):
        """
        Convert list of atomic fluctuations to list of residue
        fluctuation.

        :param atomFluctList: array 1 x N_atoms of float
        :type  atomFluctList: [float]

        :return: array 1 x N_residues of float
        :rtype: [float]

        :raise TrajError: if result length <> N_residues: 
        """
        if atomFluctList is None:
            atomFluctList = self.getFluct_global()

        ## Give all atoms of each res. the same fluct. value
        ## (the highest fluctuation of any backbone atom)
        result = self.residusMaximus( atomFluctList, self.ref.maskBB() )

        ## take first atoms only
        result = N0.take( result, self.ref.resIndex() )
##        result = N0.compress( self.ref.maskCA(), atomFluctList)

        ## check dimension
        if len( result ) != self.ref.lenResidues():
            raise TrajError(
                "getResFluct(): Length of result list (%i) <>" % len(result)+
                " number of residues (%i)." % self.ref.lenResidues() )

        return result


    def __cmpLists( self, l1, l2 ):
        """
        Compare to lists by their first, then second, etc item.

        :param l1: list
        :type  l1: [float]
        :param l2: list
        :type  l2: [float]

        :return:  result of comparison (-1 == l1[i] < l2[i])
        :rtype: [-1|0|1]
        """
        result = cmp( l1[0], l2[0] )

        if result == 0 and len(l1) > 1 and len(l2) > 1:
            return self.__cmpLists( l1[1:], l2[1:] )

        return result


    def __cmpFileNames( self, f1, f2 ):
        """
        Compare 2 file names by the numbers contained in them using the
        regular expression :class:`ex_numbers` (or as strings, if no numbers are
        found)::
           f1, f2 - strings, e.g. 'frame_1_188.pdb'

        :param f1: file name
        :type  f1: str
        :param f2: file name
        :type  f2: str

        :return: result of comparison (-1 == f1 < f2)
        :rtype: -1|0|+1
        """
        try:
            ## extract list of numbers from file names
            l1 = self.ex_numbers.findall( f1 )
            l2 = self.ex_numbers.findall( f2 )

            l1 = list(map( float, l1 ))
            l2 = list(map( float, l2 ))

            if len( l1 ) > 0 and len( l2 ) > 0:
                return self.__cmpLists( l1, l2 )

        except:
            pass

        return cmp( l1, l2 )  ## fall-back solution, default comparison


    def argsortFrames( self ):
        """
        Prepare sorting list by file names. Assuming the file names contains
        some numbers seperated by non-number characters. Sorting is done after
        Nr1 then Nr2 then..

        :return: list of frame sort order
        :rtype: numpy.array(int)
        """
        return N.argsort( self.frameNames )


    def sortFrames( self, sortList=None ):
        """
        Apply to a trajectory object to sort the frames and names according
        to the sortList.

        :param sortList: list to sort after (default: None)
        :type  sortList: [int]

        :raise TrajError: if sortList doesn't fit number of frames or names
        """
        if sortList is None:
            sortList = self.argsortFrames()

        if len(sortList) != len(self.frames) or\
           len(sortList) != len(self.frameNames) :
            raise TrajError("Unequal number of frames and names when sorting")

        ## sort according to sortList
        self.keepFrames( sortList )


    def sortAtoms( self, f_cmp=None ):
        """
        Sorts atoms **WITHIN** residues in reference and all frames.

        :param f_cmp: atom comparison function
                      C{ f_cmp( atoms[i], atoms[j]) -> -1|0|+1 }
                      (default: alphabetic ordering by atom['name'] )
        :type  f_cmp: function

        :raise PDBModelError: if sorting has changed number of atoms in
                              reference
        """
        ## Get the atom sort list from the reference structure
        atomSortList = self.ref.argsort( f_cmp )

        self.keepAtoms( atomSortList )


    def __takePca( self, indices ):
        """
        extract PCA results for certain frames.

        :param indices: frame indecies
        :type  indices: [int]

        :return: list of pca values
        :rtype: [float]        
        """
        result = copy.deepcopy( getattr(self, 'pc', None ))

        if result is not None:

            result['p'] = N0.take( result['p'], indices, 0 )

            result['u'] = N0.take( result['u'], indices, 0 )

            if result['fMask'] is not None:
                result['fMask'] = N0.take( result['fMask'], indices, 0 )

        return result


    def __compressPca( self, fMask ):
        """
        Compress PCA results using a frame mask.

        :param fMask: frame mask
        :type  fMask: [1|0]

        :return: list of pca values
        :rtype: [float]        
        """
        return self.__takePca( N0.nonzero( fMask ) )


    def getPca( self, aMask=None, fMask=None, fit=1 ):
        """
        Get the results form a principal component analysis.

        :param aMask: 1 x N_atoms of 1|0, atom mask, default: last one used
        :type  aMask: [1|0]
        :param fMask: 1 x N_frames of 1|0, frame mask, default: all
        :type  fMask: [1|0]
        :param fit: fit to average structure before doing the PC analysis
                    (default: 1)
        :type  fit: 1|0

        :return: Dictionary with results from the PC analysis::
                   dic {'p': projection of each frame in PC space,
                        'e': list of eigen values,
                        'fit':.., 'aMask':.., 'fMask':.. parameters used}
        :rtype: dict
        """
        if aMask is None:
            aMask = N0.ones( self.getRef().lenAtoms(), N0.Int32 )

        pc = getattr(self, 'pc', None)

        ## return chached result if parameters haven't changed
        if pc is not None and pc['fMask'] == fMask and pc['fit'] == fit and \
           aMask == pc['aMask']:

            return pc

        evectors, proj, evalues = self.pca( aMask, fMask, fit )

        pc = {}
        pc['aMask'] = aMask
        pc['fMask'] = fMask
        pc['fit'] = fit
        pc['p'] = proj
        pc['e'] = evalues
        pc['u'] = evectors

        self.pc = pc

        return pc


    def pca( self, atomMask=None, frameMask=None, fit=1 ):
        """
        Calculate principal components of trajectory frames.

        :param atomMask: 1 x N_atoms, [111001110..] atoms to consider
                         (default: all)
        :type  atomMask: [1|0]
        :param frameMask: 1 x N_frames, [001111..] frames to consider
                          (default all )
        :type  frameMask: [1|0]

        :return: (N_frames x N_frames), (1 x N_frames),
                 projection of each frame in PC space, eigenvalue of each PC
        :rtype: array, array, array
        """
        if frameMask is None: frameMask = N0.ones( len( self.frames ), N0.Int32 )

        if atomMask is None: atomMask = N0.ones(self.getRef().lenAtoms(),
                                                N0.Int32)

        if fit:
            self.fit( atomMask )

        refxyz = N0.average( self.frames, 0 )

        data = N0.compress( frameMask, self.frames, 0 )

        data = data - refxyz

        data = N0.compress( atomMask, data, 1 )

        ## reduce to 2D array
        data = N0.array( map( N0.ravel, data ) )

        V, L, U = LA.svd( data )

        return U, V * L, N0.power(L, 2)


    def pcMovie( self, ev, steps, factor=1., ref=0, morph=1 ):
        """
        Morph between the two extreme values of a single principal
        component.

        :param ev: EigenVector to visualize
        :type  ev: int
        :param steps: number of intermediate frames
        :type  steps: int
        :param factor: exageration factor (default: 1 = No exageration)
        :type  factor: float
        :param ref: take other eigenvecors from this frame (default: 1)
        :type  ref: int
        :param morph: morph between min and max (1) or take real values (0)
                      (default: 1)
        :type  morph: 1|0

        :return: Trajectory with frames visualizing the morphing.
        :rtype: Trajectory
        """
        fit = 1
        if self.pc is not None:
            fit = self.pc['fit']
        pc = self.getPca( fit=fit )

        ## eigenvectors (rows)
        U = pc['u']

        ## raveled and centered frames
        x_avg = N0.average(self.frames, 0)
        X = N0.array( [N0.ravel(x) for x in self.frames - x_avg] )

        ## ev'th eigenvector of reference frame
        alpha_0 = N0.dot( X[ref], U[ev] )

        ## list of deviations of ev'th eigenvector of each frame from ref
        alpha_range = N0.dot(X, U[ev]) - alpha_0

        ## get some representative alphas...
        if morph:
            a_min = factor * min(alpha_range)
            a_max = factor * max(alpha_range)
            delta = (a_max - a_min) / steps
            alpha_range = [ a_min + i*(delta) for i in range(0, steps) ]
        else:
            alpha_range = N0.sort( alpha_range )
            delta = len(alpha_range) / (steps * 1.0)
            alpha_range = [ alpha_range[ int(round( i*delta )) ]
                            for i in range(0,steps) ]

        ## scale ev'th eigenvector of ref with different alphas 
        Y = N0.array( [ X[ref] + alpha * U[ev] for alpha in alpha_range] )

        ## back convert to N x 3 coordinates
        Y = N0.reshape(Y, (Y.shape[0], -1, 3))
        Y = x_avg + Y

        result = self.__class__()
        result.ref = self.ref

        result.frames = Y
        return result


    def thin(self, step=1):
        return self.takeFrames( N.arange(0,len(self.frames), step) )

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test Adaptive clustering"""


    def test_Trajectory(self):
        """Trajectory test"""
        ## Loading
        self.traj = T.load(T.testRoot() + '/lig_pcr_00/traj.dat')

        ## sort frames after frameNames
        self.traj.sortFrames()

        ## sort atoms 
        self.traj.sortAtoms()

        ## remove waters
        self.traj = self.traj.compressAtoms(
            N0.logical_not( self.traj.ref.maskH2O()) )

        ## get fluctuation on a residue level
        r1 = self.traj.getFluct_local( verbose=self.local )

        ## fit backbone of frames to reference structure
        self.traj.fit( ref=self.traj.ref,
                       mask=self.traj.ref.maskBB(), verbose=self.local )

        self.assertAlmostEqual( N0.sum( self.traj.profile('rms') ),
                                58.101235746353879, 2 )

    def test_TrajectoryFromNCDF(self):
        self.traj2 = Trajectory(refpdb=T.testRoot('amber/md_netcdf/0.pdb'),
                                source=T.testRoot('amber/md_netcdf/heat.ncdf'))
        self.assertEqual(len(self.traj2), 6)
        
    def test_TrajectoryFromPDBs(self):
        f = T.testRoot() + '/amber/md_pdbs/'
        allfiles = os.listdir( f )
        pdbs = []
        for fn in allfiles:
            try:
                if (fn[-4:].upper() == '.PDB'):
                    pdbs += [f + fn]
            except:
                pass

        ref = pdbs[0]
        traj1 = Trajectory( pdbs[:3], ref, rmwat=False, verbose=self.local)
        
        models = [PDBModel(f) for f in pdbs]
        traj2 = Trajectory( models[:3], ref, rmwat=False, verbose=self.local)
        
        self.assertTrue(N.all(traj1.frames == traj2.frames))
        

if __name__ == '__main__':

    BT.localTest()

