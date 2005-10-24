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
"""Trajectory - Collection of coordinate frames of a molecule"""

import Numeric as N

## superposition module from M. Habeck
import rmsFit

import tools as T
import mathUtils as MU
from Errors import BiskitError
from Biskit import EHandler
from PDBModel import PDBModel
from ProfileCollection import ProfileCollection

import string
import re
import copy
import tempfile, os, types

## PCA
import LinearAlgebra as LA

class TrajError( BiskitError ):
    pass


class TrajProfiles( ProfileCollection ):
    def version( self ):
        return 'Trajectory $Revision$'


class Trajectory:
    """
    Manage many conformations of one molecule.
    Read coordinates and PDB infos from many PDB files,
    superimpose all coordinates to reference structure,
    calculate pairwise RMSD of all conformations ..
    """

    ## regExpression needed for sorting frames by their names
    ex_numbers = re.compile('\D*([0-9]+)\D*')

    def __init__( self, pdbs=None, refpdb=None, rmwat=1, castAll=0 ):
        """
        Collect coordinates into Numpy array. By default, the atom content
        of the first PDB is compared to the reference PDB to look for atoms
        that have to be removed or re-ordered. It is assumed that the other
        PDB files in the list have an identical atom content/layout and the
        same re-ordering / removing is applied to all of them. Set castAll
        to 1, in order to check each PDB seperately.
        pdbs    - [ str ], file names of all conformations
                - OR: [ PDBModel ]
        refpdb  - str, file name of reference pdb
        rmwat   - 0||1, skip all TIP3, HOH, Cl-, Na+ from all files (default 1)
        castAll - 0||1, re-analyze atom content of each frame (default 0)
        """
        self.ref = None
        self.frames = None
        self.resIndex = None
        self.frameNames = None
        self.pc = None
        self.profiles = TrajProfiles()

        if pdbs != None:
            refpdb = refpdb or pdbs[0]
            
            self.__create( pdbs, refpdb, rmwat=rmwat, castAll=castAll )

        ## version as of creation of this object
        self.initVersion = T.dateString() + ';' + self.version()


    def version( self ):
        return 'Trajectory $Revision$'


    def __create( self, pdbs, refpdb, rmwat=0, castAll=0 ):

        ## get Structure object for reference
        self.setRef( PDBModel( refpdb ) )

        ## Remove waters from ref (which will also remove it from each frame)
        if rmwat:
            wat = ['TIP3', 'HOH', 'WAT', 'Na+', 'Cl-' ]
            self.ref.remove( lambda a, wat=wat: a['residue_name'] in wat )

        ## frames x (N x 3) Array with coordinates
        self.frames = self.__collectFrames( pdbs, castAll )
        self.frames.savespace()

        ## [00011111111222233..] (continuous) residue number for each atom
        self.resIndex = self.ref.resMap()

        ## keep list of loaded files
        self.frameNames = pdbs


    def __getitem__( self, i ):
        try:
            if type( i ) is int:
                return self.getPDBModel( i )

            if type( i ) is types.SliceType:
                start = i.start or 0
                stop  = i.stop  or self.lenFrames()
                if stop > self.lenFrames():
                    stop = self.lenFrames()
                step  = i.step  or 1
                return self.takeFrames( range( start, stop, step ) )

        except MemoryError:
            raise TrajError, "out of memory, cannot extract frames %s" % str(i)

        raise TrajError, "unsupported index or slice: %s" % str( i )

    def __len__( self ):
        return self.lenFrames()

    def __setstate__(self, state ):
        """called for unpickling the object."""
        self.__dict__ = state
        ## backwards compability
        self.__defaults() 


    def __defaults(self ):
        """backwards compatibility to earlier pickled trajectories"""
        self.pc = getattr( self, 'pc', None )
        self.frameNames = getattr( self, 'frameNames', None)
        self.profiles = getattr( self, 'profiles', TrajProfiles() )
        try:
            self.frames.savespace()
        except:
            pass
        
    def __getstate__(self):
        """called before pickling the object."""
        try:
            if type( self.frames ) == list or self.frames.typecode() == 'd':
                EHandler.warning("Converting coordinates to float array.")
                self.frames = N.array( self.frames ).astype('f')
        except:
            EHandler.warning('Could not convert frames to float array.', 1)
            
        return self.__dict__


    def avgModel( self ):
        """
        -> PDBModel with average structure of trajectory (no fitting!)
        this trajectory's ref is the source of result model
        """
        result = PDBModel( self.getRef(), noxyz=1 )
        result.setXyz( N.average( self.frames ) )

        return result


    def __collectFrames( self, pdbs, castAll=0 ):
        """
        read coordinates from list of pdb files.
        pdbs - [str, str, str...]
        castAll - 0||1, analyze atom content of each frame for casting
        -> frames x (N x 3) Numpy array (of float)
        """
        frameList = []
        i = 0
        atomCast = None
        
        T.errWrite('reading %i pdbs...' % len(pdbs) )

        refNames = self.ref.atomNames()  ## cache for atom checking

        for f in pdbs:

            ## Load
            m = PDBModel(f)

            ## compare atom order & content of first frame to reference pdb
            if castAll or i==0:
                atomCast, castRef = m.compareAtoms( self.ref )

                if castRef != range( len( self.ref ) ):
                    ## we can take away atoms from each frame but not from ref
                    raise TrajError("Reference PDB doesn't match %s."
                                    %m.fileName)

                if atomCast == range( len( m ) ):
                    atomCast = None   ## no casting necessary
                else:
                    T.errWrite(' casting ')

            ## assert that frame fits reference
            if atomCast:
                m = m.take( atomCast )

            ## additional check on each 100st frame
            if i%100 == 0 and m.atomNames() <> refNames:
                raise TrajError("%s doesn't match reference pdb."%m.fileName )

            frameList.append( m.xyz )

            i += 1
            if i%10 == 0:
                T.errWrite('#')

        T.errWrite( 'done\n' )

        ## convert to 3-D Numpy Array
        return N.array(frameList).astype('f')


    def getRef( self ):
        """-> reference PDBModel"""
        return self.ref


    def setRef( self, refModel ):
        """
        Assign new reference model.
        refModel - PDBModel with same number of atoms as in the frames.
        -> old reference PDBModel
        """
        if (self.ref and refModel.equals( self.ref ) != [1,1] ) or \
           (self.frames and refModel.lenAtoms() != len( self.frames[0] ) ):

            raise TrajError(\
                "Atoms of reference model don't match trajectory frames.")
        
        old, self.ref = self.ref, refModel

        return old


    def concat( self, *traj ):
        """
        concat( traj [, traj2, traj3, ..] ) -> Trajectory
        Concatenate this with other trajectories. The ref model of the
        new Trajectory is a 'semi-deep' copy of this trajectorie's model.
        (see PDBModel.take() )
        traj - one or more Trajectory, with identical atoms as this one
        """
        if len( traj ) == 0:
            return self

        r = self.__class__()
        
        r.frames = N.concatenate( (self.frames, traj[0].frames), 0 )

        r.setRef( self.ref.clone())

        if self.frameNames and traj[0].frameNames:
            r.frameNames = self.frameNames + traj[0].frameNames

        try:
            if self.pc and traj[0].pc:
                r.pc['p'] = N.concatenate( (self.pc['p'], traj[0].pc['p']),0)
                r.pc['u'] = N.concatenate( (self.pc['u'], traj[0].pc['u']),0)
        except TypeError, why:
            EHandler.error('cannot concat PC '+str(why) )

        r.profiles = self.profiles.concat( traj[0].profiles )
    
        ## recursively add other trajectories
        return r.concat( *traj[1:] )


    def concatAtoms( self, *traj ):
        """
        concatAtoms( traj1 [traj2, traj3..]) -> Trajectory
        Concatenate 2 trajectories of same (frame) length 'horizontally', i.e.
        for each frame the coordinates of one are appended to the coordinates
        of the other. The ref model of the new trajectory is a 'semi-deep' copy
        of this trajectory's model (see PDBModel.take() ).
        """
        if len( traj ) == 0:
            return self

        r = self.__class__()
        
        r.frames = N.concatenate( (self.frames, traj[0].frames), 1 )
        r.setRef( self.ref.concat( traj[0].getRef() ) )

        r.profiles = self.profiles.clone()

        r.frameNames = self.frameNames

        return r.concatAtoms( *traj[1:] )
    

    def lenFrames( self ):
        return len( self.frames )

    def lenAtoms( self ):
        return N.shape( self.frames )[1]

    def atomMask( self, what ):
        """
        Get atom mask.
        what - funct( self.ref.atoms[i] ) -> int
             - list int ( indices )
             - mask
             - list of str (atom names)
        -> list of 1|0  1 x N_atoms
        """
        return self.getRef().mask( what )

    def resMap( self, force=0 ):
        """
        -> list of int
           [00001111222222222333..] map between atom and residue number
        """
        if not self.resIndex or not force:
            self.resIndex = self.getRef().resMap()

        return self.resIndex
    

    def takeFrames( self, indices ):
        """
        indices - list of int, positions to take
        -> copy of this Trajectory (fewer frames, semi-deep copy of ref)
        """
        ## remove out-of-bound indices
        indices = N.compress( N.less( indices, len( self.frames) ), indices )

        r = self.__class__()

        ## this step takes some time for large frames !
        r.frames = N.take( self.frames, indices, 0 )

        ## semi-deep copy of reference model
        r.setRef( self.ref.take( range( self.ref.lenAtoms() )) )

        if self.frameNames != None:
            r.frameNames = N.take( self.frameNames, indices, 0 )
            r.frameNames = map( ''.join, r.frameNames.tolist() )

        r.pc = self.__takePca( indices )

        r.profiles = self.profiles.take( indices )

        r.resIndex = self.resIndex
        
        return r


    def clone( self ):
        """
        -> Trajectory (or sub-class), copy of this trajectory
        """
        return self.takeFrames( range( self.lenFrames() ) )


    def compressFrames( self, mask ):
        """
        mask - 1 x N_frames
        -> copy of this Trajectory (fewer frames, semi-deep copy of ref)
        """
        return self.takeFrames( N.nonzero( mask ) )


    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.
        """
        self.frames = traj.frames
        self.ref = traj.ref
        self.frameNames = traj.frameNames
        self.pc = traj.pc
        self.profiles = traj.profiles
        

    def keepFrames( self, indices ):
        """in-place version of takeFrames. keepFrames( indices ) -> None"""
        r = self.takeFrames( indices )
        self.replaceContent( r )


    def removeFrames( self, indices ):
        """
        Remove given frames from this trajectory object.
        indices - lst of int, frame numbers
        """
        i = range( self.lenFrames() )
        i.remove( N.array(indices) )
        self.keepFrames( i )
        

    def takeAtoms( self, indices, returnClass=None ):
        """
        takeAtoms( indices, type=None ) -> copy of Trajectory
        indices     - [int], list of atom indices
        returnClass - class, (default: same class as this object)
        """

        returnClass = returnClass or self.__class__
        r = returnClass()

        ## copy over everything, so that child classes can preserve own fields
        r.__dict__.update( self.__dict__ )
        r.frames = r.ref = r.frameNames = r.profiles = None

        r.frames = N.take( self.frames, indices, 1 )

        r.setRef( self.ref.take( indices ) )

        r.frameNames = copy.copy( self.frameNames )
        r.resIndex = None

        r.profiles = self.profiles.clone()
        
        r.pc = self.pc   ## this is not really clean

        return r


    def compressAtoms( self, aMask, returnClass=None ):
        """
        Get copy of this trajectory with only atoms marked 1 in aMask.
        aMask - [10011100101111...], lst 1 x N_atoms of 1(keep) or 0
        returnType - class, (default: same class as this object)
        -> copy of Trajectory with fewer atoms
        """
        return self.takeAtoms( N.nonzero( aMask ), returnClass )


    def keepAtoms( self, indices ):
        """in-place version of takeAtoms. keepAtoms( indices ) -> None"""
        r = self.takeAtoms( indices )
        self.replaceContent( r )

        
    def removeAtoms( self, what ):
        """
        Remove atoms from all frames of trajectory and from reference
        structure.
        what - function( atom_dict ) -> 1 || 0    or (1..remove)
             - list of int [4, 5, 6, 200, 201..], indices of atoms to remove
             - list of int [11111100001101011100..N_atoms], mask (1..remove)
             - int, remove atom with this index

        -> N.array(1 x N_atoms_old) of 0||1,  mask used to compress the atoms
           and xyz arrays. This mask can be used to apply the same change
           to another array of same dimension as the old(!) xyz and atoms.
        """
        ## pass what on to PDBModel, collect resulting mask
        mask = N.logical_not( self.atomMask( what ) )

        self.keepAtoms( N.nonzero( mask ) )

        return mask


    def takeChains(self, chainLst, returnClass=None ):
        """
        Extract some chains from complete Trajectory.
        chainLst    - [int], chains to isolate
        returnClass - class, (default: same class as this object)
        -> Trajectory  with only those chains
        """
        chainMap = self.ref.chainMap()

        atomMask = map( lambda i, allowed=chainLst: i in allowed, chainMap)

        return self.compressAtoms( atomMask, returnClass )


    def fit( self, mask=None, ref=None, n_it=1,
             prof='rms', verbose=1, fit=1, **profInfos ):
        """
        Superimpose all coordinate frames on reference coordinates. Put rms
        values in a profile. If n_it > 1, the fraction of atoms considered
        for the fit is put into a profile called |prof|_considered
        (i.e. by default 'rms_considered').
        mask   - atoms to consider, default: [all]
        ref    - PDBModel, as reference, default: [average Structure]
        n_it   - int, number of fit iterations, kicking out outlyers on the way
                 [1] -> classic single fit, 0 -> until convergence
        prof   - str, save rms per frame in profile of this name, ['rms']
        verbose- 1|0, print progress info to STDERR [1]
        fit    - 1|0, transform frames after match, otherwise just calc rms [1]
        additional key=value pairs for rms profile info []
        """
        if ref == None:
            refxyz = N.average( self.frames, 0 )
        else:
            refxyz = ref.getXyz()
        
        if mask == None:
            mask = N.ones( len( refxyz ) )

        refxyz = N.compress( mask, refxyz, 0 )

        if verbose: T.errWrite( "rmsd fitting..." )
        
        superImposed = [] ## superimposed frames
        rms = []          ## rms value of each frame
        non_outliers = [] ## fraction of atoms considered for rms and fit
        iterations = []   ## number of iterations performed on each frame
        
        for i in range(0, len( self.frames) ):

            xyz = self.frames[i]

            if n_it != 1:
                (r, t), rmsdList = rmsFit.match( refxyz,
                                    N.compress( mask, xyz, 0), n_it)
                iterations.append( len( rmsdList ) )
                non_outliers.append( rmsdList[-1][0] )

                xyz_transformed = N.dot( xyz, N.transpose(r)) + t

                rms += [ rmsdList[-1][1] ]

            else:
                r, t = rmsFit.findTransformation( refxyz,
                                    N.compress( mask, xyz, 0) )
                
                xyz_transformed = N.dot( xyz, N.transpose(r)) + t

                d = N.sqrt(N.sum(N.power( N.compress(mask, xyz_transformed,0)\
                                         - refxyz, 2), 1))
                rms += [ N.sqrt( N.average(d**2) ) ]

            if fit:
                self.frames[i] = xyz_transformed.astype('f')

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
        rt - array( 4 x 4 ) OR array(3 x 3), array(3 x 1)
        """
        if len(rt) == 2:
            r, t = rt[0], rt[1]
        else:
            rt = rt[0]
            r, t = (rt[0:3,0:3], rt[0:3, 3])

        r = N.transpose( r )
        r = r.astype('f')
        t = t.astype('f')
                
        for i in range( len( self.frames ) ):
            self.frames[ i ] = N.array( N.dot( self.frames[i], r ) ) + t 

##         self.frames = N.array( [ N.dot( f, r )  for f in self.frames ] )
##         self.frames += t


    def blockFit2ref( self, refModel=None, mask=None, conv=1e-6 ):
        """
        Fit trajectory until convergence onto it's own average and then
        transform the average of all frames onto the reference.
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
        index - int, frame index in trajectory
        fname - name of new file
        """
        try:
            self.getPDBModel( index ).writePdb( fname )
        except:
            EHandler.error('Error writing %s.' % fname)

            
    def writePdbs( self, fname, frames=None):
        """
        Write coordinates to an NMR-style MODEL/ENDMDL pdb file.
        fname  - name of new file
        frames - list of int, frame indices (default: all)
        """

        if frames is None:
            frames = range( self.lenFrames() ) 

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
                if line[:3] <> 'END':       # filter out END
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
        fname  - str, output file name
        frames - list of int, frame indices (default: all)
        """
        if frames == None:
            frames = range( self.lenFrames() )

        template = " %7.3f" * 10 + '\n'

        ## open new file
        out = open( T.absfile(fname), 'w')
        __write = out.write  ## cache function address for speed

        __write('\n')

        n_lines = None

        for fi in frames:
            f = N.ravel( self.frames[ fi ] )

            if n_lines is None:
                n_lines = n_lines or len( f ) / 10
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
        index - int, frame index
        -> PDBModel
        """
        s = PDBModel( self.ref, noxyz=1 )
        s.setXyz( self.frames[ index ] )
        return s


    def setProfile( self, name, prof, mask=None, default=None, asarray=1,
                    comment=None, **moreInfo ):
        """
        Add/override profile.
        name - str
        prof - list of values
        mask     - list 1 x N_items of 0|1, if there are less values than
                   items, provide mask for missing values, N.sum(mask)==N_items
        default  - any, value for items masked.
        asarray  - 0|1|2, store as list (0), as array (2) or store
                   numbers as array but everything else as list (1) [1]
        comment  - str, goes into info[name]['comment']
        moreInfo - additional key-value pairs for info[name]
        !! raise ProfileError, if length of prof != N_residues
        """
        self.profiles.set( name, prof, mask, default, asarray=asarray,
                           comment=comment, **moreInfo )

    def profile( self, name, default=None ):
        """
        get( name ) -> list of values
        default - default result if no profile is found,
                  if None ..
        !! raise ProfileError, if no profile is found with |name|
        """
        return self.profiles.get( name, default )

    def profileInfo( self, name ):
        """
        profileInfo( name ) -> dict with infos about profile
        Guaranteed infos: 'version'->str, 'comment'->str, 'changed'->1||0
        !! raise ProfileError, if no profile is found with |name|
        """
        return self.profiles.getInfo( name )

    def setProfileInfo( self, name, **args ):
        """Add/Override infos about a given profile.
        e.g. setInfo('relASA', comment='new', params={'bin':'whatif'})
        !! raise ProfileError, if no profile is found with |name|
        """
        self.profiles.setInfo( name, **args )


    def profile2mask(self, profName, cutoff_min=None, cutoff_max=None ):
        """-> mask lenFrames x 1|0"""
        return self.profiles.profile2mask( profName, cutoff_min, cutoff_max )


    def plotProfile( self, *name, **args ):
        """plotProfile(name1, [name2, ..],[arg1=x, arg2=y])->biggles.FramedPlot
        """
        return self.profiles.plot( *name, **args )
    

    def pairwiseRmsd( self, aMask=None ):
        """
        Calculate rmsd between each 2 coordinate frames.
        aMask - atom mask
        -> frames x frames array of float
        """
        frames = self.frames

        if aMask != None:
            frames = N.compress( aMask, frames, 1 )
        
        result = N.zeros( (len( frames ), len( frames )), 'f' )
        
        for i in range(0, len( frames ) ):

            for j in range( i+1, len( frames ) ):

                rt, rmsdLst = rmsFit.match( frames[i], frames[j], 1 )
                result[i,j] = result[j,i] = rmsdLst[0][1]

        return result


    def getFluct_global( self, mask=None ):
        """
        Get RMS of each atom from it's average position in trajectory.
        The frames should be superimposed (fit() ) to a reference.
        mask - N x 1 list/Numpy array of 0||1, (N=atoms),
               atoms to be considered.
        -> Numpy array ( N_unmasked x 1 ) of float.
        """
        frames = self.frames
        if mask != None:
            frames = N.compress( mask, frames, 1 )

        ## mean position of each atom in all frames
        avg = N.average( frames )

        return N.average(N.sqrt(N.sum(N.power(frames - avg, 2), 2) ))


    def __resWindow( self, res, n_neighbores, rchainMap=None,
                     left_allowed=None, right_allowed=None ):
        """
        Get indices of all atoms of a residue and some atoms of its
        neighboring residues (if they belong to the same chain).
        res           - int, residue index
        n_neighbores  - int, number of residues to include right and left
        right_indices - array 1 x N_atoms of 1|0, possible neighbore atoms
        left_indices  - array 1 x N_atoms of 1|0, possible neighbore atoms
        rchainMap     - array 1 x N_residues of int, chain id of each res
        -> [ int ], [ int ] - atoms of res, atoms of neighbores
        """
        ## some defaults.. time-consuming..
        if rchainMap is None:
            rchainMap = N.take( self.chainMap(), self.resIndex() )

        left_allowed = left_allowed  or N.nonzero( self.ref.maskBB() )
        right_allowed= right_allowed or N.nonzero( self.ref.maskBB() )
        
        ## atom indices of center residue
        result = self.ref.res2atomIndices( [ res ] )

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
                        left_atoms=['C'], right_atoms=['N'] ):
        """
        Get mean displacement of each atom from it's average position after
        fitting of each residue to the reference backbone coordinates of itself
        and selected atoms of neighboring residues to the right and left.
        mask    - N_atoms x 1 array of 0||1, atoms for which fluctuation should
                  be calculated
        boder_res   - int, number of neighboring residues to use for fitting
        left_atoms  - [str], atoms (names) to use from these neighbore residues
        right_atoms - [str], atoms (names) to use from these neighbore residues
        -> Numpy array ( N_unmasked x 1 ) of float
        """

        if mask == None:
            mask = N.ones( len( self.frames[0] ), 'i' )

        T.errWrite( "rmsd fitting per residue..." )

        residues = N.nonzero( self.ref.atom2resMask( mask ) )

        ## backbone atoms used for fit
        fit_atoms_right = N.nonzero( self.ref.mask( right_atoms ) )
        fit_atoms_left  = N.nonzero( self.ref.mask( left_atoms ) )
        ## chain index of each residue
        rchainMap = N.take( self.ref.chainMap(), self.ref.resIndex() )

        result = []

        for res in residues:

            i_res, i_border = self.__resWindow(res, border_res, rchainMap,
                                          fit_atoms_left, fit_atoms_right)

            try:
                if not len( i_res ): raise PDBError, 'empty residue'

                t_res = self.takeAtoms( i_res + i_border )

                i_center = range( len( i_res ) )

                mask_BB = t_res.ref.maskBB() * t_res.ref.maskHeavy()

                ## fit with border atoms ..
                t_res.fit( ref=t_res.ref, mask=mask_BB, verbose=0 )
                ## .. but calculate only with center residue atoms
                frames = N.take( t_res.frames, i_center, 1 )

                avg = N.average( frames )

                rmsd = N.average(N.sqrt(N.sum(N.power(frames - avg, 2), 2) ))

                result.extend( rmsd )

                T.errWrite('#')

            except ZeroDivisionError:
                result.extend( N.zeros( len(i_res), 'f' ) )
                T.errWrite('?' + str( res ))

        T.errWriteln( "done" )

        return result
       

    def residusMaximus( self, atomValues, mask=None ):
        """
        Take list of value per atom, return list where all atoms of any
        residue are set to the highest value of any atom in that residue.
        (after applying mask)
        atomValues - list 1 x N, values per atom
        mask - list 1 x N, 0||1, 'master' atoms of each residue
        -> Numpy array 1 x N of float
        """
        if mask == None:
            mask = N.ones( len( self.frames[0] ), 'i' )

        ## eliminate all values that do not belong to the selected atoms
        masked = atomValues * mask

        result = []
        
        ## set all atoms of each residue to uniform value
        for res in range( 0, self.resMap()[-1]+1 ):

            ## get atom entries for this residue
            resAtoms = N.compress( N.equal( self.resMap(), res ), masked )
            
            ## get maximum value
            masterValue = max( resAtoms )

            result += resAtoms * 0.0 + masterValue

        return N.array( result )


    def getGammaFluct( self, fluctList=None ):
        """
        Set value of all atoms of each residue to fluctuation of
        its gamma atom ( CG, SG, OG ).
        fluctList - 1x N, precalculated list of values or
                    None (will be calculated new)
        -> Numpy array 1 x N of float
        """
        if fluctList == None:
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
        atomFluctList - array 1 x N_atoms of float
        -> array 1 x N_residues of float

        !! TrajError if result length <> N_residues
        """
        if atomFluctList == None:
            atomFluctList = self.getFluct_global()

        ## Give all atoms of each res. the same fluct. value
        ## (the highest fluctuation of any backbone atom)
        result = self.residusMaximus( atomFluctList, self.ref.maskBB() )

        ## take first atoms only
        result = N.take( result, self.ref.resIndex() )
##        result = N.compress( self.ref.maskCA(), atomFluctList)
        
        ## check dimension
        if len( result ) <> self.ref.lenResidues():
            raise TrajError(
                "getResFluct(): Length of result list (%i) <>" % len(result)+
                " number of residues (%i)." % self.ref.lenResidues() )

        return result


    def __cmpLists( self, l1, l2 ):
        """
        compare to lists by their first, then second, etc item.
        -> -1, 0, 1
        """
        result = cmp( l1[0], l2[0] )

        if result == 0 and len(l1) > 1 and len(l2) > 1:
            return self.__cmpLists( l1[1:], l2[1:] )

        return result


    def __cmpFileNames( self, f1, f2 ):
        """
        Compare 2 file names by the numbers contained in them (or as strings,
        if no numbers are found).
        f1, f2 - strings, e.g. 'frame_1_188.pdb'
        -> -1|0|+1
        """

        try:
            ## extract list of numbers from file names
            l1 = self.ex_numbers.findall( f1 )
            l2 = self.ex_numbers.findall( f2 )
            
            l1 = map( float, l1 )
            l2 = map( float, l2 )

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
       -> [int] 
       """
       names = self.frameNames

       result = range(0, len(names) )

       ## sort result but use items of names for the comparison
       f_cmp = lambda i,j, ns=names: self.__cmpFileNames( ns[i], ns[j]) 

       result.sort( f_cmp )

       return result


    def sortFrames( self, sortList=None ):
        """
        Apply to a trajectory object to sort the frames and names according
        to the sortList.
        !! raise TrajError if sortList doesn't fit number of frames or names
        """
        if sortList == None:
            sortList = self.argsortFrames()

        if len(sortList) != len(self.frames) or\
           len(sortList) != len(self.frameNames) :
            raise TrajError("Unequal number of frames and names when sorting")

        ## sort according to sortList
        self.keepFrames( sortList )
        

    def sortAtoms( self, f_cmp=None ):
        """
        Sorts atoms WITHIN residues in reference and all frames.
        f_cmp - atom comparison function f_cmp( atoms[i], atoms[j]) -> -1|0|+1
                (default: alphabetic ordering by atom['name'] )
        !! raise PDBModelError if sorting has changed number of atoms in
           reference
        """
        ## Get the atom sort list from the reference structure
        atomSortList = self.ref.argsort( f_cmp )

        self.keepAtoms( atomSortList )
        

    def __takePca( self, indices ):
        """
        extract PCA results for certain frames.
        """
        result = copy.deepcopy( getattr(self, 'pc', None ))

        if result != None:
            
            result['p'] = N.take( result['p'], indices, 0 )

            result['u'] = N.take( result['u'], indices, 0 )

            if result['fMask'] != None:
                result['fMask'] = N.take( result['fMask'], indices, 0 )
            
        return result


    def __compressPca( self, fMask ):

        return self.__takePca( N.nonzero( fMask ) )


    def getPca( self, aMask=None, fMask=None, fit=1 ):
        """
        aMask - 1 x N_atoms of 1||0, atom mask, default: last one used
        fMask - 1 x N_frames of 1||0, frame mask, default: all
        fit - fit to average structure before doing the PC analysis
        -> dic {'p': projection of each frame in PC space,
                'e': list of eigen values,
                'fit':.., 'aMask':.., 'fMask':.. parameters used}
        """
        if aMask == None:
            aMask = N.ones( self.getRef().lenAtoms() )
        
        pc = getattr(self, 'pc', None)

        ## return chached result if parameters haven't changed
        if pc != None and pc['fMask'] == fMask and pc['fit'] == fit and \
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
        atomMask - 1 x N_atoms, [111001110..] atoms to consider
        frameMask - 1 x N_frames, [001111..] frames to consider
                 (default all )
                 
        -> (N_frames x N_frames), (1 x N_frames)
           projection of each frame in PC space, eigenvalue of each PC
        """
        frameMask = frameMask or N.ones( len( self.frames ) )

        atomMask = atomMask or N.ones( self.getRef().lenAtoms() )

        if fit:
            self.fit( atomMask )

        refxyz = N.average( self.frames, 0 )
        
        data = N.compress( frameMask, self.frames, 0 )

        data = data - refxyz

        data = N.compress( atomMask, data, 1 )
        
        ## reduce to 2D array
        data = N.array( map( N.ravel, data ) )

        V, L, U = LA.singular_value_decomposition( data )

        return U, V * L, N.power(L, 2)


    def pcMovie( self, ev, steps, factor=1., ref=0, morph=1 ):
        """
        Morphe between the two extreme values of a single principal
        component.
        ev - int, EigenVector to visualize
        steps - int, number of intermediate frames
        factor - float, exageration factor (1 = No exageration)
        ref - int, take other eigenvecors from this frame
        morph - morph between min and max (1) or take real values (0)
        -> Trajectory
        """
        fit = 1
        if self.pc:
            fit = self.pc['fit']
        pc = self.getPca( fit=fit )

        ## eigenvectors (rows)
        U = pc['u']

        ## raveled and centered frames
        x_avg = N.average(self.frames, 0)
        X = N.array( [N.ravel(x) for x in self.frames - x_avg] )

        ## ev'th eigenvector of reference frame
        alpha_0 = N.dot( X[ref], U[ev] )

        ## list of deviations of ev'th eigenvector of each frame from ref
        alpha_range = N.dot(X, U[ev]) - alpha_0

        ## get some representative alphas...
        if morph:
            a_min = factor * min(alpha_range)
            a_max = factor * max(alpha_range)
            delta = (a_max - a_min) / steps
            alpha_range = [ a_min + i*(delta) for i in range(0, steps) ]
        else:
            alpha_range = N.sort( alpha_range )
            delta = len(alpha_range) / (steps * 1.0)
            alpha_range = [ alpha_range[ int(round( i*delta )) ]
                            for i in range(0,steps) ]

        ## scale ev'th eigenvector of ref with different alphas 
        Y = N.array( [ X[ref] + alpha * U[ev] for alpha in alpha_range] )

        ## back convert to N x 3 coordinates
        Y = N.reshape(Y, (Y.shape[0], -1, 3))
        Y = x_avg + Y

        result = self.__class__()
        result.ref = self.ref

        result.frames = Y
        return result


###################
## TESTING
###################

if __name__ == '__main__':

    import os, time

##     f = T.testRoot() + '/lig_pc2_00/pdb/'
##     allfiles = os.listdir( f )
##     pdbs = []
##     for fn in allfiles:
##         try:
##             if (fn[-7:].upper() == '.PDB.GZ'):
##                 pdbs += [f + fn]
##         except:
##             pass
    
##     ref = pdbs[0]
##     traj = Trajectory( pdbs[:3], ref, rmwat=0 )
 
    print "Loading"
    traj = T.Load(T.testRoot() + '/lig_pc2_00/traj.dat')

    t0 = time.time()
    
    traj.sortFrames()
    
    ## sort atoms 
    traj.sortAtoms()

    traj = traj.compressAtoms( N.logical_not( traj.ref.maskH2O()) )
    r1 = traj.getFluct_local()

    t_res = traj.takeAtoms( traj.ref.res2atomIndices(
        traj.ref.atom2resIndices( [0] ) ) )

    t_res.fit( ref=t_res.ref, mask=t_res.ref.maskBB() )

    print "done in %f s." % (time.time() - t0)
    
## ## eigenvectors (rows)
    
## U = u

## ## frames (centered)

## x_avg = N.average(traj.frames, 0)
## X = N.array([N.ravel(x) for x in traj.frames - x_avg])

## ## reference structure j, eigenvector of interest U[i]

## j = 0
## i = 0

## alpha_0 = N.dot(X[j], U[i])
## alpha_range = N.dot(X, U[i]) - alpha_0
## alpha_range = N.arange(-200., 200., 40.)
## Y = N.array([X[j] + alpha * U[i] for alpha in alpha_range])
## Y = N.reshape(Y, (Y.shape[0], -1, 3))
## Y = x_avg + Y

##     from Biskit.EnsembleTraj import traj2ensemble

##     T.flushPrint("Loading...")
##     traj = T.Load('~/interfaces/a11/com_pcr_00/traj_0.dat')
##     traj = traj.thin( step=50 )
##     traj = traj2ensemble( traj )
##     T.flushPrint("done\n")
