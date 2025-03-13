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
Multi-copy trajectory
"""

import types

from biskit import EHandler
import biskit.core.oldnumeric as N0
import biskit.mathUtils as M
import biskit.tools as T

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.md; __package__ = "biskit.md"

from .trajectory import Trajectory, TrajError

try:
    import biggles
except:
    biggles = 0


def traj2ensemble( traj, members=10 ):
    """
    Cast normal Trajectory into EnsembleTrajectory. Resort frames
    alphabetically by their frame name. This must lead to a sorting first by
    time step and then by member.

    :param members: number of member trajectories in EnsembleTraj (default: 10)
    :type  members: int
    
    :return: converted trajectory
    :rtype: EnsembleTraj
    """
    result = EnsembleTraj( n_members=members )

    result.__dict__.update( traj.__dict__ )

    result.sortFrames()

    return result


class EnsembleTrajError( TrajError ):
    pass

class EnsembleTraj( Trajectory ):
    """
    Multiple copy trajectory. The order of the frame names is assumed to be
    as follows::
      time0_member1, time0_member2, .., time1_member1, time1_member2, ..
      i.e. frames must be sorted first by time then by member. 
    """

    def __init__(self, pdbs=None, refpdb=None, n_members=10, **kw ):
        """
        Create EnsembleTraj from PDB files or pickled PDBModels. The file names
        are used for alphabetical sorting and must lead to an ordering first
        by time then by member (i.e. frame_10_01.pdb where 10 is the time step
        and 01 is the member index). Unlike normal alphabetical order,
        frame_2_2 sorts before frame_10_10, i.e. we attempt to interprete
        numbers contained in the frame name strings.

        :param pdbs: file names of all conformations OR PDBModels
        :type  pdbs: [ str ] OR [ PDBModel ]
        :param refpdb: file name of reference pdb
        :type  refpdb: str
        :param n_members: number or tjajectories in the ensemble
        :type  n_members: int
        :param kw: optional key=value pairs, see :class:`Biskit.Trajectory`.
        :type  kw: key=value

        :raise EnsembleTrajError: if member trajectories don't have
                                  equal number of frame
        """

        Trajectory.__init__( self, pdbs, refpdb, **kw )

        self.n_members = n_members

        if self.frameNames is not None:
            self.sortFrames()

        if self.n_members and self.frames is not None and \
               len( self ) % self.n_members != 0:
            raise EnsembleTrajError('Member trajectories must have equal number of frames.')


    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.

        :param traj: trajectory
        :type  traj: trajectory
        
        Note: Overrides Trajectory method.
        """
        Trajectory.replaceContent( self, traj )
        self.n_members = traj.n_members


    def thin( self, step=1 ):
        """
        Keep only each step'th frame from trajectory with 10 ensemble members.
        
        :param step: 1..keep all frames, 2..skip first and every second, ..
                     (default: 1)
        :type  step: int
        
        :return: reduced EnsembleTraj
        :rtype: EnsembleTraj
        """
        T.ensure( step, int, forbidden=[0] )

        ## 10 x lenFrames/10, frame indices of each member
        mI = [ self.memberIndices( i ) for i in range(self.n_members) ]

        mI = N0.array( mI )

        mI = N0.take( mI, range( -1, N0.shape( mI )[1], step )[1:], 1 )

        mI = N0.transpose( mI )

        return self.takeFrames( N0.ravel( mI ))


    def memberList( self, step = 1 ):
        """
        Extract ensemble member Trajectories.
        
        @todo: cast member Trajectories back to normal Trajectory object
        
        :param step: leave out each step-1 frame (for each ensemble member)
                     (default: 1)
        :type  step: int

        :return: [EnsembleTraj] with only 1 member each
                 i.e. [ traj_member1, traj_member2, .. traj_member10 ]
        :rtype: [trajectory]
        """
        ## set of n_members trajectories, one for each ensemble member
        tm = [self.takeFrames( range(i, self.lenFrames(),self.n_members*step ))
               for i in range(0,self.n_members) ]

        for t in tm:
            t.n_members = 1

        return tm


    def memberIndices( self, member, step=1 ):
        """
        List of frame indices for this member::
          memberIndices( int_member, [int_step] )

        :param member: member trajectory
        :type  member: int
        :param step: return only every i'th frame (default: 1)
        :type  step: int

        :return: indices for members
        :rtype: [int]
        """
        r = list(range( member, self.lenFrames(), self.n_members))
        if step != 1:
            r = N0.take( r, range( 0, len( r ), step ) ).tolist()
        return r


    def memberMask( self, member ):
        """
        Get mask for all frames belonging to a given ensemble member.
        
        :param member: member index starting with 0
        :type  member: int
        
        :return: member mask, N0.array( N_frames x 1) of 1||0
        :rtype: [1|0]
        """
        result = N0.zeros( self.lenFrames() )

        if isinstance( member, (int, N0.integer)):
            N0.put( result, self.memberIndices( member ), 1 )

        if type( member ) == list:
            for m in member:
                N0.put( result, self.memberIndices( m ), 1 )

        return result


    def takeFrames( self, indices ):
        """
        Return a copy of the trajectory containing only the specified frames.
        
        :param indices: positions to take
        :type  indices: [int]
        
        :return: copy of this Trajectory (fewer frames, semi-deep copy of ref)
        :rtype: Trajectory
        """
        r = Trajectory.takeFrames( self, indices )

        r.n_members = self.n_members
        return r


    def takeMember( self, i ):
        """
        Take all frames belonging to member::
          takeMember( int ) -> EnsembleTraj, single member trajectory

        :param i: member trajectory
        :type  i: int

        :return: trajectory containing member i
        :rtype: Trajectory        
        """
        r = self.takeFrames( self.memberIndices( i ) )
        r.n_members = 1
        return r


    def takeMembers( self, mIndices ):
        """
        Take all frames belonging to the members in mIndices::
          takeMembers( mIndices ) -> EnsembleTraj with frames of given members
        
        :param mIndices: list of member indices
        :type  mIndices: [int] OR array('i')
        
        :return: EnsembleTraj with specified members
        :rtype: EnsembleTraj
        
        @todo: return self.__class__ instead of EnsembleTraj
        """
        try:
            ## assumes that each member traj has same number of frames
            fi = N0.array( [ self.memberIndices( i ) for i in mIndices ] )
            fi = N0.ravel( N0.transpose( fi ) )

            n_members = len( mIndices )

            ## has wrong n_members and member order
            t = self.takeFrames( fi )

            result = EnsembleTraj( n_members=n_members )

            result.__dict__.update( t.__dict__ )
            result.n_members = n_members
            result.resetFrameNames()

            return result

        except TypeError:
            raise EnsembleTrajError('takeMembers TypeError '+\
                  str(mIndices)+\
                  "\nlenFrames: %i; n_members: %i" %(len(self), self.n_members))


    def concatEnsembles( self, *traj ):
        """
        Concatenate this with other trajectories in a zig zac manner,
        resulting in an ensembleTraj with additional members.
        The ref model of the new Trajectory is a 'semi-deep' copy of this
        trajectorie's model.(see :class:`PDBModel.take()` )::
          concat( traj [, traj2, traj3, ..] ) -> Trajectory
        
        :param traj: with identical atoms as this one
        :type  traj: one or more EnsembleTrajectory

        @todo: fix so that pc, and profiles are not lost
        """
        if len( traj ) == 0:
            return self

        r = self.__class__( n_members = self.n_members + traj[0].n_members )

        min_members = min( self.n_members, traj[0].n_members )
        min_frames = min( self.lenFrames(), traj[0].lenFrames() )

        steps = self.lenFrames()//self.n_members + \
                traj[0].lenFrames()//traj[0].n_members

        def __everyOther( traj_0, traj_1, list_0, list_1, minMembers,
                          minFrames, loops ):
            result = []
            for j in range( 0, minMembers/2 ):

                for i in range( j*loops , j*loops + minFrames*2/minMembers ):
                    result += [ list_0[i] ]
                    result += [ list_1[i] ]

                while i < j*traj_0.n_members:
                    result += [ list_0[i] ]

                while i < j*traj_1.n_members:
                    result += [ list_1[i] ]

            return result

        frames = __everyOther( self, traj[0], self.frames,
                               traj[0].frames, min_members,
                               min_frames, steps )

        r.frames = N0.array(frames) 
        r.setRef( self.ref.clone())

        if self.frameNames and traj[0].frameNames:
            r.frameNames =  __everyOther( self, traj[0], self.frameNames,
                                          traj[0].frameNames, min_members,
                                          min_frames, steps )
        try:
            # NOT TESTED!!
            if self.pc and traj[0].pc:
                r.pc['p'] =  __everyOther( self, traj[0], self.pc['p'],
                               traj[0].pc['p'], min_members, steps )

                r.pc['u'] =  __everyOther( self, traj[0], self.pc['u'],
                               traj[0].pc['u'], min_members, steps )

#                r.pc['p'] = N0.concatenate( (self.pc['p'], traj[0].pc['p']),0)
#                r.pc['u'] = N0.concatenate( (self.pc['u'], traj[0].pc['u']),0)
        except TypeError as why:
            EHandler.error('cannot concat PC '+str(why) )

#        r.profiles = self.profiles.concat( traj[0].profiles )

        ## recursively add other trajectories
        return r.concat( *traj[1:] )


    def compressMembers( self, mask ):
        """
        Apply mask to member trajectories.
        
        :param mask: positions in trajectory list to keep or remove
        :type  mask: [1|0]
        
        :return: compressed EnsembleTraj 
        :rtype: EnsembleTraj
        """
        return self.takeMembers( N0.nonzero( mask ) )


    def keepMembers( self, indices ):
        """
        in-place version of takeMembers. keepMembers( indices ) -> None
        
        :param indices: member numbers
        :type  indices: [int]        
        """
        r = self.takeMembers( indices )
        self.replaceContent( r )


    def removeMembers( self, indices ):
        """
        Remove given member trajectories from this ensemble.
        
        :param indices: trajectory (member) numbers
        :type  indices: [int]
        """
        i = list(range( self.n_members))
        i.remove( N0.array(indices) )
        self.keepMembers( i )


    def resetFrameNames( self ):
        """
        Reset frame names to t_00_m_00 .. t_50_m_10.
        """
        n = self.n_members
        steps = self.lenFrames() // n

        time_member = []
        for s in range( steps ):
            for m in range( n ):
                time_member += [ (s,m) ]

        self.frameNames = [ 't_%03i_m_%02i' % tm for tm in time_member ]


    def argsortMember( self, inverse_time=0, inverse_member=0, step=1 ):
        """
        Get list of frame indices sorted first by member, then by time (unlike
        the normal sorting of EnsembleTraj, which is first by time then by
        member). Default is ascending order (first frame is first time step
        of member 0).
        
        :param inverse_time: descending time order (last frame first)
                             (default: 0)
        :type  inverse_time: 1|0
        :param inverse_member: descending member order (last member first)
                               (default: 0)
        :type  inverse_member: 1|0
        :param step: take only every step frame (default: 1)
        :type  step: int
        
        :return: length is N_frames (only if step is 1 )
        :rtype: [int]
        """
        r = []

        a, e, j = 0, self.n_members, 1
        if inverse_member:
            a, e, j = self.n_members, 0, -1

        for i in range( a, e, j ):

            mi = self.memberIndices( i, step )

            if inverse_time: mi.reverse()

            r.extend( mi )

        return r


    def plotMemberProfiles( self, *name, **arg ):
        """
        Plot profiles of all member trajectories seperately::
          plotMemberProfiles( name1, [name2, .. ],[ arg1=x,..])
            -> biggles.Table

        :param name: profile name(s)
        :type  name: str
        :param arg: pairs for biggles.Curve() and/or xlabel=..,ylabel=..
        :type  arg: key=value

        :return: biggles plot object
        :rtype: biggles.FramedArray()   
        """
        if not biggles:
            raise ImportError('biggles module could not be imported.')
        
        rows = self.n_members // 2 + self.n_members % 2
        page = biggles.FramedArray( rows , 2)

        biggles.configure('fontsize_min', 1)
        colors = T.colorSpectrum( len( name ) , '00FF00', 'FF00FF') 

        ml = self.memberList()

        i = 0
        minV = maxV = None

        for t in ml:

            for j in range( len(name)):

                p = t.profile( name[j] )

                if minV is None or minV > min( p ):
                    minV = min(p)
                if maxV is None or maxV < max( p ):
                    maxV = max(p)

                page[i//2, i%2].add( biggles.Curve( list(range( len(p))), list(p),
                                                   color=colors[j], **arg ) )

                page[i//2, i%2].add( biggles.PlotLabel( 0.8, 0.8-j/8.0, name[j],
                                                       color=colors[j]) )

            page[i//2, i%2].add( biggles.PlotLabel( 0.1, 0.9, 'Traj %i' % i))

            i += 1

        if self.n_members % 2 != 0:
            line = biggles.Line( (0,minV), (len(p),maxV) )
            page[ self.n_members//2, 1 ].add( line )

        page.uniform_limits = 1
        page.xlabel = arg.get('xlabel',None)
        page.ylabel = arg.get('ylabel',None)

        return page


    def fitMembers( self, refIndex=None, refModel=None, mask=None, n_it=1,
                    prof='rms', verbose=1, **profInfos ):
        """
        RMSD-fit each member trajectory seperately onto one of its frames
        or its average structure.
    
        :param refIndex: index of reference frame within member traj.
                         If None -> fit to average coordinates
                         (if refModel is None)
        :type  refIndex: int OR None
        :param refModel: fit to this structure (default: None)
        :type  refModel: PDBModel
        :param mask: atoms to consider (default: None, all heavy)
        :type  mask: [1|0] OR None
        :param n_it: number of fit iterations, kicking out outliers
                     on the way (default: 1)::
                       1 -> classic single fit,
                       0 -> until convergence
        :type  n_it: 1|0
        :param prof: save rms per frame in profile of this name (default: rms)
        :type  prof: str
        :param profInfos: description key-value pairs for profile
        :type  profInfos: key=value
        """
        ml = self.memberList()
        n = self.n_members

        for m in ml:
            if refIndex is None:
                if refModel==None:
                    ref = None
                else:
                    ref = refModel
            else:
                ref = m[ refIndex ]

            m.fit( ref=ref, mask=mask, n_it=n_it, prof=prof,
                   verbose=verbose, **profInfos )

        ## ToDo for memory efficiency: replace frame chunks on the fly
        
        self.replaceContent( ml[0].concat( *ml[1:] ) )
        self.n_members = n
        self.sortFrames()


    def blockFit( self, ref=None, mask=None ):
        """
        RMSD-fit the average of each member trajectory (i.e. the trajectory
        en block) onto the overall average (default) or a given structure.
        
        :param ref: reference structure (default: average structure)
        :type  ref: PDBModel
        :param mask: atoms to consider (default: None, all heavy)
        :type  mask: [1|0] OR None
        """
        ref = ref or self.avgModel()

        for m in range( self.n_members ):

            indices = self.memberIndices( m )

            ## get a copy of this member's Trajectory
            traj = self.takeFrames( indices )
            m_avg = traj.avgModel()

            r, t = m_avg.transformation( ref, mask )
            traj.transform( r, t )

            ## replace original frames of this member
            N0.put( self.frames, indices, traj.frames )



    def outliers( self, z=1.0, mask=None, prof='rmsCA_last', 
                  last=10, step=1, verbose=1  ):
        """
        Identify outlier trajectories. First we calculate the CA-RMS of every
        |step|th frame to the last frame. Outliers are member trajectories for
        which the **slope** of this rms profile is z standard deviations below
        the mean of all members.
        
        :param z: z-value threshold
        :type  z: float
        :param mask: atom mask used (default: ref.maskCA())
        :type  mask: [int]
        :param prof: name of pre-calculated profile to use
                     (default: 'rmsCA_last')
        :type  prof: str
        :param last: skip |last| last frames from linear regression
        :type  last: int
        :param step: frame offset
        :type  step: int
        
        :return: member mask of outlier trajectories
        :rtype: [0|1]
        """
        if mask is None: mask = self.ref.maskCA()

        traj = self.compressAtoms( mask )
        if step != 1:
            traj = traj.thin( step )

        if not prof in traj.profiles:
            traj.fitMembers( refIndex=-1, prof=prof, verbose=verbose )

        p_all = traj.profiles[ prof ]
        n = traj.n_members
        l = len( traj )

        pm = [ p_all[ member : l : n ][:-last] for member in range( n ) ]

        slopes = [ M.linfit( list(range( l//n - last)), p )[0] for p in pm ]

        mean, sd = N0.average( slopes ), M.SD( slopes )

        return [ r - mean < - z * sd for r in slopes ]


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """EnsembleTraj test"""

    def prepare(self):
        self.tr = T.load( T.testRoot() + '/lig_pcr_00/traj.dat')
        self.t50 = T.load( T.testRoot() + 'md/traj_ens50_peptide.dat')
                
    def test_EnsembleTraj( self ):
        """EnsembleTraj.fit/fitMembers/plotMembers test """
        ## The second part of the test will fail with the slimmed
        ## down test trajectory of T.testRoot(). To run the full
        ## test pease select a larger trajectory.    

        self.tr = traj2ensemble( self.tr )

        mask = self.tr.memberMask( 1 )

        self.tr.fit( ref=self.tr.ref,
                     mask=self.tr.ref.maskCA(),
                     prof='rms_CA_ref',
                     verbose=self.local )

        self.tr.fitMembers( mask=self.tr.ref.maskCA(),
                            prof='rms_CA_0', refIndex=0,
                            verbose=self.local )
        
        self.tr.fitMembers( mask=self.tr.ref.maskCA(),
                            prof='rms_CA_av',
                            verbose=self.local )

        if biggles:
            
            self.p = self.tr.plotMemberProfiles( 'rms_CA_av', 'rms_CA_0',
                                            'rms_CA_ref', xlabel='frame' )
            if self.local or self.VERBOSITY > 2:
                self.p.show()

        self.assertAlmostEqual( 26.19851,
                                 N0.sum( self.tr.profile('rms_CA_av') ), 2 )

    def test_outliers(self, traj=None):
        """EnsembleTraj.outliers/concat test"""

        self.t2 = self.tr.concat( self.tr )

        self.o = self.t2.outliers( z=1.2, mask=self.tr.ref.maskCA(),
                              verbose=self.local )
        if self.local:
            print(self.o)

        self.t = self.t2.compressMembers( N0.logical_not( self.o ) )
        
        if biggles:
            self.p2 = self.t.plotMemberProfiles( 'rms', xlabel='frame' )
            
            if self.local or self.VERBOSITY > 2:
                self.p2.show()

        self.assertEqual( self.o, 10 * [False] )

    def test_fit50members(self):
        """EnsembleTraj.fitMembers test"""
        self.t50.fitMembers( mask=self.t50.ref.maskCA(), prof='rmsCA_avg', verbose=0)
        self.assertEqual(self.t50.n_members, 50)
        self.assertTrue( 'rmsCA_avg' in self.t50.profiles )


if __name__ == '__main__':
   
    BT.localTest()
