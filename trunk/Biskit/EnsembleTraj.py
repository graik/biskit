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
"""Multiple copy trajectory"""

from Trajectory import Trajectory, TrajError
from Biskit import EHandler

import Biskit.mathUtils as M

import types
import Numeric as N
import Biskit.tools as T
import biggles

def traj2ensemble( traj, members=10 ):
    """
    Cast normal Trajectory into EnsembleTrajectory. Resort frames
    alphabetically by their frame name. This must lead to a sorting first by
    time step and then by member.
    -> EnsembleTraj
    """
    result = EnsembleTraj( n_members=members )

    result.__dict__.update( traj.__dict__ )

    result.sortFrames()

    return result

class EnsembleTrajError( TrajError ):
    pass

class EnsembleTraj( Trajectory ):
    """
    Multiple copy trajectory. The order of the frames is assumed to be
    as follows:
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
        """

        Trajectory.__init__( self, pdbs, refpdb, **kw )

        self.n_members = n_members

        if self.frameNames != None:
            self.sortFrames()

        if self.n_members and self.frames is not None and \
               len( self ) % self.n_members != 0:
            raise EnsembleTrajError,\
                  'Member trajectories must have equal number of frames.'


    def version( self ):
        return Trajectory.version(self) + '; EnsembleTraj $Revision$'

    def replaceContent( self, traj ):
        """
        Replace content of this trajectory by content of given traj.
        No deep-copying, only references are taken.
        Overrides Trajectory method.
        """
        Trajectory.replaceContent( self, traj )
        self.n_members = traj.n_members
        
    def thin( self, step=1 ):
        """
        Keep only each step'th frame from trajectory with 10 ensemble members.
        step - int, 1..keep all frames, 2..skip first and every second, ..
        -> EnsembleTraj
        """
        T.ensure( step, int, forbidden=[0] )

        ## 10 x lenFrames/10, frame indices of each member
        mI = [ self.memberIndices( i ) for i in range(self.n_members) ]

        mI = N.array( mI )

        mI = N.take( mI, range( -1, N.shape( mI )[1], step )[1:], 1 )

        mI = N.transpose( mI )
        
        return self.takeFrames( N.ravel( mI ))


    def memberList( self, step = 1 ):
        """
        Extract ensemble member Trajectories.
        step - leave out each step-1 frame (for each ensemble member)

        ToDo: cast member Trajectories back to normal Trajectory object

        -> [EnsembleTraj] i.e. [ traj_member1, traj_member2, .. traj_member10 ]
        """
        ## set of 10 trajectories, one for each ensemble member
        tm = [self.takeFrames( range(i, self.lenFrames(),self.n_members*step ))
               for i in range(0,self.n_members) ]

        return tm


    def memberIndices( self, member, step=1 ):
        """
        memberIndices( int_member, [int_step] ) -> list of frame indices for
        this member
        """
        r = range( member, self.lenFrames(), self.n_members )
        if step != 1:
            r = take( r, range( 0, len( r ), step ) ).tolist()
        return r
    

    def memberMask( self, member ):
        """
        Get mask for all frames belonging to a given ensemble member.
        member - int, member index starting with 0
        -> N.array( N_frames x 1) of 1||0
        """
        result = N.zeros( self.lenFrames() )
         
        if type( member ) == types.IntType:
            N.put( result, self.memberIndices( member ), 1 )
            
        if type( member ) == types.ListType:
            for m in member:
                N.put( result, self.memberIndices( m ), 1 )

        return result


    def takeFrames( self, indices ):
        """
        indices - list of int, positions to take
        -> copy of this Trajectory (fewer frames, semi-deep copy of ref)
        """
        r = Trajectory.takeFrames( self, indices )

        r.n_members = self.n_members
        return r
    

    def takeMember( self, i ):
        """
        takeMember( int ) -> EnsembleTraj, single member trajectory
        """
        r = self.takeFrames( self.memberIndices( i ) )
        r.n_members = 1
        return r


    def takeMembers( self, mIndices ):
        """
        takeMembers( mIndices ) -> EnsembleTraj with frames of given members
        mIndices - [int] OR array('i')
        -> EnsembleTraj
        @todo: return self.__class__ instead of EnsembleTraj
        """
        try:
            ## assumes that each member traj has same number of frames
            fi = N.array( [ self.memberIndices( i ) for i in mIndices ] )
            fi = N.ravel( N.transpose( fi ) )

            n_members = len( mIndices )

            ## has wrong n_members and member order
            t = self.takeFrames( fi )

            result = EnsembleTraj( n_members=n_members )

            result.__dict__.update( t.__dict__ )
            result.n_members = n_members
            result.resetFrameNames()

            return result

        except TypeError:
            raise EnsembleTrajError, 'takeMembers TypeError '+\
                  str(mIndices)+\
                  "\nlenFrames: %i; n_members: %i" %(len(self), self.n_members)


    def concatEnsembles( self, *traj ):
        """
        concat( traj [, traj2, traj3, ..] ) -> Trajectory
        Concatenate this with other trajectories in a zig zac manner,
        resulting in an ensembleTraj with additional members.
        The ref model of the new Trajectory is a 'semi-deep' copy of this
        trajectorie's model.(see PDBModel.take() )
        traj - one or more EnsembleTrajectory, with identical atoms as this one

        # TODO fix so that pc, and profiles are not lost
        """
        if len( traj ) == 0:
            return self
        
        r = self.__class__( n_members = self.n_members + traj[0].n_members )
        
        min_members = min( self.n_members, traj[0].n_members )
        min_frames = min( self.lenFrames(), traj[0].lenFrames() )

        steps = self.lenFrames()/self.n_members + \
                traj[0].lenFrames()/traj[0].n_members

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
        
        r.frames = N.array(frames) 
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
                
#                r.pc['p'] = N.concatenate( (self.pc['p'], traj[0].pc['p']),0)
#                r.pc['u'] = N.concatenate( (self.pc['u'], traj[0].pc['u']),0)
        except TypeError, why:
            EHandler.error('cannot concat PC '+str(why) )

#        r.profiles = self.profiles.concat( traj[0].profiles )
    
        ## recursively add other trajectories
        return r.concat( *traj[1:] )

    
    def compressMembers( self, mask ):
        """
        Apply mask to member trajectories.
        mask - [ 1|0 ], positions in trajectory list to keep or remove
        -> EnsembleTraj
        """
        return self.takeMembers( N.nonzero( mask ) )

        
    def keepMembers( self, indices ):
        """in-place version of takeMembers. keepMembers( indices ) -> None"""
        r = self.takeMembers( indices )
        self.replaceContent( r )


    def removeMembers( self, indices ):
        """
        Remove given member trajectories from this ensemble.
        indices - [ int ], trajectory numbers
        """
        i = range( self.n_members )
        i.remove( N.array(indices) )
        self.keepMembers( i )


    def resetFrameNames( self ):
        """
        Reset frame names to t_00_m_00 .. t_50_m_10.
        """
        n = self.n_members
        steps = self.lenFrames() / n

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
        inverse_time   - 1|0, descending time order (last frame first) [0]
        inverse_member - 1|0, descending member order (last member first) [0]
        step           - int, take only every step frame [1]
        -> [ int ], length is N_frames (only if step is 1 )
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
        Plot profiles of all member trajectories seperately.
        plotMemberProfiles( name1, [name2, .. ],[ arg1=x,..]) -> biggles.Table
        name - str, profile name(s)
        arg  - key=value, pairs for biggles.Curve() and/or xlabel=..,ylabel=..
        """
        rows = self.n_members / 2 + self.n_members % 2
        page = biggles.FramedArray( rows , 2)

        biggles.configure('fontsize_min', 1)
        colors = colorSpectrum( len( name ) , '00FF00', 'FF00FF') 

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
            
                page[i/2, i%2].add( biggles.Curve( range( len(p) ), p,
                                                   color=colors[j], **arg ) )

                page[i/2, i%2].add( biggles.PlotLabel( 0.8, 0.8-j/8.0, name[j],
                                                       color=colors[j]) )

            page[i/2, i%2].add( biggles.PlotLabel( 0.1, 0.9, 'Traj %i' % i))

            i += 1

        if self.n_members % 2 != 0:
            line = biggles.Line( (0,minV), (len(p),maxV) )
            page[ self.n_members/2, 1 ].add( line )
            
        page.uniform_limits = 1
        page.xlabel = arg.get('xlabel',None)
        page.ylabel = arg.get('ylabel',None)

        return page


    def fitMembers( self, refIndex=None, refModel=None, mask=None, n_it=1,
             prof='rms', **profInfos ):
        """
        RMSD-fit each member trajectory seperately onto one of its frames
        or its average structure.
        refIndex - int, index of reference frame within member traj
                   [None] -> fit to average coordinates (if refModel == None)
        refModel - PDBModel, fit to this structure [None]
        mask     - array of 1||0, atoms to consider,  [None] -> all heavy
        n_it - int, number of fit iterations, kicking out outliers on the way
               [1] -> classic single fit, 0 -> until convergence
        prof     - str, save rms per frame in profile of this name, ['rms']
        profInfos- description key-value pairs for profile, []
        """
        ml = self.memberList()

        for m in ml:
            if refIndex == None:
                if refModel==None:
                    ref = None
                else:
                    ref = refModel
            else:
                ref = m[ refIndex ]
            
            m.fit( ref=ref, mask=mask, n_it=n_it, prof=prof, **profInfos )

        ## ToDo for memory efficiency: replace frame chunks on the fly

        self.replaceContent( ml[0].concat( *ml[1:] ) )
        self.sortFrames()


    def blockFit( self, ref=None, mask=None ):
        """
        RMSD-fit the average of each member trajectory (i.e. the trajectory
        en block) onto the overall average (default) or a given structure.
        ref  - PDBModel, reference structure (default: average structure)
        mask - array of 1||0, atoms to consider,  [None] -> all heavy
        conv - float, iteratively fit / recalculate average until
               (rms_i - rms_i-1) < conv
               [None -> only 1 iteration] 1e-6 is a good value otherwise
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
            N.put( self.frames, indices, traj.frames )

            

    def outliers( self, z=1.0, mask=None, prof='rmsCA_last', 
                  last=10, step=1  ):
        """
        Identify outlier trajectories. First we calculate the CA-RMS of every
        |step|th frame to the last frame. Outliers are member trajectories for
        which the slope of this rms profile is z standard deviations below the
        mean of all members.
        z    - float, z-value threshold
        mask - [int], atom mask used                      [ ref.maskCA() ]
        prof - str, name of pre-calculated profile to use ['rmsCA_last']
        last - int, skip |last| last frames from linear regression
        step - int, frame offset
        -> [ 0|1 ], member mask of outlier trajectories
        """
        mask = mask or self.ref.maskCA()

        traj = self.compressAtoms( mask )
        if step != 1:
            traj = traj.thin( step )
        
        if not prof in traj.profiles:
            traj.fitMembers( refIndex=-1, prof=prof )

        p_all = traj.profiles[ prof ]
        n = traj.n_members
        l = len( traj )

        pm = [ p_all[ member : l : n ][:-last] for member in range( n ) ]

        slopes = [ M.linfit( range( l/n - last ), p )[0] for p in pm ]
        
        mean, sd = N.average( slopes ), M.SD( slopes )

        return [ r - mean < - z * sd for r in slopes ]


#######TEST############

if __name__ == '__main__':

    import time
    
    tr = T.Load( T.testRoot() + '/lig_pc2_00/traj.dat')

    ## The second part of the test will fail with the slimmed
    ## down test trajectory of T.testRoot(). To run the full
    ## test pease select a larger trajectory.    

    ## tr = T.Load( '~/interfaces/c11/lig_pcr_00/traj.dat')

    t0 = time.time()

    tr = traj2ensemble( tr )

    mask = tr.memberMask( 1 )

    tr.fit( ref=tr.ref, mask=tr.ref.maskCA(), prof='rms_CA_ref' )
    
    tr.fitMembers(mask=tr.ref.maskCA(), prof='rms_CA_0', refIndex=0)
    tr.fitMembers(mask=tr.ref.maskCA(), prof='rms_CA_av')
    
    p = tr.plotMemberProfiles( 'rms_CA_av', 'rms_CA_0', 'rms_CA_ref',
                               xlabel='frame' )

    p.show()

    print "Outliers..."

    o = tr.outliers( z=1.2, mask=tr.ref.maskCA() )
    print o

    t = tr.compressMembers( N.logical_not( o ) )

    p2 = t.plotMemberProfiles( 'rms_CA_av', 'rms_CA_0', 'rms_CA_ref',
                               xlabel='frame' )
    p2.show()

    print "done in ", time.time() - t0, 's'
