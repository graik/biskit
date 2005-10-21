## Parallize calculation of pairwise rmsd between the frames of a trajectory
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## last $Date$
## $Revision$

import Biskit.hosts as hosts

import Biskit.tools as T
import Biskit.settings as settings
import Biskit.mathUtils as mathUtils
from Biskit.Errors import BiskitError
from Biskit.LogFile import ErrLog, LogFile
from Biskit.EnsembleTraj import EnsembleTraj, traj2ensemble

import tempfile
import Numeric as N
import os, time

## PVM imports
from Biskit.PVM.TrackingJobMaster import TrackingJobMaster


class FlexError( BiskitError ):
    pass

class TrajFlexMaster(TrackingJobMaster):
    """
    Parallize calculation of pairwise rmsd between the frames of one or two
    trajectories.
    """

    slave_script =  T.projectRoot() + '/Biskit/TrajFlexSlave.py'

    def __init__(self, traj1, traj2=None,
                 aMask       =None,
                 hosts       =hosts.cpus_all,
                 niceness    =hosts.nice_dic,
                 show_output =0,
                 add_hosts   =0,
                 log=None, slaveLog=None, verbose=1,
                 only_off_diagonal=1, only_cross_member=0):
        """
        traj1, traj2 - Trajectory or EnsembleTraj, traj1 and 2 must have the
                       same atom content. If only traj1 is given, the pairwise
                       rms is calculated between its frames.
        aMask        - [ 0|1 ], atom mask, consider only subset of atoms [all]
        hosts        - [ str ], slave hosts to be used [Biskit.hosts.cpus_all]
        niceness     - { str:int, 'default':int }, nice value for each host []
        show_output  - 0|1, open xterm window for each slave [0]
        add_hosts    - 1|0, add hosts to PVM before starting [0]
        log          - LogFile, log file for master [None-> StdErr]
        slaveLog     - LogFile, slave log [None->'TrajFlexSlave_errors.log']
        verbose      - 0|1, print progress infos [1]
        the following two options are only considered for a single trajectory:
        only_off_diagonal - 0|1, don't calculate self-rms of frames [1]
        only_cross_member - 0|1, don't calculate rms between frames from same
                            member trajectory (requires EnsembleTraj) [0]
        """

        ## create temporary folder accessible to all slaves
        tmp = tempfile.tempdir
        tempfile.tempdir = settings.tempDirShared

        self.outFolder = tempfile.mktemp('trajFlex_')
        tempfile.tempdir = tmp
        os.mkdir( self.outFolder )

        self.log      = log or ErrLog()
        self.slaveLog = slaveLog or LogFile('TrajFlexSlave_errors.log')
        self.verbose  = verbose
        self.hosts    = hosts
        
        self.traj_1 = traj1
        self.traj_2 = traj2

        self.only_off_diagonal = traj2 is None and only_off_diagonal
        self.only_cross_member = traj2 is None and only_cross_member

        self.trajMap = None
        if traj2 is None:
            self.trajMap = self.memberMap( traj1 )

        ## pickle chunks of coordinate frames
        frame_files_1 = self.__dumpFrames( traj1, self.outFolder, 'traj1' )
        ## None if traj2 is None
        frame_files_2 = self.__dumpFrames( traj2, self.outFolder, 'traj2' )

        ## assemble job dict
        self.tasks = self.__taskDict( frame_files_1, frame_files_2)

        chunk_size = 1
        TrackingJobMaster.__init__(self, self.tasks, chunk_size,
                                   hosts, niceness, self.slave_script,
                                   show_output=show_output,
                                   add_hosts=add_hosts)


    def getInitParameters(self,  slave_tid):
        """
        hand over parameters to slave once.
        """
        return {'ferror':self.slaveLog.fname,
                'trajMap':self.trajMap,
                'only_off_diagonal':self.only_off_diagonal,
                'only_cross_member':self.only_cross_member }


    def __windowSize( self, n_per_node, n_nodes, n_frames ):
        """
        n_per_node - int, how many chunks should be generated per node
        n_nodes    - int, number of slave nodes
        n_frames   - int, length of trajectory
        -> int, calculate number of frames per chunk
        """
        r = int(round( n_frames * 1.0 / N.sqrt(n_per_node * n_nodes) ))
        if r > 25:
            return r

        return 25

    def cleanup( self ):
        T.tryRemove( self.outFolder, verbose=self.verbose, tree=1 )


    def __getFrameWindows( self, traj, n_frames ):
        """
        Divide frame indices into chunks.
        n_frames  - int, number of frames per chunk
        -> [ (int, int) ], list with start and stop frame index of each chunk
        """
        if traj is None:
            return None

        l = len( traj )
        ## number of windows
        n, rest = l / n_frames, l % n_frames
        if rest:
            n += 1

        ## get start and end frame for all windows
        i_windows = []

        for i in range(0,n):

            start, stop = i*n_frames, i*n_frames + n_frames

            if stop > l: stop = l
            if i == n-1:
                stop = l

            i_windows.append( (start, stop) )

        return i_windows
        

    def __taskDict( self, f_frames_1, f_frames_2 ):
        """
        f_frames - { (int, int) : str }, file name of each chunk of frames
        n_frames - int, jobs will be sent out in chunks of n_frames x n_frames
        -> { ((int, int),(int, int)) : (str, str) }
        """
        intra_traj = f_frames_2 is None
        if intra_traj:
            if self.verbose:
                self.log.add('Intra-trajectory calculation requested.')
            f_frames_2 = f_frames_1

        if self.verbose: self.log.add_nobreak('setting up task list...')

        i_windows = f_frames_1.keys()
        j_windows = f_frames_2.keys()

        ## assemble dict with all combinations of windows
        tasks = {}

        for i in range( len(i_windows) ):
            for j in range( i * intra_traj, len(j_windows) ):

                start_i, stop_i = i_windows[i]
                start_j, stop_j = j_windows[j]

                key = ((i_windows[i], j_windows[j]))

                tasks[key] = ( f_frames_1[ i_windows[i] ],
                               f_frames_2[ j_windows[j] ] )
                
        if self.verbose: self.log.add('done')

        return tasks


    def __dumpFrames(self, traj, outFolder, prefix ):
        """
        traj      - Trajectory
        outFolder - str, folder for pickled arrays
        nFrames   - int, number of frames per window (chunk)
        prefix    - str, file name prefix
        -> { (int,int) : str } OR None, if traj is None
        """
        if traj is None:
            return None
        
        if self.verbose: self.log.add_nobreak('dumping frame chunks...')

        n_frames = self.__windowSize( 20, len( self.hosts ), len( traj ) )

        i_windows = self.__getFrameWindows( traj, n_frames )

        r = {}

        for i in range( len(i_windows) ):

            w = i_windows[i]
            
            a = traj.frames[ w[0]:w[1] ]
            f = outFolder + '/%s_%i_to_%i.dat' % ((prefix,) + w)
            T.Dump( a, f )
            r[w] = f

            if self.verbose and i % (len(i_windows)/50 + 1) == 0:
                self.log.add_nobreak('#')

        if self.verbose: self.log.add('done')

        return r


    def memberMap(self, traj):
        """
        Tell which traj frame belongs to which member trajectory.
        -> [ int ], member index of each frame
        -> None, if traj is not a EnsembleTraj
        """
        if not isinstance( traj, EnsembleTraj ):
            return None

        r = N.zeros( len(traj), 'i' )
        
        for i in range( traj.n_members ):

            mi = traj.memberIndices( i )
            N.put( r, mi, i )

        return r.tolist()


    def getResult( self, mirror=0 ):
        """
        Get result matrix ordered such as input trajectory.
        mirror - 1|0, mirror the matrix at diagonal [1] (only for intra-traj)
        -> array( (n_frames, n_frames), 'f'), matrix of pairwise rms
        """
        if self.verbose:   self.log.add_nobreak('assembling result matrix...')
        
        intra_traj = self.traj_2 is None

        n1 = n2 = len( self.traj_1 )
        if self.traj_2 is not None:
            n2 = len( self.traj_2 )

        a  = N.zeros( (n1,n2), 'f' )

        if self.verbose: self.log.add_nobreak('#')
        
        for key, value in self.result.items():
            i_start, i_stop = key[0]
            j_start, j_stop = key[1]

            window = N.reshape( value, (i_stop-i_start, j_stop-j_start) )
            window = window.astype('f')

            a[i_start:i_stop, j_start:j_stop] = window

        if self.verbose: self.log.add_nobreak('#')

        if intra_traj:
            for i in range( N.shape(a)[0] ):
                for j in range( i, N.shape(a)[1] ):
                    if a[j,i] == 0:
                        a[j,i] = a[i,j]
                    else:
                        a[i,j] = a[j,i]

        if self.verbose: self.log.add_nobreak('#')

        if intra_traj and not mirror:
            for i in range( N.shape(a)[0] ):
                for j in range( i, N.shape(a)[1] ):
                    a[j,i] = 0.

        if self.verbose:   self.log.add('done')

        return a


    def rmsMatrixByMember( self, mirror=0, step=1 ):
        """
        Get result matrix ordered first by member then by time. (requires
        EnsembleTraj)
        mirror - 0|1, mirror matrix at diagonal (only for intra-traj. rms) [0]
        step   - int, take only every step frame [1]
        """
        intra_traj = self.traj_2 is None
        
        m = self.getResult( mirror=intra_traj )

        i1 = i2 = self.traj_1.argsortMember( step=step )

        if self.traj_2 is not None:
            i2 = self.traj_2.argsortMember( step=step )

        a = N.take( m, i1, 0 )
        a = N.take( a, i2, 1 )

        if intra_traj and not mirror:
            for i in range( N.shape(a)[0] ):
                for j in range( i, N.shape(a)[1] ):
                    a[j,i] = 0.

        return a


    def rmsList( self ):
        """
        -> [float], list of all calculated pairwise rms values
        !! FlexError, if there are no results yet
        """
        r = []
        for v in self.result.values():
            r.extend( v )
            
        if not r:
            raise FlexError, "No results yet."

        return r

          
    def averageRms( self ):
        """
        -> (float, float), average pairwise rmsd and it's standard deviation
        !! FlexError, if there are no results yet
        """
        r = self.rmsList()
        return N.average(r), mathUtils.SD(r)


if __name__ == '__main__':

    from Biskit.MatrixPlot import MatrixPlot
    from RandomArray import random

    traj_1 = T.Load( T.testRoot() + '/lig_pc2_00/traj.dat' )
    traj_1 = traj2ensemble( traj_1 )

    ## create fake second trajectory by adding noise to first
    frames = []
    for i in range( len( traj_1 ) ):
        f = traj_1.frames[i]
        d = N.zeros( N.shape( f ), 'f')
        if i > 0:
            d = random( N.shape( f ) ) * ((i / 10) + 1) 
        frames += [f + d]

    traj_2 = traj_1.clone()
    traj_2.frames = frames


    ## Bigger example (use all nodes!):
##     T.flushPrint("Loading...")
##     traj = T.Load('~/interfaces/a11/com_pcr_00/traj_0.dat')
##     traj = traj.thin( step=50 )
##     traj = traj2ensemble( traj )
##     T.flushPrint("done\n")

    master = TrajFlexMaster( traj_1, traj_2,
                             hosts=hosts.cpus_all,
                             show_output=0,
                             add_hosts=1,
                             log=None, slaveLog=None,
                             verbose=0,
                             only_cross_member=0)

    r = master.calculateResult( mirror=0 )

    p = MatrixPlot( r, palette='plasma2', legend=1 )
    p.show()
