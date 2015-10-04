## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
##
## last $Author$
## last $Date$
## $Revision$

"""
Parallize calculation of pairwise rmsd between the frames of a trajectory
"""
import Biskit.PVM.hosts as hosts

import Biskit.tools as T
import Biskit.settings as settings
import Biskit.mathUtils as mathUtils
from Biskit.Errors import BiskitError
from Biskit.LogFile import ErrLog, LogFile
from Biskit.EnsembleTraj import EnsembleTraj, traj2ensemble

import tempfile
import numpy as N
import os

## PVM imports
## from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.PVM import TrackingJobMaster


class FlexError( BiskitError ):
    pass


class TrajFlexMaster(TrackingJobMaster):
    """
    Parallize calculation of pairwise rmsd between the frames of one or two
    trajectories.
    """

    slave_script =  T.projectRoot() + '/Biskit/TrajFlexSlave.py'

    def __init__(self, traj1, traj2=None, hosts=hosts.cpus_all,
                 niceness=hosts.nice_dic, show_output=0, add_hosts=0,
                 log=None, slaveLog=None, verbose=1,
                 only_off_diagonal=1, only_cross_member=0):
        """
        @param traj1: Trajectory or EnsembleTraj, traj1 and 2 must have the
                       same atom content. If only traj1 is given, the pairwise
                       rms is calculated between its frames.
        @type  traj1: Trajectory OR EnsembleTraj
        @param traj2: see traj1
        @type  traj2: Trajectory OR EnsembleTraj
##        @param aMask: atom mask, consider only subset of atoms (default: all)
##        @type  aMask: [0|1]
        @param hosts: slave hosts to be used
                      (default: L{Biskit.PVM.hosts.cpus_all})
        @type  hosts: [str]
        @param niceness: { str:int, 'default':int }, nice value for each
                         host 
        @type  niceness: dict
        @param show_output: open xterm window for each slave (default: 0)
        @type  show_output: 0|1
        @param add_hosts: add hosts to PVM before starting (default: 0)
        @type  add_hosts: 1|0
        @param log: log file for master (default: None-> StdErr )
        @type  log: LogFile
        @param slaveLog: slave log (default: None->'TrajFlexSlave_errors.log')
        @type  slaveLog: LogFile
        @param verbose: print progress infos (default: 1)
        @type  verbose: 0|1
        @param only_off_diagonal: Don't calculate self-rms of frames.
                                  Only considered for a single trajectory
                                  (default: 1)
        @type  only_off_diagonal: 0|1
        @param only_cross_member: Don't calculate rms between frames from
                                  same member trajectory only considered for
                                  a single trajectory(requires EnsembleTraj)
                                  (default: 0)
        @type  only_cross_member: 0|1
        """
        ## create temporary folder accessible to all slaves
        self.outFolder = tempfile.mktemp('trajFlex_',
                                         dir=settings.tempDirShared )
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
                                   show_output=show_output, verbose=verbose,
                                   add_hosts=add_hosts)


    def getInitParameters(self, slave_tid):
        """
        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}
        """
        return {'ferror':self.slaveLog.fname,
                'trajMap':self.trajMap,
                'only_off_diagonal':self.only_off_diagonal,
                'only_cross_member':self.only_cross_member }


    def __windowSize( self, n_per_node, n_nodes, n_frames ):
        """
        @param n_per_node: how many chunks should be generated per node
        @type  n_per_node: int
        @param n_nodes: number of slave nodes
        @type  n_nodes: int
        @param n_frames: length of trajectory
        @type  n_frames: int

        @return: calculate number of frames per chunk
        @rtype: int
        """
        r = int(round( n_frames * 1.0 / N.sqrt(n_per_node * n_nodes) ))
        if r > 25:
            return r

        return 25


    def cleanup( self ):
        """
        Tidy up.
        """
        T.tryRemove( self.outFolder, verbose=self.verbose, tree=1 )


    def __getFrameWindows( self, traj, n_frames ):
        """
        Divide frame indices into chunks.

        @param n_frames: number of frames per chunk
        @type  n_frames: int

        @return: list with start and stop frame index of each chunk
        @rtype: [(int, int)]
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
        @param f_frames_1: file name of chunk 1 of frames
        @type  f_frames_1: {(int, int) : str}
        @param f_frames_2: file name of chunk 2 of frames
        @type  f_frames_2: {(int, int) : str}

        @return: { ((int, int),(int, int)) : (str, str) }
        @rtype:  {((int, int),(int, int)) : (str, str)}
        """
        intra_traj = f_frames_2 is None
        if intra_traj:
            if self.verbose:
                self.log.add('Intra-trajectory calculation requested.')
            f_frames_2 = f_frames_1

        if self.verbose: self.log.write('setting up task list...')

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
        @param traj: Trajectory
        @type  traj: Trajectory
        @param outFolder: folder for pickled arrays
        @type  outFolder: str
        @param prefix: file name prefix
        @type  prefix: str
        @return: { (int,int) : str } OR None, if traj is None
        @rtype: {(int,int) : str}
        """
        if traj is None:
            return None

        if self.verbose: self.log.write('dumping frame chunks...')

        n_frames = self.__windowSize( 20, len( self.hosts ), len( traj ) )

        i_windows = self.__getFrameWindows( traj, n_frames )

        r = {}

        for i in range( len(i_windows) ):

            w = i_windows[i]

            a = traj.frames[ w[0]:w[1] ]
            f = outFolder + '/%s_%i_to_%i.dat' % ((prefix,) + w)
            T.dump( a, f )
            r[w] = f

            if self.verbose and i % (len(i_windows)/50 + 1) == 0:
                self.log.write('#')

        if self.verbose: self.log.add('done')

        return r


    def memberMap(self, traj):
        """
        Tell which traj frame belongs to which member trajectory.

        @param traj: Trajectory
        @type  traj: Trajectory

        @return: member index of each frame OR None if traj is
                 not a EnsembleTraj
        @rtype: [ int ] OR None
        """
        if not isinstance( traj, EnsembleTraj ):
            return None

        r = N.zeros( len(traj), N.int )

        for i in range( traj.n_members ):

            mi = traj.memberIndices( i )
            N.put( r, mi, i )

        return r.tolist()


    def getResult( self, mirror=0 ):
        """
        Get result matrix ordered such as input trajectory.

        @param mirror: mirror the matrix at diagonal (default: 1)
                       (only for intra-traj)
        @type  mirror: 1|0

        @return: array( (n_frames, n_frames), 'f'), matrix of pairwise rms
        @rtype: array
        """
        if self.verbose:   self.log.write('assembling result matrix...')

        intra_traj = self.traj_2 is None

        n1 = n2 = len( self.traj_1 )
        if self.traj_2 is not None:
            n2 = len( self.traj_2 )

        a  = N.zeros( (n1,n2), N.float32 )

        if self.verbose: self.log.write('#')

        for key, value in self.result.items():
            i_start, i_stop = key[0]
            j_start, j_stop = key[1]

            window = N.reshape( value, (i_stop-i_start, j_stop-j_start) )
            window = window.astype(N.float32)

            a[i_start:i_stop, j_start:j_stop] = window

        if self.verbose: self.log.write('#')

        if intra_traj:
            for i in range( N.shape(a)[0] ):
                for j in range( i, N.shape(a)[1] ):
                    if a[j,i] == 0:
                        a[j,i] = a[i,j]
                    else:
                        a[i,j] = a[j,i]

        if self.verbose: self.log.write('#')

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

        @param mirror: mirror matrix at diagonal (only for intra-traj. rms)
                       (default: 0)
        @type  mirror: 0|1

        @param step: take only every step frame [1]
        @type  step: int
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
        @return: list of all calculated pairwise rms values
        @rtype: [float]

        @raise FlexError: if there are no results yet
        """
        r = []
        for v in self.result.values():
            r.extend( v )

        if not r:
            raise FlexError, "No results yet."

        return r


    def averageRms( self ):
        """
        @return: average pairwise rmsd and it's standard deviation
        @rtype: (float, float)

        @raise FlexError: if there are no results yet
        """
        r = self.rmsList()
        return N.average(r), mathUtils.SD(r)


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test Adaptive clustering"""

    TAGS = [ BT.PVM ]

    def test_FlexMaster(self):
        """TrajFlexMaster test"""
        from Biskit.MatrixPlot import MatrixPlot
        from numpy.oldnumeric.random_array import random

        assert len(hosts.cpus_all) > 0,\
               'Master requires at least 1 PVM node for initialisation.'

        traj_1 = T.load( T.testRoot() + '/lig_pcr_00/traj.dat' )
        traj_1 = traj2ensemble( traj_1 )

        ## create fake second trajectory by adding
        ## increasing noise to first
        frames = []
        for i in range( len( traj_1 ) ):
            f = traj_1.frames[i]
            d = N.zeros( N.shape( f ), N.float32)
            if i > 0:
                d = random( N.shape( f ) ) * ((i / 10) + 1) 
            frames += [f + d]

        traj_2 = traj_1.clone()
        traj_2.frames = frames

        master = TrajFlexMaster( traj_1, traj_2,
                                 hosts=hosts.cpus_all,
                                 show_output= self.local,
                                 add_hosts=1,
                                 log=None,
                                 slaveLog=None,
                                 verbose= self.local,
                                 only_cross_member=0 )

        r = master.calculateResult( mirror=0 )

        if self.local:
            p = MatrixPlot( r, palette='plasma2', legend=1 )
            p.show()            

if __name__ == '__main__':

    BT.localTest()
