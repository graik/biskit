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
Parallize calculation of pairwise rmsd between the frames of a trajectory.
"""

from Biskit.PVM import JobSlave
import Biskit.tools as T
import Biskit.rmsFit as rmsFit
import numpy as N
from Biskit.LogFile import ErrLog, LogFile

import os, time

class TrajFlexSlave( JobSlave ):
    """
    Calculate the pairwise rmsd between frames of a trajectory.
    """

    def initialize(self, params):
        """
        Expects::
          {'ferror':str,
           'trajMap':[int],
           'only_off_diagonal':1|0,
           'only_cross_member':1|0}
        
        @param params: parameters passed over from the L{TrajFlexMaster}
        @type  params: dict
        """

        self.__dict__.update( params )

        if not self.ferror:
            self.errorLog = ErrLog()
        else:
            self.errorLog = LogFile( self.ferror, mode='a' )

        self.frame_cache = {}


    def __getFrames(self, f ):
        """
        Load coordinate frames from file or take them from own cache

        @param f: file name
        @type  f: str

        @return: coordiante frames
        @rtype: array        
        """
        if f in self.frame_cache:
            return self.frame_cache[ f ]

        self.frame_cache[f] = T.load(f)

        return self.frame_cache[ f ]


    def requested( self, i, j ):
        """
        Checks if the rmsd of two frames i and j are to be calculated::
          requested( int_f1, int_f2 ) -> 1|0

        @param i: frame number
        @type  i: int
        @param j: frame number
        @type  j: int
        
        @return: if the rms of the two frames is supposed to be calculated
        @rtype: 1|0
        """
        if self.only_off_diagonal and i == j:
            return 0

        if self.trajMap is None or not self.only_cross_member:
            return 1

        return self.trajMap[i] != self.trajMap[j]


    def calcRmsd( self, window, f1, f2 ):
        """
        Calulate the rmsd between two frame chunks.
        
        @param window: start and end of two frame chunks within the
                       whole trajectory
        @type  window: ((int, int),(int,int))
        @param f1: frame chunk
        @type  f1: array
        @param f2: frame chunk
        @type  f2: array
        
        @return: the rms between the frames
        @rtype: [float]
        """
        try:
            i_start, i_stop = window[0]
            j_start, j_stop = window[1]

            a = N.zeros( (i_stop-i_start, j_stop-j_start), N.float )

            i = j = -1

            ## block on the diagonal, only calculate one half of it
            S = (self.only_off_diagonal and window[0] == window[1])

            for i in range( i_start, i_stop ):
                for j in range( S * i - S * j_start + j_start, j_stop ):

                    if self.requested( i, j ):

                        rt, rmsdLst = rmsFit.match( f1[i-i_start],
                                                    f2[j-j_start], 1 )
                        a[i-i_start,j-j_start] = rmsdLst[0][1]

            return N.ravel(a).tolist()

        except Exception, why:
            self.reportError( 'ERROR '+str(why), (i,j) )
            return


    def go(self, jobs):
        """
        Run job.
        
        @param jobs: { ((int,int),(int,int)) : (str, str) }, maps start and end
                      position of two chunks of coordinate frames to the
                      files where the two chunks are pickled.
        @type  jobs: {((int,int),(int,int)) : (str, str)}

        @return: the rms between the frames
        @rtype: [float]              
        """

        result = {}
        frames = None

        startTime = time.time()

        try:

            for i, frames in jobs.items():

                T.flushPrint( str(i) )

                f1 = self.__getFrames( frames[0])
                f2 = self.__getFrames( frames[1])

                result[ i ] = self.calcRmsd(i, f1, f2 )

            print "\navg time for last %i complexes: %f s" %\
                  ( len(jobs), (time.time()-startTime)/len(jobs))

        except IOError, why:
            self.reportError("Cannot open temporary frame file "+\
                             "(can happen if slave catches exit signal)",
                             frames )

        return result


    def reportError(self, msg, window ):
        """
        Report errors.
        
        @param msg: error message
        @type  msg: str
        @param window: start and end of two frame chunks within the
                       whole trajectory
        @type  window: ((int, int),(int,int))        
        """
        try:
            s = '%s on %s, frames %s \n' % \
                (msg, os.uname()[1], str(window) )
            s += '\nErrorTrace:\n' + T.lastErrorTrace() + '\n'

            try:
                print s
            except:
                pass

            self.errorLog.add( s )

        except Exception, why:
            f = open('ErrorReportError_TrajFlexSlave','a')
            f.write( str(why) )
            try:
                f.write( T.lastErrorTrace() )
            except:
                pass
            f.close()

##############
## empty test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, the Slave is tested in L{Biskit.TrajFlexMaster}"""
    pass


if __name__ == '__main__':

    import sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = TrajFlexSlave()
    slave.start()

