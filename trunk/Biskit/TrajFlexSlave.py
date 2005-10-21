## Parallize calculation of pairwise rmsd between the frames of a trajectory
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## last $Date$
## $Revision$

from Biskit.PVM.dispatcher import JobSlave
import Biskit.tools as T
import Biskit.rmsFit as rmsFit
import Numeric as N
from Biskit.LogFile import ErrLog, LogFile

import os, time

class TrajFlexSlave( JobSlave ):
    """
    Calculate the pairwise rmsd between frames of a trajectory.
    """

    def initialize(self, params):
        """
        expects
        {'ferror':str,
        'trajMap':[int],
        'only_off_diagonal':1|0,
        'only_cross_member':1|0}
        """

        self.__dict__.update( params )

        if not self.ferror:
            self.errorLog = ErrLog()
        else:
            self.errorLog = LogFile( self.ferror, mode='a' )

        self.frame_cache = {}


    def __getFrames(self, f ):
        """Load coordinate frames from file or take them from own cache"""
        if f in self.frame_cache:
            return self.frame_cache[ f ]

        self.frame_cache[f] = T.Load(f)

        return self.frame_cache[ f ]


    def requested( self, i, j ):
        """
        requested( int_f1, int_f2 ) -> 1|0
        -> 1, if the rms of the two frames is supposed to be calculated
        """
        if self.only_off_diagonal and i == j:
            return 0

        if self.trajMap is None or not self.only_cross_member:
            return 1

        return self.trajMap[i] != self.trajMap[j]


    def calcRmsd( self, window, f1, f2 ):
        """
        window - ((int, int),(int,int)), start and end of two frame chunks
                 within the whole trajectory
        f1, f2 - array, the two frame chunks
        -> [ float ], the rms between the frames
        """
        try:
            i_start, i_stop = window[0]
            j_start, j_stop = window[1]
            
            a = N.zeros( (i_stop-i_start, j_stop-j_start), 'f' )

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
        jobs - { ((int,int),(int,int)) : (str, str) }, maps start and end
               position of two chunks of coordinate frames to the files
               where the two chunks are pickled.
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

if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:
        
        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = TrajFlexSlave()
    slave.start()

