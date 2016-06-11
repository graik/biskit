##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
Add some extra functionality to JobMaster
"""

from Biskit.PVM.dispatcher import JobMaster
import pvm
from Biskit.PVM.Status import Status
import Biskit.tools as T

from threading import Thread, RLock, _RLock, Condition, _Condition
import time
import copy

class TrackingJobMaster( JobMaster ):
    """    
    This class extends JobMaster with the following extras:
      - reporting of the average time each slave spends on a job
      - automatic adding of slave computers to PVM
      - different ways to be notified of a completed calculation
      - restarting of interrupted calculations

    The calculation is performed non-blocking in a thread after a call
    to master.start().
    The end of calculation is signalled on master.lock / master.lockMsg.
    The result can then be obtained with getResult().

    Alternatively, a callback method can be registered that is called
    after the calculation finished (master.setCallback()).

    The perhaps easiest (but also least flexible) way is to instead use the
    calculateResult() method. This starts the calculation and blocks execution
    until the result is returned.

    Consider overriding cleanup(), done() and getResult().

    An interrupted calculation can be restarted from a restart file:
       - during calculation, pickle the result of getRst() to a file
       - call the script Biskit/restartPVM -i |file_name|

    Manual restart is possible as follows:
      1. pickle master.data, master.result, master.status.objects
      2. master.exit() / Exception / kill, etc.
      3. initialize master with same parameters as before
      4. unpickle and re-assign master.data, master.result,
         master.status.objects
      5. master.start()

    @note: The master sends out an exit signal to all slaves but doesn't
      wait for a response (there isn't any) and continues in the finish()
      method. Since, at the end, the same job is distributed to several slaves,
      some of them might still be running when cleanup() or done() are
      executed. The slave script must tolerate errors that, e.g., happen
      if cleanup() is called while it is running.

    @todo: try finding some solution to the problem where the master
           sends out an exit signal to all slaves but doesn't wait for
           a response (see note)
    @todo: test restart function
    @todo: restart data are not automatically saved (e.g. in intervals)
    """

    def __init__(self, data={}, chunk_size=5,
                 hosts=[], niceness={'default':20},
                 slave_script='', verbose=1,
                 show_output=0, add_hosts=1, redistribute=1 ):
        """
        @param data: dict of items to be processed
        @type  data: {str_id:any}
        @param chunk_size: number of items that are processed per job
        @type  chunk_size: int
        @param hosts: list of host-names
        @type  hosts: [str]
        @param niceness: host niceness dictionary {str_host-name: int_niceness}
        @type  niceness: {str:int}
        @param slave_script: absolute path to slave-script
        @type  slave_script: str
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        @param show_output: display one xterm per slave (default: 0)
        @type  show_output: 1|0
        @param add_hosts: add hosts to PVM before starting (default: 1)
        @type  add_hosts: 1|0
        @param redistribute: at the end, send same job out several times
                             (default: 1)
        @type  redistribute: 1|0
        """
        if add_hosts:
            if verbose: T.errWrite('adding %i hosts to pvm...' % len(hosts) )
            pvm.addHosts( hosts=hosts )
            if verbose: T.errWriteln('done')

        JobMaster.__init__( self, data, chunk_size, hosts, niceness,
                            slave_script, show_output=show_output,
                            redistribute=redistribute, verbose=verbose )

        self.progress = {}

        self.disabled_hosts = []
        self.slow_hosts = {}

        self.verbose = verbose

        ## end of calculation is signalled on lockMsg
        self.lock = RLock()
        self.lockMsg = Condition( self.lock )

        ## this method is called when everything is calculated
        self.call_done = None


    def hostnameFromTID( self, slave_tid ):
        """
        Get nickname of host from TaskID.

        @param slave_tid: slave task tid
        @type  slave_tid: int           
        """
        nickname = self.nicknameFromTID( slave_tid )
        return nickname.split('_')[0]


    def is_valid_slave( self, slave_tid ):
        """
        Override JobMaster method to disable slow nodes on the fly

        @param slave_tid: slave task tid
        @type  slave_tid: int           
        """
        return self.hostnameFromTID( slave_tid ) not in self.disabled_hosts


    def mark_slow_slaves( self, host_list, slow_factor ):
        """
        @param host_list: list of hosts
        @type  host_list: [str]
        @param slow_factor: factor describing the calculation speed of a node
        @type  slow_factor: float
        """
        for h in host_list:
            self.slow_hosts[h] = slow_factor


    def start_job( self, slave_tid ):
        """
        Overriding JobMaster method

        @param slave_tid: slave task tid
        @type  slave_tid: int           
        """
        host = self.nicknameFromTID( slave_tid )

        d = {'given':0, 'done':0, 'time':0 }
        if self.progress.has_key( host ):
            d = self.progress[ host ]

        d['given'] += 1
        d['timeStart'] = time.time()

        self.progress[ host ] = d


    def job_done( self, slave_tid, result ):
        """
        Overriding JobMaster method

        @param slave_tid: slave task tid
        @type  slave_tid: int
        @param result: slave result dictionary
        @type  result: dict
        """
        host = self.nicknameFromTID( slave_tid )

        self.progress[host]['done'] += 1
        self.progress[host]['time'] = time.time() \
            - self.progress[host]['timeStart']


    def reportProgress( self ):
        """
        Report how many jobs were processed in what time per host.
        """
        if self.verbose:
            print 'host                     \tgiven\tdone\t  time'
            for host in self.progress:

                d = self.progress[host]
                print '%-25s\t%i\t%i\t%6.2f s' %\
                      (host, d['given'], d['done'], d['time'])


    def setCallback( self, funct ):
        """
        Register function to be called after calculation is finished.
        @param funct: will be called with an instance of the master
                      as single argument
        @type  funct: function
        """
        self.call_done = funct


    def cleanup( self ):
        """
        Called after exit. Override.
        """
        pass


    def done( self ):
        """
        Called by finish() after exit(), cleanup(), and reportProgress(), but
        before thread notification (notifyAll() ) and before executing
        the callBack method. Override.
        """
        pass


    def notifyAll( self ):
        """
        Notify thread waiting on self.lockMsg that master has finished.
        """
        self.lock.acquire()
        self.lockMsg.notifyAll()
        self.lock.release()


    def finish(self):
        """
        Called one time, after last job result has been received. It should
        not be necessary to override this further. Override done() instead.
        """
        self.exit()
        self.cleanup()

        self.reportProgress()

        self.done()

        self.notifyAll()

        if self.call_done:
            self.call_done( self )


    def getResult( self, **arg ):
        """
        Return result dict, if it is available.
        Override to return something else - which will also be the return value
        of calculateResult().

        @param arg: keyword-value pairs, for subclass implementations
        @type  arg: {key:value}

        @return: {any:any}
        @rtype: {any:any}
        """
        return self.result


    def calculateResult( self, **arg ):
        """
        Convenience function that is starting the parallel calculation and
        blocks execution until it is finished.

        @param arg: keyword-value pairs, for subclass implementations
        @type  arg: {key:value}

        @return: array( (n_frames, n_frames), 'f'), matrix of pairwise rms
        @rtype: array
        """
        self.start()

        self.lock.acquire()
        self.lockMsg.wait()
        self.lock.release()

        return self.getResult( **arg )


    def getRst( self ):
        """
        Get data necessary for a restart of the running calculation.
        Locks, file handles and private data are *NOT* saved.
        Override if necessary but call this method in child method.

        @return: {..}, dict with 'pickleable' fields of master
        @rtype: dict
        """
        self.status.lock.acquire()

        ## collect master parameters that can be pickled
        rst = {}
        for k,v in self.__dict__.items():

            skip = 0
            for t in [ Thread, _RLock, _Condition, Status, file ]:
                if isinstance( v, t ):
                    skip = 1

            if str(k)[0] == '_':
                skip = 1

            if not skip:
                rst[k] = copy.copy( v )

        rst['status_objects'] = copy.deepcopy( self.status.objects )
        rst['master_class'] = self.__class__

        self.status.lock.release()

        return rst


    def saveRst( self, fname ):
        """
        Pickle data necessary for a restart of the running calculation.

        @param fname: file name
        @type  fname: str
        """
        T.dump( self.getRst(), fname )


    def setRst( self, rst_data ):
        """
        Prepare this master for restart, called by restart().
        Override if necessary but call in child.

        @param rst_data: {..}, parameters for master.__dict__ + some
                         special fields
        @type  rst_data: dict

        @return: {..}, parameters for master.__dict__ without special fields
        @rtype: dict
        """
        self.__class__ = rst_data['master_class']
        self.status.objects = rst_data['status_objects']

        del rst_data['master_class']
        del rst_data['status_objects']

        return rst_data


def restart( rst_data, **params ):
    """
    @param rst_data: restart data
    @type  rst_data: dict
    @param params: key-value pairs, for subclass implementations
    @type  params: {key:value}
    """
    ## create empty master
    master = TrackingJobMaster( **params )

    ## switch to required subclass and handle special information
    rst_data = master.setRst( rst_data )

    ## set all remaining fields of master
    master.__dict__.update( rst_data )

    return master

if __name__ == '__main__':

    import pypvm

    try:

        pypvm.hostinfo()

    except:
        print T.lastError()
