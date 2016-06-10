##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2002-2005 Wolfgang Rieping
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

## Contributions: Wolfgang Rieping, Raik Gruenberg

"""
Manage Master/Slave tasks.
"""


from PVMThread import PVMMasterSlave
from Biskit import ExeConfigCache
from Status import Status
import Biskit.settings as settings
import socket, pvm
import pypvm

MSG_JOB_START = 1
MSG_JOB_DONE = 2

class JobMaster(PVMMasterSlave):

    def __init__(self, data, chunk_size, hosts, niceness, slave_script,
                 show_output = 0, result = None, redistribute=1, verbose=1 ):
        """
        @param data: dict of items to be proessed {id:object}.
        @type  data: dict
        @param chunk_size: number of items that are processed by a job
        @type  chunk_size: int
        @param hosts: list of host-names
        @type  hosts: [str]
        @param niceness: nice dictionary [host-name: niceness]
        @type  niceness: dict
        @param slave_script: absolute path to slave-script
        @type  slave_script: str
        @param result: items which have already been processed
                       (ie. they are contained in result) are not
                       processed again.
        @type  result: dict
        @param redistribute: at the end, send same job out several times
                             (default: 1)
        @type  redistribute: 1|0
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        """
        PVMMasterSlave.__init__(self, verbose=verbose)

        ## change names of multiple hosts
        d = {}

        for host in hosts:
            if host in d:
                d[host] += 1
            else:
                d[host] = 1

        unique_list = []

        for host, number in d.items():
            if number > 1:

                for i in range(number):
                    nickname = host + '_%d' % i

                    d = {'host': host,
                         'nickname': nickname}

                    unique_list.append(d)
            else:
                d = {'host': host,
                     'nickname': host}

                unique_list.append(d)

        exe = ExeConfigCache.get('xterm')
        exe.validate()
        self.xterm_bin = exe.bin
        
        self.hosts = unique_list
        self.niceness = niceness
        self.data = data
        self.slave_script = slave_script
        self.chunk_size = chunk_size

        self.current_pos = 0
        self.show_output = show_output

        self.verbose = verbose
        
        if result is None:
            result = {}

        self.result = result

        ## set-up status for results
        items = []

        for key in data.keys():
            if not key in self.result:
                items.append(key)

        self.status = Status(items, redistribute=redistribute )

        self.__finished = 0

        if verbose: print 'Processing %d items ...' % len(items)


    def start(self):
        """
        Start slave job
        """
        assert len(self.hosts) > 0,\
               'Master needs at least 1 pvm node to start calculations.'
        self.finished = 0

        PVMMasterSlave.start(self)

        self.startMessageLoop()
        self.spawnAll(self.niceness, self.show_output)


    def getInitParameters(self, slave_tid):
        """
        Override to collect slave initiation parameters.
        """
        return None


    def done(self):
        """
        Override to do something after last job result has been received.
        """
        print "Done"


    def finish( self ):
        """
        This method is called one time, after the master has
        received the last missing result.
        """
        self.done()


    def __finish(self):
        """
        Call finish(); but only call it once.
        """
        self.status.lock.acquire()

        if not self.__finished:
            self.finish()
            self.__finished = 1

        self.status.lock.release()


    def bindMessages(self, slave_tid):
        """
        @param slave_tid: slave task tid
        @type  slave_tid: int
        """
        self.bind(MSG_JOB_DONE, slave_tid, self.__job_done)


    def spawn(self, host, nickname, niceness, show_output = 0):
        """
        Spawn a job.

        @param host: host name
        @type  host: str
        @param nickname: host nickname (for uniqueness, i.e more than
                         on job vill be run on a multiple cpu machine)
        @type  nickname: str
        @param niceness: nice dictionary [host-name: niceness]
        @type  niceness: int

        @return: slave task tid
        @rtype: int
        """
        if show_output:
            display = socket.gethostname() + ':0.0'

            command = self.xterm_bin
            argv = ['-title', str(host), '-geometry', '30x10',
                    '-display', display, '-e',
                    settings.python_bin, '-i', self.slave_script,
                    str(niceness)]
        else:
            command = settings.python_bin
            argv = ['-i', self.slave_script, str(niceness)]

        args = (command, argv, pypvm.spawnOpts['TaskHost'], host, 1)

        return PVMMasterSlave.spawn(self, args, nickname)


    def spawnAll(self, niceness, show_output = 0):
        """
        Spawn many jobs.

        @param niceness: nice dictionary [host-name: niceness]
        @type  niceness: dict
        """
        self.slaves = {}

        for d in self.hosts:
            host = d['host']
            nickname = d['nickname']

            try:
                nice = niceness[host]
            except:
                nice = niceness.get('default', 0)

            slave_tid = self.spawn(host, nickname, nice, show_output)

            if slave_tid <= 0:
                print 'error spawning', host
                try:
                    print '\t', pvm.pvmerrors[ slave_tid ]
                except Exception, error:
                    print 'unknown error', error

            else:
                self.bindMessages(slave_tid)
                self.slaves[slave_tid] = d
                if self.verbose: print slave_tid, nickname, 'spawned.'


    def initializationDone(self, slave_tid):
        """
        is called by a slave that has been initialized and
        is now ready for start-up.

        @param slave_tid: slave task tid
        @type  slave_tid: int        
        """
        ## start processing
        self.__start_job(slave_tid)


    def start_job(self, slave_tid):
        """
        Called when a new job is about to be started. Override to add
        other startup tasks than the default, see L{__start_job}

        @param slave_tid: slave task tid
        @type  slave_tid: int
        """
        pass


    def get_slave_chunk(self, data_keys ):
        """
        Assemble task dictionary that is send to the slave for a single job.
        Override this, if the values of self.data are to be changed/created on
        the fly.
        
        @param data_keys: subset of keys to self.data
        @type  data_keys: [any]
        
        @return: dict mapping the data keys to the actual data values
        @rtype: {any:any}
        """
        chunk = {}

        for id in data_keys:
            chunk[id] = self.data[id]

        return chunk


    def __start_job(self, slave_tid):
        """
        Tasks performed befor the job is launched.
        
        @param slave_tid: slave task tid
        @type  slave_tid: int
        """
        self.start_job(slave_tid)

        ## get items that have not been processed
        queue, n_left = self.status.next_chunk( self.chunk_size)

        if not queue:
            return

        nickname = self.slaves[slave_tid]['nickname']
        if self.verbose:
            print '%d (%s) %d items left'%(slave_tid, nickname, n_left)

        chunk = self.get_slave_chunk( queue )

        self.send(slave_tid, MSG_JOB_START, (chunk,))


    def is_valid_slave(self, slave_tid):
        """
        Checked each time before a new job is given to a slave, if 0, the
        job is given to another slave. Override.

        @param slave_tid: slave task tid
        @type  slave_tid: int        
        """
        return 1


    def job_done(self, slave_tid, result):
        """
        Override to add tasks to be preformend when the job is done
        (other than the default, see L{__job_done}).
        
        @param slave_tid: slave task tid
        @type  slave_tid: int
        """
        pass


    def __job_done(self, slave_tid, result):
        """
        Tasks that are preformed when the job is done.
        
        @param slave_tid: slave task tid
        @type  slave_tid: int
        @param result: slave result dictionary
        @type  result: dict
        """
        ## synchronize on internal lock of Status to avoid the distribution
        ## of new items while processed ones are not yet marked "finished"
        self.status.lock.acquire()

        self.job_done(slave_tid, result)

        self.result.update(result)

        ## mark result as finished.
        for item in result.keys():
            self.status.deactivate(item)

        self.status.lock.release()

        ## once again
        if not self.is_valid_slave(slave_tid):
            return

        if not self.status.done():
            self.__start_job(slave_tid)
        else:
            self.__finish()


class JobSlave(PVMMasterSlave):

    def __init__(self):
        """
        """
        PVMMasterSlave.__init__(self)
        self.setMessageLoopDelay(0.5)


    def start(self):
        """
        """
        PVMMasterSlave.start(self)

        self.bindMessages()
        self.startMessageLoop()


    def bindMessages(self):
        """
        """
        parent = self.getParent()
        self.bind(MSG_JOB_START, parent, self.__go)


    def initialize(self, params):
        """
        Automatically invoked by parent after slave's
        message-loop is up. Override to use.
        """
        pass


    def go(self, *args, **kw):
        """
        Must be overridden in order to do the actual work.
        Result should be returned.
        Default tasks are defined in L{__go}.

        @param args: arguments
        @type  args: (any)
        @param kw: dictionary with key=value pairs
        @type  kw: {key:value}
        """
        pass


    def __go(self, *args, **kw):
        """
        Startup tasks.
        
        @param args: arguments
        @type  args: (any)
        @param kw: dictionary with key=value pairs
        @type  kw: {key:value}        
        """
        result = self.go(*args, **kw)

        ## send result back to parent
        my_tid = self.getTID()

        self.send(self.getParent(), MSG_JOB_DONE, (my_tid, result))

################
## empty test ##
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, dispatcher is tested by other PVM tests."""
    pass


if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = JobSlave()
    slave.start()
