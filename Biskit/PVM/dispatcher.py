##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2002-2005 Wolfgang Rieping; All rights reserved
##
## Contributions: Wolfgang Rieping, Raik Gruenberg
## Created:  May 25, 2002
##

from PVMThread import PVMMasterSlave
import Biskit.settings as settings

MSG_JOB_START = 1
MSG_JOB_DONE = 2

class JobMaster(PVMMasterSlave):

    def __init__(self, data, chunk_size, hosts, niceness, slave_script,
                 show_output = 0, result = None, redistribute=1 ):
        """
        data: dict of items to be proessed {id:object}.
        chunk_size: number of items that are processed by a job
        hosts: list of host-names
        niceness: dictionary [host-name: niceness]
        slave_script: absolute path to slave-script
        result: dict; items which have already been processed
        (ie. they are contained in result) are not processed again.
        redistribute - 1||0, at the end, send same job out several times [1]
        """

        from Status import Status

        PVMMasterSlave.__init__(self)

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

        self.hosts = unique_list
        self.niceness = niceness
        self.data = data
        self.slave_script = slave_script
        self.chunk_size = chunk_size
        
        self.current_pos = 0
        self.show_output = show_output
        

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

        print 'Processing %d items ...' % len(items)
        
    def start(self):

        self.finished = 0

        PVMMasterSlave.start(self)
        
        self.startMessageLoop()
        self.spawnAll(self.niceness, self.show_output)

    def getInitParameters(self, slave_tid):
        return None


    def done(self):
        """
        Override to do something after last job result has been received.
        """
        print "Done"


    def finish( self ):
        """
        This method is called one time, after the master has received the
        last missing result.
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

        self.bind(MSG_JOB_DONE, slave_tid, self.__job_done)

    def spawn(self, host, nickname, niceness, show_output = 0):

        import pvm, socket

        if show_output:

            display = socket.gethostname() + ':0.0'

            command = settings.xterm_bin
            argv = ['-title', str(host), '-geometry', '30x10',
                    '-display', display, '-e',
                    settings.python_bin, '-i', self.slave_script,
                    str(niceness)]
        else:

            command = settings.python_bin
            argv = ['-i', self.slave_script, str(niceness)]

        args = (command, argv, pvm.spawnOpts['TaskHost'], host, 1)

        return PVMMasterSlave.spawn(self, args, nickname)

    def spawnAll(self, niceness, show_output = 0):

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
                print 'return code:', slave_tid

            else:

                self.bindMessages(slave_tid)
                
                self.slaves[slave_tid] = d

                print slave_tid, nickname, 'spawned.'

    def initializationDone(self, slave_tid):
        """
        is called by a slave that has been initialized and
        is now ready for start-up.
        """

        ## start processing

        self.__start_job(slave_tid)


    def start_job(self, slave_tid):
        """
        called when a new job is about to be started
        """
        pass


    def get_slave_chunk(self, data_keys ):
        """
        Assemble task dictionary that is send to the slave for a single job.
        Override this, if the values of self.data are to be changed/created on
        the fly.
        data_keys - [any], subset of keys to self.data
        -> { any:any }, dict mapping the data keys to the actual data values
        """
        chunk = {}

        for id in data_keys:
            chunk[id] = self.data[id]

        return chunk


    def __start_job(self, slave_tid):

        self.start_job(slave_tid)

        ## get items that have not been processed

        queue, n_left = self.status.next_chunk( self.chunk_size)

        if not queue:
            return

        nickname = self.slaves[slave_tid]['nickname']
        print '%d (%s) %d items left' % (slave_tid, nickname, n_left)
        
        chunk = self.get_slave_chunk( queue )
        
        self.send(slave_tid, MSG_JOB_START, (chunk,))


    def is_valid_slave(self, slave_tid):
        """
        Checked each time before a new job is given to a slave, if 0, the
        job is given to another slave. Override.
        """
        return 1

    def job_done(self, slave_tid, result):
        pass

    def __job_done(self, slave_tid, result):

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

        PVMMasterSlave.__init__(self)
        self.setMessageLoopDelay(0.5)

    def start(self):

        PVMMasterSlave.start(self)

        self.bindMessages()
        self.startMessageLoop()

    def bindMessages(self):

        parent = self.getParent()

        self.bind(MSG_JOB_START, parent, self.__go)

    def initialize(self, params):
        """
        automatically invoked by parent after slave's
        message-loop is up.
        """

        pass

    def go(self, *args, **kw):
        """
        must be overridden in order to do
        the actual work. result should be returned.
        """

        pass

    def __go(self, *args, **kw):

        result = self.go(*args, **kw)

        ## send result back to parent

        my_tid = self.getTID()

        self.send(self.getParent(), MSG_JOB_DONE, (my_tid, result))

if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:
        
        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = JobSlave()
    slave.start()
