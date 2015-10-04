## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

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

## last $Author$
## last $Date$
## $Revision$

"""
Binding incoming pvm-messages to methods.
"""

import pypvm as P
from Biskit.PVM import pvm
from threading import Thread

MSG_PING = 999999
MSG_MESSAGELOOP_UP = 999998
MSG_INITIALIZE = 999997
MSG_INITIALIZATION_DONE = 999996
MSG_EXIT = 999995

class Logfile:

    def __init__(self, filename = None):

        if filename is None:

            import os

            filename = '/tmp/pvm_thread.log.' + str(os.getpid())

        self.filename = filename
        self.lines = []


    def post(self, message):

        import time

        self.lines.append('\n' + str(time.time()) + ':')
        self.lines.append(message)

        file = open(self.filename, 'w')
        file.write('\n'.join(self.lines))
        file.close()


    def write(self):
        print '\n'.join(self.lines)


    def rm(self):
        import os
        if os.path.exists(self.filename):
            os.unlink(self.filename)


## This is the old implementation of the PVMThread class

class PVMThread(Thread):
    """
    a simple class binding incoming pvm-messages to methods.
    """


    def __init__(self, log = None):
        from threading import Event

        Thread.__init__(self)

        self.__mytid = P.mytid()

        try:
            ## get process ID of parent
            self.__parent = P.parent()
        except:
            self.__parent = None

        self.__bindings = {}
        self.__stop = 0
        self.__tasks = {}

        self.setMessageLoopDelay(0.1)
        self.__loopEvent = Event()

        if log: self.__log = Logfile()
        else: self.__log = None

        ## add ping-mechanism

        self.setPingTimeout(5.)
        self.bind(MSG_PING, -1, self._ping)


    def getTID(self):
        return self.__mytid


    def getParent(self):
        """
        @return: process ID of parent
        @rtype: int
        """
        return self.__parent


    def getBindings(self):
        return self.__bindings


    def bind(self, message_tag, tid, method):
        self.getBindings()[(tid, message_tag)] = method


    def unbind(self, tid_message_tag):
        del self.getBindings()[tid_message_tag]


    def setMessageLoopDelay(self, delay):
        """
        time is in seconds
        """
        self.__delay = delay


    def getMessageLoopDelay(self):
        return self.__delay


    def stopMessageLoop(self):
        self.__loopEvent.clear()


    def startMessageLoop(self):
        self.__loopEvent.set()


    def messageLoopIsStopped(self):
        return not self.__loopEvent.isSet()


    def setPingTimeout(self, t):
        self.__pingTimeout = t


    def getPingTimeout(self):
        return self.__pingTimeout


    def stop(self):
        self.__stop = 1
        ## otherwise this would pause the message-loop
        self.startMessageLoop()


    def isStopped(self):
        return self.__stop


    def getTasks(self):
        return self.__tasks


    def nicknameFromTID(self, tid):
        for nickname, TID in self.getTasks().items():
            if TID == tid:
                return nickname

        raise NameError, 'tid not known'


    def spawn(self, pvm_task, nickname = None):
        child_tid = P.spawn(*pvm_task)[0]

        if child_tid > 0:

            if nickname is None:
                nickname = child_tid

            self.__tasks[nickname] = child_tid

        return child_tid


    def send_primitive(self, tid, msg_tag, value):
        result = pvm.pack_and_send(tid, msg_tag, value)

        self.post_message_sent(msg_tag, tid, value)


    def send(self, task, msg_tag, value = None):
        """
        if 'task' is a tuple, msg_tag is send to all the tids
        given in that list.
        """
        from time import time

        if not type(task) == type(()):
            task = (task,)

        hash = self.getTasks()

        time_stamp = time()

        for t in task:
            ## if 'task' is a nickname, try to convert it into
            ## a valid tid

            try:
                tid = hash[t]
            except:
                tid = t

            self.send_primitive(tid, msg_tag, (time_stamp, value))


    def sendToAll(self, msg_tag, value):
        from time import time

        tids = tuple(self.getTasks().values())

        time_stamp = time()

        for tid in tids:
            self.send_primitive(tid, msg_tag, (time_stamp, value))

    def run(self):
        import time
        from numpy import argsort

        while not self.isStopped():

            bindings = self.getBindings()

            ## loop through all messages and check for
            ## incoming events

            incoming = {}

            for tid, message in bindings.keys():

                if P.probe(tid, message):

                    P.recv(tid, message)

                    ## parameters must be tuple

                    parameters = pvm.unpack()

##    self.post_message_received(message, tid, parameters)

                    value = (message, parameters[0], parameters[1])

                    try:
                        incoming[tid].append(value)
                    except:
                        incoming[tid] = [value]

            ## for every incoming message call
            ## bound method

            ## TODO: do this in a specific order!
            ## probably pvm supports a sort of time-stamp or so.

            for tid, values in incoming.items():

                time_stamps = map(lambda v: v[1], values)
                indices = argsort(time_stamps)

                for i in indices:

                    message = values[i][0]
                    parameters = values[i][2]

##    self.post_execute_method(message, tid, parameters)

                    if parameters is None:
                        bindings[(tid, message)]()
                    else:
                        bindings[(tid, message)](*parameters)

            ## wait some time

            time.sleep(self.getMessageLoopDelay())

            ## if message-loop is stopped, wait until
            ## it is continued

            self.__loopEvent.wait()


    def ping(self, nickname):
        self.send(nickname, MSG_PING, (self.getTID(),))

        tid = self.getTasks()[nickname]

        result = P.trecv(self.getPingTimeout(), MSG_PING, tid)

        try:
            pvm.unpack()
        except:
            pass

        if result <= 0:
            return 0
        else:
            return 1


    def _ping(self, sender_tid):
        self.send(sender_tid, MSG_PING, None)


    ## additional functionality for better debuging

    def log(self):
        if self.__log:
            self.__log.write()


    def post(self, message):
        if self.__log:
            self.__log.post(message)


    def post_message_received(self, msg_tag, tid, params):
        log_msg = 'PVMThread: received and unpacked: ' + \
                'msg_tag = %d, tid = %d\n' %(msg_tag, tid)
        log_msg += '[params = %s]' %str(params)
        self.post(log_msg)


    def post_message_sent(self, msg_tag, tid, params):
        log_msg = 'PVMThread: message sent: ' + \
                'msg_tag = %d, tid = %d\n' %(msg_tag, tid)
        log_msg += '[params = %s]' %str(params)
        self.post(log_msg)


    def post_execute_method(self, msg_tag, tid, params):
        log_msg = 'PVMThread: execute method bound to: ' + \
                'msg_tag = %d, tid = %d\n' %(msg_tag, tid)
        log_msg += '[params =  %s]' %str(params)
        self.post(log_msg)

    def rm_log(self):
        if self.__log:
            self.__log.rm()


class PVMMasterSlave(PVMThread):

    def __init__(self, verbose=1, *arg, **kw):
        PVMThread.__init__(self, *arg, **kw)

        self.verbose= verbose

        ## bind messages

        parent = self.getParent()

        if parent is not None:

            self.bind(MSG_INITIALIZE, parent, self.__initialize)
            self.bind(MSG_EXIT, parent, self.exit)


    def spawn(self, pvm_task, nickname = None):
        child_tid = PVMThread.spawn(self, pvm_task, nickname)

        ## bind messages being nessecary for automatic
        ## slave start-up determination

        self.bind(MSG_MESSAGELOOP_UP, child_tid, self.messageLoopIsUp)

        ## this message could also be bound later, before we send
        ## the MSG_INITIALIZE message.

        self.bind(MSG_INITIALIZATION_DONE, child_tid, self.initializationDone)

        return child_tid


    def startMessageLoop(self):
        PVMThread.startMessageLoop(self)

        ## send message to parent to indicate that
        ## we are now ready for message transfer.

        parent = self.getParent()

        if parent is not None:
            self.send(parent, MSG_MESSAGELOOP_UP, (self.getTID(),))


    def messageLoopIsUp(self, slave):
        """
        called by slave when its messge-loop has been started
        """
        ## initialize slave
        ## get init parameters for slave

        init_params = self.getInitParameters(slave)

        if init_params is not None:
            init_params = (init_params,)

        self.send(slave, MSG_INITIALIZE, init_params)


    def getInitParameters(self, slave_tid):
        return None


    def initialize(self, parameters):
        pass


    def __initialize(self, parameters):

        self.initialize(parameters)

        ## send message to parent, that are initialized properly
        parent = self.getParent()

        if parent is not None:
            self.send(parent, MSG_INITIALIZATION_DONE, (self.getTID(),))


    def initializationDone(self, slave_tid):
        pass


    def exit(self):
        parent = self.getParent()

        ## if we are the master, kill slaves (but don't send signal to self)
        if parent is None:

            for nickname, tid in self.getTasks().items():

                if tid != self.getTID():

                    self.send(nickname, MSG_EXIT, None)
                    if self.verbose: print nickname, 'shutting down...'
        else:
            P.kill(self.getTID())

        ## stop message-loop

        self.stop()

################
## empty test ##
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, PVMThread is tested indirectly by other PVM tests."""
    pass
