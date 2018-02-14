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
## Contributions: Raik Gruenberg

"""
Thread - save job handling for JobMaster (see L{dispatcher}).
"""

from threading import RLock

class Status:
    """
    Keep track of objects that are processed by JobMaster.
    
    Thread-savety:
     - next_chunk() is synchronized.
     - deactivate() must be synchronized to the same Status.lock
     - activate() is not any longer supposed to be called from outside.
    """

    def __init__(self, objects, redistribute=1):
        """
        @param objects: e.g. job IDs
        @type  objects: any
        @param redistribute: send out started but not yet finished jobs
                             to idle slaves (to not wait for slow slaves)
                             (default: 1)
        @type  redistribute: 1|0
        """
        self.objects = {}
        map(lambda k, o = self.objects: o.update({k: None}), objects)

        self.redistribute = redistribute

        self.lock = RLock()


    def __activate(self, item):
        """
        not thread-save, synchronize on self.lock!
        """
        if self.objects[item] is None:
            self.objects[item] = 1
        else:
            self.objects[item] += 1


    def deactivate(self, item):
        """
        mark item as processed.
        not thread-save, synchronize on self.lock!
        """
        if self.objects[item] is None:
            raise '%s has never been activated.' % item

        self.objects[item] = 0


    def not_done(self):
        """
        Not yet finished items. The not-yet-started come first, then
        come the active items ordered by how often they have already been
        distributed.
        
        @return: list of not yet finished items (started or not)
        @rtype: [object]
        """
        o = self.objects
        r = [ k for k in o.keys() if (o[k] != 0) ]

        ## sort by number of times that the item has already been distributed
        pairs = [(o[k], k) for k in r ]
        pairs.sort()
        r = [ x[1] for x in pairs ]

        return r


    def not_started(self):
        """
        @return: list of not yet started items
        @rtype: [object]
        """
        r = filter(lambda key, s = self: s.objects[key] is None,
                      self.objects.keys())

        return r


    def activ(self):
        """
        @return: list of currently active (started) items
        @rtype: [object]
        """
        o = self.objects
        r = [ k for k in o.keys() if (o[k] is not None) and (o[k] > 0) ]

        return r


    def done(self):
        """
        @return: 1 if all items have been processed and finished
        @rtype: 1|0
        """
        r = not len(filter(lambda v: v <> 0, self.objects.values()))

        return r


    def next_chunk(self, nmax ):
        """
        Get next chunk of at most nmax items that need to be processed.
        Thread-save.

        @param nmax: size of chunk
        @type  nmax: int        

        @return: chunk of items to be processed, number of unproc OR
                   None if all items have been processed
        @rtype: ([object], int) OR None
        """
        self.lock.acquire()

        queue = self.not_started()

        ## at the end, re-distribute running jobs
        if not queue and self.redistribute:
            queue = self.not_done()

        n_left = len(queue)

        if n_left > nmax:
            queue = queue[:nmax]

        for item in queue:
           self.__activate(item)

        self.lock.release()

        return queue, n_left


    def __str__(self):
        names = { None: 'not_started',
                  0: 'finished',
                  1: 'active'}

        l = []

        for key, value in self.objects.items():

            key = str(key)

            l.append(key + ': ' + names.get(value, 'active+%s'%str(value)))

        l.append('Status: ' + 'not ' * (not self.done()) + 'done')

        return '\n'.join(l)

    __repr__ = __str__

################
## empty test ##
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, Status is tested indirectly in other PVM applications."""
    pass
