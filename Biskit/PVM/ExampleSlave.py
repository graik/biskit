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
An example of a Slave/Master setup.
"""

from Biskit.PVM import JobSlave
import Biskit.tools as T

import time

class Slave(JobSlave):

    def initialize(self, params):
        """
        Parameters handed over by master.
        """
        self.params = params
        ## or perhaps easier:
        ## self.__dict__.update( params )


    def go(self, dict):
        """
        Perform slave task.
        """
        d = {}

        T.flushPrint( self.params['progress_str'] )
        for id, val in dict.items():

            d[id] = val+1

            for i in range(1, 1000):
                f = 1.0 * i / 1200232112312.11

            T.flushPrint( str(id) + ' ' )
            time.sleep(0.5)

        print "Done."
        return d

################
## empty test ##
import Biskit.test as BT

from threading import Thread, RLock, Condition
import time

class Flusher( Thread ):

    def __init__( self, *f ):
        Thread.__init__( self )
        self.setDaemon( True )

        self.files = f
        self.stop = False

    def setStop( self ):
        self.stop = True

    def run( self ):

        while not self.stop:
            for f in self.files:
                f.flush()
            time.sleep( 2 )


class Test(BT.BiskitTest):
    """Mock test, Slave is tested by ExampleMaster."""
    pass


if __name__ == '__main__':

    import tempfile, os, sys

    f_out = open( tempfile.mktemp( '.out', 'slave_', os.getcwd() ), 'w' )

    sys.stdout = f_out
    sys.stderr = f_out
    flusher = Flusher( f_out )
    flusher.start()

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = Slave()
    slave.start()
