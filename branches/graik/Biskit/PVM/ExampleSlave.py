from Biskit.PVM.dispatcher import JobSlave
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
import Biskit.tools as T
import subprocess

import time

class Slave(JobSlave):

    def initialize(self, params):

        self.params = params
        ## or perhaps easier:
        ## self.__dict__.update( params )

    def go(self, dict):

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
    
if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:
        
        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = Slave()
    slave.start()
