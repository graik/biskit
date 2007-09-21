##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
from Biskit.PVM.TrackingJobMaster import *

from Biskit.hosts import nodes_all
from Biskit.tools import projectRoot

class Master(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/PVM/ExampleSlave.py'

    def __init__(self, *args, **kw):

        TrackingJobMaster.__init__(self, *args, **kw)

    def getInitParameters(self, slave_tid):
        """ hand over parameters to slave onnce. """
        return {'progress_str':'slave calculating..'}

    def cleanup( self ):
        print "Cleaning up..."

    def done( self ):
        print "Now we are done."
    

if __name__ == '__main__':

    import time
    from Biskit.PVM.ExampleMaster import Master as Master

    hosts = nodes_all[12:17]

    niceness = {'default': 20}

    data = {}
    for i in range(100):
        data[i]=i+1

    ## data - dictionary of some_id : some_object pairs
    ## 10 - how many items to pass to each slave in one go
    ## niceness - dict {'hostname':nice_value, .. or 'default':nice_value}
    ## slave_script - name of slave.py
    ## in the end get results from master.result
    ## -> dictionary of some_id : result_object pairs

    master = Master(data, 10, hosts, niceness, Master.slave_script,
                    show_output=1, redistribute=1 )

    ## blocking call
    r = master.calculateResult()

    ## Example for non-blocking call with saving of restart info
##     master.start()
##     time.sleep( 10 )
##     rst = master.getRst()