##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## $Revision$
## last $Date$
## last $Author$

"""
An example of a Master/Slave setup
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster

from Biskit.hosts import nodes_all, cpus_all
from Biskit.tools import projectRoot


class Master(TrackingJobMaster):

    ## Slave script that ges with this master
    slave_script =  projectRoot() + '/Biskit/PVM/ExampleSlave.py'

    
    def __init__(self, verbose=1, *args, **kw):
        """
        Parameters nedded by master and/or slave.
        """
        self.verbose = verbose

        TrackingJobMaster.__init__(self, verbose=verbose, *args, **kw)

                         
    def getInitParameters(self, slave_tid):
        """
        Hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}   
        """
        return {'progress_str':'slave calculating..'}


    def cleanup( self ):
        """
        Tidy up tasks.
        """
        if self.verbose: print "Cleaning up..."


    def done( self ):
        """
        Called when master is done.
        """
        if self.verbose: print "Now we are done."



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype:  int
        """
        import time
        from Biskit.PVM.ExampleMaster import Master

        hosts = cpus_all[:8]

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
        
        master = Master( data=data, chunk_size=10,
                         hosts=hosts, niceness=niceness,
                         slave_script=Master.slave_script,
                         show_output=local, redistribute=1,
                         verbose=local )

        ## blocking call
        r = master.calculateResult()

        ## Example for non-blocking call with saving of restart info
        ##     master.start()
        ##     time.sleep( 10 )
        ##     rst = master.getRst()
        
        if local:
            globals().update( locals() )
            
        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    

if __name__ == '__main__':

    test = Test()
    
    assert test.run( local=1 ) ==  test.expected_result()
