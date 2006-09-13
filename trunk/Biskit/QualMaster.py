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
##
## last $Author$
## last $Date$
## $Revision$

"""
Parallizes calculation and plotting of ensemble trajectory RMSDs
"""

from Biskit.hosts import cpus_all, nice_dic, nodes_all
import Biskit.tools as T
from Biskit.PVM.TrackingJobMaster import TrackingJobMaster


class QualMaster(TrackingJobMaster):

    def __init__(self, trajFiles, n_hosts=20, **kw):
        """
        @param dat: data dictionary
        @type  dat: dict
        @param hosts: list of host-names
        @type  hosts: [str]
        """
        dat = {}
        i = 0
        for f in trajFiles:
            dat[i] = T.absfile( f )
            i += 1

        niceness = nice_dic
        hosts = nodes_all[ :n_hosts ]

        slave_script = T.projectRoot() +'/Biskit/QualSlave.py'

        TrackingJobMaster.__init__(self, dat, 1, hosts, niceness,
                                   slave_script, **kw)


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}
        """
        return 1


    def done(self):
        self.exit()


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
        @rtype: int
        """
        import os.path
        
        ## a minimal list of trajectories
        traj_list = [ T.testRoot() + '/lig_pcr_00/traj.dat' ]

        master = QualMaster( traj_list,
                             show_output=local,
                             verbose=local )

        ## run and wait for result
        r = master.calculateResult()
        #master.start()

        f_plot = '%s/rms_traj.eps'%os.path.dirname(r[0])
        
        if local:
            print 'A RMSD plot is writen to: %s'%f_plot
            globals().update( locals() )
            
        ## cleanup
        T.tryRemove( f_plot )
        
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

    assert test.run( local=1 ) == test.expected_result()

