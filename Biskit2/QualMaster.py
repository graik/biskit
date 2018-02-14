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
##

"""
Parallizes calculation and plotting of ensemble trajectory RMSDs
"""

from Biskit.PVM.hosts import cpus_all, nice_dic, nodes_all
import Biskit.tools as T
from Biskit.PVM import TrackingJobMaster


class QualMaster(TrackingJobMaster):

    def __init__(self, trajFiles, n_hosts=20, **kw):
        """
        @param trajFiles: list of trajectory files
        @type  trajFiles: [ str ]
        @param n_hosts: number of hosts to use
        @type  n_hosts: int
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
        return {'verbose':self.verbose}


    def done(self):
        self.exit()


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.PVM ]

    def test_QualMaster(self):
        """QualMaster test"""
        import os.path

        ## a minimal list of trajectories
        traj_list = [ T.testRoot() + '/lig_pcr_00/traj.dat' ]

        self.master = QualMaster( traj_list,
                                  show_output=self.local,
                                  verbose=self.local )

        ## run and wait for result
        self.r = self.master.calculateResult()
        #master.start()

        self.f_plot = '%s/rms_traj.eps'%os.path.dirname( self.r[0] )

        if self.local:
            print 'A RMSD plot is writen to: %s'% self.f_plot

    def cleanUp(self):
        try:
            T.tryRemove( self.f_plot )
        except:
            pass


if __name__ == '__main__':

    BT.localTest()
