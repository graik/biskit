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
Parallizes the convertion of PDB files into pickled PDBModel objects.
"""

import Biskit.PVM.hosts as hosts
import Biskit.tools as T
## from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.PVM import TrackingJobMaster


class StructMaster(TrackingJobMaster):

    def __init__(self, dat, chunk, hosts, outFolder, skipWat=0, amber=0,
                 sort=0, **kw):
        """
        @param dat: data dictionary
        @type  dat: dict
        @param chunk: chunk size passed to slave
        @type  chunk: int
        @param hosts: list of host-names
        @type  hosts: [str]
        @param outFolder: alternative output folder
        @type  outFolder: str
        """
        niceness = {'default': 0}
        slave_script = T.projectRoot() + '/Biskit/StructureSlave.py'

        TrackingJobMaster.__init__(self, dat, chunk, hosts, niceness,
                                   slave_script, **kw)

        self.options = {}
        self.options['out'] = outFolder

        self.options['skipRes'] = None
        if skipWat:
            self.options['skipRes'] = ['WAT','TIP3','H2O','WWW','Na+','Cl-']

        if kw.has_key('show_output'):
            self.options['report'] = not kw['show_output']

        self.options['amber'] = amber

        self.options['sort'] = sort


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}
        """
        return self.options


    def done(self):

        self.exit()



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.PVM ]

    def prepare(self):
        import tempfile
        self.out_folder = tempfile.mkdtemp( '_test_StructureMaster' )

    def test_StructureMaster(self):
        """StructureMaster test"""
        import os

        pdbs = {T.testRoot() + '/lig/1A19.pdb':'',
                T.testRoot() + '/rec/1A2P.pdb':'',
                T.testRoot() + '/com/1BGS.pdb':''}

        self.master = StructMaster( pdbs,
                                    2,
                                    hosts=hosts.cpus_all,
                                    outFolder= self.out_folder,
                                    show_output=self.local,
                                    verbose=self.local,
                                    add_hosts=1 )

        ## run and wait for result
        self.r = self.master.calculateResult()

        if self.local:
            print 'The converted pdb files has been written to %s' \
                  % self.out_folder

        self.assert_( os.path.exists( self.out_folder + '/1A19.model') )

    def cleanUp(self):
        T.tryRemove( self.out_folder, tree=1 )

if __name__ == '__main__':

    BT.localTest()
