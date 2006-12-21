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
Parallizes the convertion of PDB files into pickled PDBModel objects.
"""

import Biskit.hosts as hosts
import Biskit.tools as T
## from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.PVM import TrackingJobMaster


class StructMaster(TrackingJobMaster):

    def __init__(self, dat, chunk, hosts, outFolder, skipWat=0, amber=0,
                 sort=0, add_hosts=1, **kw):
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
        import tempfile
        out_folder = tempfile.mkdtemp( '_test_StructureMaster' )
        
        pdbs = {T.testRoot() + '/lig/1A19.pdb':'',
                T.testRoot() + '/rec/1A2P.pdb':'',
                T.testRoot() + '/com/1BGS.pdb':''}

        master = StructMaster( pdbs,
                               2,
                               hosts=hosts.cpus_all,
                               outFolder=out_folder,
                               show_output=local,
                               verbose=local,
                               add_hosts=1 )

        ## run and wait for result
        r = master.calculateResult()

        if local:
            print 'The converted pdb files has been written to %s'%out_folder
            globals().update( locals() )

        ## cleanup
        T.tryRemove( out_folder, tree=1 )

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

