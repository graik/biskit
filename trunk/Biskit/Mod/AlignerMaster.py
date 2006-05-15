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

## Contributions: Olivier PERIN, Raik Gruenberg
## last $Author$
## last $Date$
## $Revision$
"""
Parallelise Sequence Alignment
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster

import Biskit.hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os, glob

class AlignerMaster(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/Mod/AlignerSlave.py'

    def __init__(self, hosts, folders, pdbFolder=None, fastaTemplates=None,
                 fastaSequences=None, fastaTarget=None, ferror=None, **kw):
        """
        @param hosts: list of host-names
        @type  hosts: [str]
        @param folders: list of project directories (path)
        @type  folders: [str]
        @param pdbFolder: pdb directory (*.alpha for Aligner)
        @type  pdbFolder: str
        @param fastaTemplates: path to find 'templates.fasta'
        @type  fastaTemplates: str
        @param fastaSequences: path to find 'nr.fasta'
        @type  fastaSequences: str
        @param fastaTarget: path to find 'target.fasta'
        @type  fastaTarget: str
        @param ferror: filename to output errors from the Slave
        @type  ferror: str
        @param kw: additional TrackingJobMaster arguments::
            chunk_size   - int, number of items that are processed per job
            niceness     - {str_host-name: int_niceness}
            slave_script - str, absolute path to slave-script
            show_output  - 1|0, display one xterm per slave [0]
            add_hosts    - 1|0, add hosts to PVM before starting [1]
        @type kw: param=value
        """

        self.folders = folders
        self.pdbFolder = pdbFolder 
        self.fastaTemplates = fastaTemplates 
        self.fastaSequences = fastaSequences 
        self.fastaTarget = fastaTarget
        self.ferror = ferror or 'AlignSlaveErrors.out'

        data = self.setupJobs()

        TrackingJobMaster.__init__(self, data=data, chunk_size=1,
                                   hosts=hosts,
                                   slave_script=self.slave_script,
                                   redistribute=0, **kw)


    def __dir_or_none( self, folder, filename ):
        if filename is None:
            return None
        return os.path.join( folder, filename )


    def setupJobs(self):
        """
        Prepare the job dictionnary for 'AlignerSlave'
        
        @return: input informations for aligner for each project
        @rtype: {{str}}
        """
        r = {}

        for f in self.folders:
            aligner_input = {}
            aligner_input['outFolder'] = T.absfile(f)
            aligner_input['fastaTemplates'] = \
                                  self.__dir_or_none( f, self.fastaTemplates )
            aligner_input['fastaSequences'] = \
                                  self.__dir_or_none( f, self.fastaSequences )
            aligner_input['fastaTarget'] = \
                                  self.__dir_or_none( f, self.fastaTarget )

            pdb_list = []
            if self.pdbFolder is None:
                pdb_list = None
            else:
                pdbfiles = os.listdir(f + self.pdbFolder)
                for pdb in pdbfiles:
                    pdb_list.append(f + self.pdbFolder + '/%s'%pdb)

            aligner_input['pdbFiles'] = pdb_list

            r[T.absfile(f)] = aligner_input

        return r


    def getInitParameters(self, slave_tid):
        """
        Hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}        
        """
        return {'progress_str':'slave calculating..',
                'ferror':self.ferror, 'os.environ':os.environ }


    def cleanup( self ):
        print "Cleaning up..."


    def done( self ):
        print "Done aligning."


if __name__ == '__main__':

    mask     = T.absfile(T.testRoot()+'/Mod/project/validation/1MQ1')
    projects = glob.glob( mask )

    master = AlignerMaster(folders=projects,
                      ferror=T.testRoot()+'/Mod/project/AlignErrors.out',
                      hosts=hosts.cpus_all[ : 10 ],
                      show_output=1)

    r = master.calculateResult()






