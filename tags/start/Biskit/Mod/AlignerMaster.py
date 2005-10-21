##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004, Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: Olivier PERIN, Raik Gruenberg
##
## last $Author$
## last $Date$

from Biskit.PVM.TrackingJobMaster import *

import Biskit.hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os, glob

class AlignerMaster(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/Mod/AlignerSlave.py'

    def __init__(self, hosts, folders, pdbFolder=None, fastaTemplates=None,
                 fastaSequences=None, fastaTarget=None,ferror=None, **kw):
        """
        folders         - [str], list of project directories (path)
        pdbFolder       - str, pdb directory (*.alpha for Aligner)
        fastaTemplates  - str, path to find 'templates.fasta'
        fastaSequences  - str, path to find 'nr.fasta'
        fastaTarget     - str, path to find 'target.fasta'
        ferror          - str, filename to output errors from the Slave
        
        ### TrackingJobMaster arguments ###
        chunk_size   - int, number of items that are processed per job
        hosts        - [ str ], list of host-names
        niceness     - {str_host-name: int_niceness}
        slave_script - str, absolute path to slave-script
        show_output  - 1||0, display one xterm per slave [0]
        add_hosts    - 1||0, add hosts to PVM before starting [1]
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
        -> {{str}}, input informations for aligner for each project
        """

        r = {}

        for f in self.folders:
            aligner_input = {}
            aligner_input['outFolder'] = T.absfile(f)
            aligner_input['fastaTemplates'] = self.__dir_or_none( f, self.fastaTemplates )
            aligner_input['fastaSequences'] = self.__dir_or_none( f, self.fastaSequences )
            aligner_input['fastaTarget'] = self.__dir_or_none( f, self.fastaTarget )

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
        """ hand over parameters to slave once. """
        
        return {'progress_str':'slave calculating..',
                'ferror':self.ferror, 'os.environ':os.environ }

    def cleanup( self ):
        print "Cleaning up..."

    def done( self ):
        print "Now we are done."
    

if __name__ == '__main__':

    mask     = T.absfile(T.testRoot()+'/Mod/project/validation/1MQ1')
    projects = glob.glob( mask )

    master = AlignerMaster(folders=projects,
                      ferror=T.testRoot()+'/Mod/project/AlignErrors.out',
                      hosts=hosts.cpus_all[ : 10 ],
                      show_output=1)

    r = master.calculateResult()

   
    


   
