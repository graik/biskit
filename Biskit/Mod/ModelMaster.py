##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2005, Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: Olivier PERIN, Raik Gruenberg
##
## last $Author$
## last $Date$
"""
Parallelize Modeller runs
"""

from Biskit.PVM.TrackingJobMaster import *

import Biskit.hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os, glob

class ModelMaster(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/Mod/ModelSlave.py'

    def __init__(self, hosts, folders, fastaTarget=None, f_pir=None, template_folder=None, fout=None, starting_model=1, ending_model=10, ferror=None, **kw):
        """
        folders         - [str], list of project directories (path)
        fastaTarget     - str, path to find 'target.fasta'
        f_pir           - str, alignment filename
        template_folder - str, directories for input atom files
        fout            - str, default modeller input file 'modeller.top'
        starting_model  - int, index of the first model
        ending_model    - int, index of the last model
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
        self.fastaTarget = fastaTarget
        self.f_pir = f_pir
        
        self.template_folder = template_folder
        self.fout = fout
        self.starting_model = starting_model
        self.ending_model = ending_model
        self.ferror = ferror or 'ModelSlaveErrors.out'
        
        data = self.setupJobs()
       

        TrackingJobMaster.__init__(self, data=data, chunk_size=1, hosts=hosts, slave_script = self.slave_script, redistribute=0, **kw)
     

    def __dir_or_none( self, folder, filename ):

        if filename is None:
            return None
        return os.path.join( folder, filename )

    def setupJobs(self):
        """
        Prepare the job dictionnary for 'ModelSlave'
        -> {{str}}, input informations for modeller for each project
        """

        r = {}

        for f in self.folders:
            modeller_input = {}
            modeller_input['outFolder'] = T.absfile(f)
            modeller_input['fastaTarget'] = self.__dir_or_none( f, self.fastaTarget )
            modeller_input['f_pir'] = self.__dir_or_none( f, self.f_pir )
            modeller_input['template_folder'] = self.__dir_or_none( f, self.template_folder )
            modeller_input['starting_model'] = self.starting_model
            modeller_input['ending_model'] = self.ending_model

            r[T.absfile(f)] = modeller_input

        return r


    def getInitParameters(self, slave_tid):
        """ hand over parameters to slave once. """
        
        return {'progress_str':'slave calculating..', 'ferror':self.ferror}

    def cleanup( self ):
        print "Cleaning up..."

    def done( self ):
        print "Now we are done."
    

if __name__ == '__main__':

    options = T.cmdDict({'d':glob.glob(T.absfile('~/Homstrad_modeller/*'))})

    d = options['d']
    
    master = ModelMaster(folders = d, hosts=hosts.cpus_all[ : 10 ])

    r = master.calculateResult()
