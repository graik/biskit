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

## Contributions: Olivier PERIN, Raik Gruenberg
"""
Parallelize Modeller runs
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster

import Biskit.PVM.hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os

class ModelMaster(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/Mod/ModelSlave.py'

    def __init__(self, hosts, folders, fastaTarget=None, f_pir=None,
                 template_folder=None, fout=None, starting_model=1,
                 ending_model=10, ferror=None, verbose=1, **kw):
        """
        @param hosts: list of host-names
        @type  hosts: [str]
        @param folders: list of project directories (path)
        @type  folders: [str]
        @param fastaTarget: path to find 'target.fasta'
        @type  fastaTarget: str
        @param f_pir: alignment filename
        @type  f_pir: str
        @param template_folder: directories for input atom files
        @type  template_folder: str
        @param fout: default modeller input file 'modeller.top'
        @type  fout: str
        @param starting_model: index of the first model (default: 1)
        @type  starting_model: int
        @param ending_model: index of the last model (default: 10)
        @type  ending_model: int
        @param ferror: filename to output errors from the Slave
        @type  ferror: str
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        @param kw: additional TrackingJobMaster arguments::
            chunk_size   - int, number of items that are processed per job
            niceness     - {str_host-name: int_niceness}
            slave_script - str, absolute path to slave-script
            show_output  - 1|0, display one xterm per slave [0]
            add_hosts    - 1|0, add hosts to PVM before starting [1]
        @type kw: param=value
        """
        self.verbose = verbose
        self.folders = folders
        self.fastaTarget = fastaTarget
        self.f_pir = f_pir

        self.template_folder = template_folder
        self.fout = fout
        self.starting_model = starting_model
        self.ending_model = ending_model
        self.ferror = ferror or 'ModelSlaveErrors.out'

        data = self.setupJobs()


        TrackingJobMaster.__init__(self, data=data, chunk_size=1, hosts=hosts,
                                   slave_script = self.slave_script,
                                   redistribute=0, verbose=verbose, **kw)


    def __dir_or_none( self, folder, filename ):
        if filename is None:
            return None
        return os.path.join( folder, filename )


    def setupJobs(self):
        """
        Prepare the job dictionnary for 'ModelSlave'

        @return: input informations for modeller for each project
        @rtype: {{str}}
        """
        r = {}

        for f in self.folders:
            modeller_input = {}
            modeller_input['outFolder'] = T.absfile(f)
            modeller_input['fastaTarget'] = \
                          self.__dir_or_none( f, self.fastaTarget )
            modeller_input['f_pir'] = self.__dir_or_none( f, self.f_pir )
            modeller_input['template_folder'] = \
                          self.__dir_or_none( f, self.template_folder )
            modeller_input['starting_model'] = self.starting_model
            modeller_input['ending_model'] = self.ending_model

            r[T.absfile(f)] = modeller_input

        return r


    def getInitParameters(self, slave_tid):
        """
        Hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionary with init parameters
        @rtype: {param:value}
        """
        return {'progress_str':'slave calculating..', 'ferror':self.ferror}


    def cleanup( self ):
        if self.verbose: print "Cleaning up..."


    def done( self ):
        if self.verbose: print "Done modeling."




#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """

    TAGS = [BT.PVM, BT.LONG, BT.EXE]

    def prepare(self):
        import tempfile
        import shutil

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_Modeller' )
        os.mkdir( self.outfolder +'/templates' )
        os.mkdir( self.outfolder +'/t_coffee' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/modeller',
                         self.outfolder + '/templates/modeller' )

        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     self.outfolder + '/t_coffee' )    

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     self.outfolder  )


    def test_ModelMaster( self):
        """Mod.ModelMaster test"""

        nodes = hosts.cpus_all[ : 10 ]

        self.master = ModelMaster(folders = [self.outfolder],
                                  hosts=nodes,
                                  show_output = self.local,
                                  verbose = self.local )

        assert len(nodes) > 0, 'master needs at least 1 pvm node.'

        self.r = self.master.calculateResult()

        if self.debug and self.local:
            print 'The models result can be found in %s/modeller'%\
                  self.outfolder

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )


if __name__ == '__main__':

    BT.localTest()
