##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
Parallelize Modeller runs
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster

import Biskit.hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os, glob

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
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0, run=0, model_testRoot=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0        
        @param run: run the full test (call external application) or not
        @type  run: 1|0
        @param model_testRoot: align the full validation project in testRoot
        @type  model_testRoot: 1|0

        @return: 1
        @rtype:  int
        """
        import tempfile
        import shutil

        if model_testRoot:
            
            projRoot =  T.testRoot()+'/Mod/project'
            projects = glob.glob( projRoot + '/validation/*' )
            
            master = ModelMaster(folders = projects,
                                 hosts=hosts.cpus_all[ : 10 ])
            
            r = master.calculateResult()
            
            if local:
                globals().update( locals() )
            
            return 1


        ## collect the input files needed
        outfolder = tempfile.mkdtemp( '_test_Modeller' )
        os.mkdir( outfolder +'/templates' )
        os.mkdir( outfolder +'/t_coffee' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/modeller',
                         outfolder + '/templates/modeller' )

        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     outfolder + '/t_coffee' )    

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                         outfolder  )

        master = ModelMaster(folders = [outfolder],
                             hosts=hosts.cpus_all[ : 10 ],
                             show_output = local,
                             verbose = local )

        if run:
            r = master.calculateResult()
            if local: print 'The models result can be found in %s/modeller'%outfolder
            
        if local:
            globals().update( locals() )
            
        ## cleanup
        T.tryRemove( outfolder, tree=1 )
        
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
    
    assert test.run( run=1, local=1 ) ==  test.expected_result()

