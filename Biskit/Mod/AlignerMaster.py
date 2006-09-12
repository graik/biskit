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



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0, run=0, align_testRoot=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0        
        @param run: run the full test (call external application) or not
        @type  run: 1|0
        @param align_testRoot: align the full validation project in testRoot
        @type  align_testRoot: 1|0

        @return: 1
        @rtype:  int
        """
        import tempfile
        import shutil

        ## a full test case that demands a valid testRoot project
        if align_testRoot:
            projRoot =  T.testRoot()+'/Mod/project'
            projects = glob.glob( projRoot + '/validation/*' )

            master = AlignerMaster(folders = projects,
                                   ferror = projRoot+'/AlignErrors.out',
                                   hosts = hosts.cpus_all[ : 10 ],
                                   show_output = 1)
            
            r = self.master.calculateResult()

            if local:
                globals().update( locals() )
            
            return 1

        ## collect the input files needed
        outfolder = tempfile.mkdtemp( '_test_AlignerMaster' )
        os.mkdir( outfolder +'/templates' )
        os.mkdir( outfolder +'/sequences' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/t_coffee',
                         outfolder + '/templates/t_coffee' )

        shutil.copy( T.testRoot() + '/Mod/project/templates/templates.fasta',
                     outfolder + '/templates' )

        shutil.copy( T.testRoot() + '/Mod/project/sequences/nr.fasta',
                     outfolder + '/sequences/' )

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     outfolder  )

        master = AlignerMaster(folders=[outfolder],
                               ferror=outfolder+'/AlignErrors.out',
                               hosts=hosts.cpus_all[ : 5 ],
                               show_output=1)

        if run:
            r = master.calculateResult()
            print 'The alignment result can be found in %s/t_coffee'%outfolder
            
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
    
    assert test.run( run=1, local=1 ) ==  test.expected_result()



