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
Parallelise Sequence Alignment
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster

from Biskit.PVM import hosts as hosts
from Biskit.tools import projectRoot
import Biskit.tools as T
import os, glob

class AlignerMaster(TrackingJobMaster):

    slave_script =  projectRoot() + '/Biskit/Mod/AlignerSlave.py'

    def __init__(self, hosts, folders, pdbFolder=None, fastaTemplates=None,
                 fastaSequences=None, fastaTarget=None, ferror=None,
                 verbose=1, **kw):
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
        self.pdbFolder = pdbFolder 
        self.fastaTemplates = fastaTemplates 
        self.fastaSequences = fastaSequences 
        self.fastaTarget = fastaTarget
        self.ferror = ferror or 'AlignSlaveErrors.out'

        data = self.setupJobs()

        TrackingJobMaster.__init__(self, data=data, chunk_size=1,
                                   hosts=hosts,
                                   slave_script=self.slave_script,
                                   redistribute=0,
                                   verbose=verbose, **kw)


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
import Biskit.test as BT

class FullProjectTest(BT.BiskitTest):
    """
    Performs full cross-validation alignments on the test project.
    Requires files that are not checked in but can be generated
    from the sequence in test/Mod/project.
    """

    TAGS = [ BT.LONG, BT.PVM, BT.EXE ]

    ## rename to test_fullProject to perform this test
    def t_fullProject( self ):
        """
        Mod.AlignerMaster full project test
        """
        ## a full test case that demands a valid testRoot project
        projRoot =  T.testRoot()+'/Mod/project'
        self.projects = glob.glob( projRoot + '/validation/*' )

        self.master = AlignerMaster(folders = self.projects,
                                    ferror = projRoot+'/AlignErrors.out',
                                    hosts = hosts.cpus_all[ : 10 ],
                                    show_output = self.local)

        self.r = self.master.calculateResult()


class TestBase( BT.BiskitTest ):

    def prepare(self):
        import tempfile
        import shutil

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_AlignerMaster' )
        os.mkdir( self.outfolder +'/templates' )
        os.mkdir( self.outfolder +'/sequences' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/t_coffee',
                         self.outfolder + '/templates/t_coffee' )

        shutil.copy( T.testRoot() + '/Mod/project/templates/templates.fasta',
                     self.outfolder + '/templates' )

        shutil.copy( T.testRoot() + '/Mod/project/sequences/nr.fasta',
                     self.outfolder + '/sequences/' )

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     self.outfolder  )


    def t_AlignerMaster(self, run=True):

        nodes = hosts.cpus_all[ : 5 ]

        self.master = AlignerMaster( folders=[self.outfolder],
                                     ferror=self.outfolder+'/AlignErrors.out',
                                     hosts=nodes,
                                     show_output=self.local,
                                     verbose=self.local )

        if run:
            assert len(nodes) > 0, 'master needs at least 1 pvm node.'

            self.r = self.master.calculateResult()
            if self.local and self.DEBUG:
                self.log.add('The alignment result is in %s/t_coffee'%\
                             self.outfolder)

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )


class TestDry( TestBase ):
    """Dry run test case that does not actually run t-coffee"""

    TAGS = [ BT.PVM ]

    def test_AlignerMaster(self):
        """Mod.AlignerMaster limited dry run test"""
        return self.t_AlignerMaster(run=False)

class TestReal( TestBase ):
    """Full AlignerMaster test case"""

    TAGS = [BT.PVM, BT.EXE, BT.LONG]

    def test_AlignerMaster(self):
        """Mod.AlignerMaster full test"""
        return self.t_AlignerMaster(run=True)


if __name__ == '__main__':

    BT.localTest()

