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
                 ending_model=10, ferror=None, **kw):
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
        @param kw: additional TrackingJobMaster arguments::
            chunk_size   - int, number of items that are processed per job
            niceness     - {str_host-name: int_niceness}
            slave_script - str, absolute path to slave-script
            show_output  - 1|0, display one xterm per slave [0]
            add_hosts    - 1|0, add hosts to PVM before starting [1]
        @type kw: param=value
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
        print "Cleaning up..."


    def done( self ):
        print "Done modeling."


if __name__ == '__main__':

    options = T.cmdDict({'d':glob.glob(T.absfile('~/Homstrad_modeller/*'))})

    d = options['d']

    master = ModelMaster(folders = d, hosts=hosts.cpus_all[ : 10 ])

    r = master.calculateResult()
