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

from Biskit.PVM import JobSlave
import Biskit.tools as T

from Biskit.Mod.Aligner import Aligner
from Biskit.Mod.TemplateCleaner import TemplateCleaner as TC
from Biskit.Mod.ValidationSetup import ValidationSetup as VS

from Biskit import LogFile
import os

class AlignerSlave(JobSlave):
    """
    See also: Aligner.py, AlignerMaster.py
    """

    def initialize(self, params):
        """
        Initialize AlignerSlave.

        @param params: dictionary with init parameters
        @type  params: {param:value}          
        """
        self.__dict__.update( params )
        self.params = params

        ## Only the PATH must be updated from the master to run properly
        os.environ["PATH"]=self.params['os.environ']["PATH"]

        self.errorLog = LogFile( self.ferror, mode='a' )



    def reportError(self, msg, d ):
        """
        Report error.

        @param msg: error message
        @type  msg: str
        @param d: error data
        @type  d: any
        """
        try:
            s = '%s on %s, job %r\n' % (msg, os.uname()[1], d)
            s += '\nErrorTrace:\n' + T.lastErrorTrace() + '\n'

            self.errorLog.add( s )

            try:
                print msg
            except:
                pass
        except Exception, why:
            f = open('ErrorReportError_XRefineSlave','a')
            f.write( str(why) )
            try:
                f.write( T.lastErrorTrace() )
            except:
                pass
            f.close()


    def prepareT_coffee(self, input_file):
        """
        Prepare list of coordinate files (.alpha)

        @param input_file: file with list of .alpha files
        @type  input_file: str

        @return: list of file names
        @rtype: [str]
        """
        alpha_index = open(T.absfile('%s'%input_file,'a+'))

        string_lines = alpha_index.readlines()

        alpha_path = []

        for line in string_lines:

            alpha_path.append(line[:-1])

        return alpha_path


    def go(self, dict):
        """
        Run alignment job.

        @param dict: dictionary with run parameters
        @type  dict: {param:value} 
        """
        d = {}
        val = None

        try:

            T.flushPrint( self.progress_str )
            for id, val in dict.items():

                aligner_log = LogFile( '%s/Aligner.log' %val["outFolder"] )

                d[id] = val

                aligner_log.add('Slave aligns %s on %s' % (id,os.uname()[1]) )

                a = Aligner( outFolder= val["outFolder"], log=aligner_log)

                ## For the cross validation
                if not os.path.exists(val["outFolder"] + TC.F_COFFEE):

                    input_file = val["outFolder"] + VS.F_TCOFFEE

                    alpha_path = self.prepareT_coffee(input_file)

                    a.align_for_modeller_inp( pdbFiles=alpha_path,
                              fasta_templates=val["fastaTemplates"],
                              fasta_sequences=val["fastaSequences"],
                              fasta_target=val["fastaTarget"])

                ## For a classic project folder    
                else:
                    a.align_for_modeller_inp(pdbFiles=val["pdbFiles"],
                              fasta_templates=val["fastaTemplates"],
                              fasta_sequences=val["fastaSequences"],
                              fasta_target=val["fastaTarget"])

                a.go()


        except Exception, why:
            self.reportError( 'ERROR '+str(why), val )

        print "Done."

        return d

##############
## empty test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, the Slave is tested in L{Biskit.Mod.AlignerMaster}"""
    pass

if __name__ == '__main__':

    slave = AlignerSlave()
    slave.start()
