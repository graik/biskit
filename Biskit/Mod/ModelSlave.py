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

from Biskit.PVM.dispatcher import JobSlave
import Biskit.tools as T

from Biskit.Mod import Modeller as M

from Biskit import LogFile
import os
import subprocess
import settings
import socket


class ModelSlave(JobSlave):

    def initialize(self, params):
        """
        Initialize AlignerSlave.

        @param params: dictionary with init parameters
        @type  params: {param:value}          
        """
        self.__dict__.update( params )
        self.params = params

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
                f.write( t.lastErrorTrace() )
            except:
                pass
            f.close()


    def go(self, dict):
        """
        Run Modeller job.

        @param dict: dictionary with run parameters
        @type  dict: {param:value} 
        """
        d = {}
        val = None

        try:

            T.flushPrint( self.params['progress_str'] )
            for id, val in dict.items():

                modeller_log = LogFile( '%s/Modeller.log' %val["outFolder"] )     

                d[id] = val

                m = M( outFolder= val["outFolder"], log=modeller_log)

                m.prepare_modeller(fasta_target=val["fastaTarget"], f_pir=val["f_pir"], template_folder=val["template_folder"],starting_model=val["starting_model"], ending_model=val["ending_model"])

                m.go()

                m.postProcess()

        except Exception, why:
            self.reportError( 'ERROR '+str(why), val )

        print "Done."

        return d


if __name__ == '__main__':

    slave = ModelSlave()
    slave.start()
