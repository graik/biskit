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
##
##

"""
Parallellized AmberEntropist calculation.
"""

from Biskit.PVM import JobSlave
import Biskit.tools as T
import Biskit.settings as settings

from Biskit import LogFile
from Biskit.AmberEntropist import AmberEntropist
from Biskit.AmberCrdEntropist import EntropistError

import os, time

class AmberEntropySlave( JobSlave ):
    """
    Collect AmberEntropist jobs from AmberEntropyMaster and return result.
    """

    def initialize(self, params):
        """
        expects::
          {'nice':int, 'ferror':str, .. }

        @param params: initialisation parameters passed from the master
        @type  params: dict
        """
        self.__dict__.update( params )
        self.errorLog = LogFile( self.ferror, mode='a' )


    def reportError(self, msg, id ):
        try:
            try:
                print msg
            except:
                pass

            msg = 'trouble with ' + msg
            s = '%s on %s, run %s\n' % (msg, os.uname()[1], id)
            s += '\Error:' + T.lastError()
            s += '\nErrorTrace:\n' + T.lastErrorTrace() + '\n'
            s += '\n'

            self.errorLog.add( s )

        except Exception, why:
            f = open('ErrorReportError_AmberEntropySlave','a')
            f.write( str(type(why)) )
            try:
                f.write( T.lastErrorTrace() )
            except:
                pass
            f.close()


    def go(self, jobs):
        """
        The calculation.
        
        @param jobs: dictionary with { int_id : str_protocol }
        @type  jobs: dict

        @return: result from AmberEntropist.run()
        @rtype: dict
        """
        result = {}

        startTime = time.time()

        for id, protocol in jobs.items():

            try:
                T.flushPrint( "%s " % str(id) )

                protocol.update( {'nice':self.nice} )

                x = None  ## free memory from previous run

                x = AmberEntropist( **protocol )

                x.run()

                r = x.result

                if r:
                    r['__version_AmberEntropist'] = x.version()

                    result[ id ] = r
                else:
                    result[ id ] = None

            except EntropistError, why:
                self.reportError( str(type(why)), id )
            except IOError, why:
                self.reportError( str(why), id )
            except Exception, why:
                self.reportError( 'ERROR '+str(type(why)), id )

        print "\navg time for last %i jobs: %f s" %\
              ( len(jobs), (time.time()-startTime)/len(jobs))

        return result


if __name__ == '__main__':

    import sys

    if len(sys.argv) == 2:

        nice = int(sys.argv[1])
        os.nice(nice)

    slave = AmberEntropySlave()
    slave.start()


