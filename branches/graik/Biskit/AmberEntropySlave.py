##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

from Biskit.PVM.dispatcher import JobSlave
import Biskit.tools as T
import Biskit.settings as settings

from Biskit import LogFile
from Biskit.AmberEntropist import AmberEntropist
from Biskit.AmberCrdEntropist import EntropistError

import os, time

class AmberEntropySlave( JobSlave ):

    def initialize(self, params):
        """
        expects
        {'nice':int, 'ferror':str, .. }
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
        """jobs - { int_soln : file_name_complex }"""

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

    import os, sys

    if len(sys.argv) == 2:
        
        nice = int(sys.argv[1])
        os.nice(nice)

    slave = AmberEntropySlave()
    slave.start()


