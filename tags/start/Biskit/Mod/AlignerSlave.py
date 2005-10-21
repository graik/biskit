##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2005, Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: Olivier PERIN, Raik Gruenberg
##
## last $Author$
## last $Date$
"""
Parallelise Sequence Alignment
"""

from Biskit.PVM.dispatcher import JobSlave
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

        self.__dict__.update( params )
        self.params = params

        ## Only the PATH must be updated from the master to run properly
        os.environ["PATH"]=self.params['os.environ']["PATH"]
        
        self.errorLog = LogFile( self.ferror, mode='a' )

     

    def reportError(self, msg, d ):
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


    def prepareT_coffee(self, input_file):

        alpha_index = open(T.absfile('%s'%input_file,'a+'))

        string_lines = alpha_index.readlines()

        alpha_path = []

        for line in string_lines:

            alpha_path.append(line[:-1])

        return alpha_path
    

    def go(self, dict):

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
    
if __name__ == '__main__':

    slave = AlignerSlave()
    slave.start()
