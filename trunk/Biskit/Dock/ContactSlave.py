## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## last $Date$
## last $Author$
## $Revision$

"""
Calculate contact matrix and some scores for complexes.
"""

from Biskit.PVM import JobSlave
import Biskit.tools as T
from Biskit import mathUtils as MU
from Biskit.LogFile import StdLog, LogFile
import numpy as N
from Complex import Complex
import os.path
import time


class ContactSlave(JobSlave):
    """
    Calculate contact matrix and some scores for complexes.
    ContactMaster creates several instances of this class, each on one node.
    """

    def version( self ):
        """
        Version of Dock.Complex

        @return: version of class
        @rtype: str
        """
        return 'ContactSlave $Revision$'


    def initialize(self, params):
        """
        Copy the parameters that ContactMaster is passing in as dict into
        fields of this class.

        @param params: defined in L{ContactMaster}
        @type  params: dict
        """
        self.ferror = params['ferror']
        self.log    = StdLog()
        if params['log']:        ## log to same file as master
            self.log = LogFile( self.flog )
        self.verbose = params['verbose']

        ## reference complex data
        self.c_ref_res_4_5 = self.c_ref_atom_4_5 = None
        self.c_ref_atom_10 = None
        self.mask_rec = self.mask_lig = None
        self.mask_interface = self.mask_interface_bb = None
        self.ref_interface_model = self.ref_interface_model_bb = None

        ## reduced models and reduced reference contacts
        self.reduced_recs = self.reduced_ligs = None
        self.c_ref_ratom_10 = None

        ## reference residue / atom contact matrices
        if params.get( 'c_ref_res_4_5', None):
            self.c_ref_res_4_5= MU.unpackBinaryMatrix(params['c_ref_res_4_5'])

            self.c_ref_atom_4_5 = MU.unpackBinaryMatrix(
                params['c_ref_atom_4_5'])
            self.c_ref_atom_10 = MU.unpackBinaryMatrix(
                params['c_ref_atom_10'])

        ## atom masks for casting
        if params.get('mask_rec', None):
            self.mask_rec = MU.unpackBinaryMatrix( params['mask_rec'] )
            self.mask_lig = MU.unpackBinaryMatrix( params['mask_lig'] )

        ## for interface rms calculation
        if params.get('mask_interface', None):
            self.mask_interface = MU.unpackBinaryMatrix(
                params['mask_interface'])
            self.mask_interface_bb = MU.unpackBinaryMatrix(
                params['mask_interface_bb'] )

            self.ref_interface_model =  params['ref_interface_model']
            self.ref_interface_model_bb =  params['ref_interface_model_bb']

        ## reduced rec and lig models indexed by source
        self.reduced_recs = params['reduced_recs']
        self.reduced_ligs = params['reduced_ligs']

        if params.get( 'c_ref_ratom_10', None):
            self.c_ref_ratom_10= MU.unpackBinaryMatrix(
                params['c_ref_ratom_10'])

        ## only calculate certain values
        self.force = params.get('force', [] )


    def reportError(self, msg, soln ):
        """
        Report any errors to log

        @param msg: error message
        @type  msg: str
        @param soln: solution number for complex giving the error
        @type  soln: int
        """
        try:
            s = '%s on %s, soln %i\n' % (msg, os.uname()[1], soln)
            s += '\t' + T.lastError() + '\n'
            s += 'TraceBack: \n' + T.lastErrorTrace() + '\n'
            f = open( self.ferror, 'a' )
            f.write( s )
            f.close()
        except:
            f = open('ErrorReportError_ContactSlave','a')
            f.write('')
            f.close()


    def __containsAny( self, lst_or_dic, *items ):
        """
        Check if dictionary or list contains items.

        @param lst_or_dic: lokk for items in this
        @type  lst_or_dic: list OR dict
        @param items: items to look for
        @type  items: any

        @return: result of test
        @rtype: 1|0
        """
        for i in items:
            if i in lst_or_dic:
                return 1
        return 0


    def requested( self, c, *keys ):
        """
        Determine the keys in an info dictionary of a complex that
        need to be calculated or updated

        @param c: Complex
        @type  c: Complex
        @param keys: key or keys for c.infos dict
        @type  keys: str OR [str]

        @return: 1 if the given value should be calculated for the
                 given complex
        @rtype: 1|0
        """
        ## force update
        if self.force and self.__containsAny( self.force, *keys):
            return 1

        ## fill empty or nonexisting fields
        for k in keys:
            if not self.force and c.get(k, None) is None:
                return 1

        return 0


    def pickleError( self, o ):
        """
        Pickle object to disc.

        @param o: object to pickle
        @type  o: any
        """
        try:
            fname = self.ferror + '_dat'
            if not os.path.exists( fname ):
                T.dump( o, fname )
        except Exception, why:
            f = open('ErrorReportError_ContactSlave','a')
            f.write('Could not pickle error infos\n')
            f.write( str( why ) )
            f.close()


    def calcContacts( self, soln, c ):
        """
        Calculate contact matrices and fraction of native contacts, residue-
        and atom-based, with different distance cutoffs.

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex
        """
        try:
            if self.requested(c, 'fnac_4.5') and self.c_ref_atom_4_5 != None:
                ## cache pairwise atom distances for following calculations
                contacts = c.atomContacts( 4.5, self.mask_rec, self.mask_lig,
                                           cache=1, map_back=0 )
                ref = N.ravel( self.c_ref_atom_4_5 )

                c['fnac_4.5'] = N.sum( N.ravel(contacts) * ref )\
                 / float( N.sum(ref))

            if self.requested(c, 'fnac_10') and self.c_ref_atom_10 != None:

                contacts = c.atomContacts( 10., self.mask_rec, self.mask_lig,
                                           cache=1, map_back=0 )

                ref = N.ravel( self.c_ref_atom_10 )
                c['fnac_10'] = N.sum( N.ravel(contacts) * ref ) \
                 / float( N.sum(ref))

            if self.requested(c, 'c_res_4.5') \
               or ( self.c_ref_res_4_5 != None \
                    and (self.requested(c,'fnrc_4.5','fnSurf_rec'))):

                res_cont = c.resContacts( 4.5,
                                          cache=self.requested(c, 'c_res_4.5'))

                if self.c_ref_res_4_5 != None \
                   and self.requested(c, 'fnrc_4.5' ):
                    ref = N.ravel( self.c_ref_res_4_5 )
                    c['fnrc_4.5'] = N.sum(N.ravel(res_cont)*ref) \
                     /float(N.sum(ref))

                if self.c_ref_res_4_5 != None \
                   and self.requested(c, 'fnSurf_rec'):
                    r, l = c.fractionNativeSurface(res_cont,
                                                   self.c_ref_res_4_5 )
                    c['fnSurf_rec'] = r
                    c['fnSurf_lig'] = l

        except:
            m1 = m2 = s = 0
            try:
                m1, m2, s = c.get('model1',0), c.get('model2',0),\
                  c.get('soln',0)
            except:
                pass
            self.reportError('contact error (r %i : l %i, #%i)'%\
                             (m1,m2,s), soln)
##             self.pickleError( {'com':c, 'mrec':self.mask_rec,
##                                'mlig':self.mask_lig } )


    def calcReducedContacts( self, soln, c ):
        """
        Get contact matrices and/or fnarc from reduced-atom models.

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex       
        """
        if not (self.reduced_recs and self.reduced_ligs):
            return

        if not self.requested(c,'c_ratom_10','fnarc_10'):
            return

        try:
            ## create Complex with same orientation but reduced coordinates
            red_rec = self.reduced_recs[ c.rec_model.source ]
            red_lig = self.reduced_ligs[ c.lig_model.source ]
            red_com = Complex( red_rec, red_lig, c.ligandMatrix )

            contacts = red_com.atomContacts( 10.0, cache=1 )

            if self.requested(c, 'c_ratom_10'):
                c['c_ratom_10'] = MU.packBinaryMatrix(contacts)

            if self.c_ref_ratom_10 is not None:
                ref = N.ravel( self.c_ref_ratom_10 )
                c['fnarc_10'] = N.sum( N.ravel(contacts) * ref )\
                 / float( N.sum(ref))

        except:
            self.reportError('reduced contacts error', soln)
##             self.pickleError({'com':c, 'red_recs':self.reduced_recs,
##                               'red_ligs':self.reduced_ligs})


    def calcInterfaceRms( self, soln, c ):
        """
        RMS between this and reference interface atoms after superposition.
          - rms_if_bb, considers backbone of interface residues (10 A cutoff)
            (same as used for CAPRI)
          - rms_if, considers all atoms of more tightly defined interf.
            residues (correlates better with fraction of native contacts)

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex
        """
        try:
            if self.requested(c, 'rms_if_bb') and self.ref_interface_model_bb:

                m = c.model().compress(self.mask_interface_bb)
                c['rms_if_bb'] = self.ref_interface_model_bb.rms( m )

            if self.requested(c, 'rms_if') and self.ref_interface_model:

                m = c.model().compress( self.mask_interface )
                c['rms_if'] = self.ref_interface_model.rms( m )

        except:
            self.reportError('interface rms error', soln)


##     def calcProsa( self, soln, c ):
##         """ProsaII energy score"""
##         if self.requested( c, 'eProsa'):
##             try:
##                 prosaE = c.prosaEnergy()
##                 c['eProsa'] = prosaE
##             except:
##                 self.reportError('Prosa Error', soln )
##                 c['eProsa'] = None


    def calcProsa( self, soln, c ):
        """
        Prosa2003 energy score

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex        
        """
#        import socket
        if self.requested( c, 'eProsa'):
            try:
                prosaE = c.prosa2003Energy()
                c['eProsa'] = prosaE
            except:
                c['eProsa'] = None
#                c['eProsa'] = socket.gethostname()
                self.reportError('Prosa Error', soln )


    def calcPairScore( self, soln, c ):
        """
        calculate contact pair score

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex        
        """
        if self.requested( c,'ePairScore'):
            try:
                pairScore = c.contPairScore(cutoff=6.0, log=self.log,
                                            verbose=self.verbose )
                c['ePairScore'] = pairScore
            except:
                c['ePairScore'] = None
                self.reportError('PairScore Error', soln )


    def calcConservation( self, soln, c, method ):
        """
        calculate conservation score

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex
        @param method: scoring method to use see L{Complex.conservationScore}
        @type  method: str       
        """
        if self.requested( c, method):
            try:
                c[method] = c.conservationScore( method, log=self.log,
                                                 verbose=self.verbose )
            except:
                self.reportError('Conservation score Error', soln)


    def calcFoldX( self, soln, c ):
        """
        calculate fold-X binding energies

        @param soln: solution number
        @type  soln: int
        @param c: Complex
        @type  c: Complex        
        """
        if self.requested( c, 'foldX' ):
            try:
                c['foldX'] = c.foldXEnergy()
            except:
                self.reportError('FoldX Error', soln)


    def go(self, cmplxDic):
        """
        Obtain contact matrix for all complexes.

        @param cmplxDic: dictionary of complexes::
                         {soln:Complex, soln:Complex, ...} 
        @type  cmplxDic: {int:Complex}

        @return: similar dictionary with Complex.info['soln'] as keys and
                 file names of matrices as value::
                 { soln : fname, soln : fname ....}
        @rtype: {int:str}

        """
        result = {}

        startTime = time.time()

        for soln, c in cmplxDic.items():
            T.flushPrint( "%i," % soln )

##             if not os.path.exists( T.absfile('~/debug.dic') ):
##                 T.dump( cmplxDic,  T.absfile('~/debug.dic') )

            self.calcContacts( soln, c )

            self.calcInterfaceRms( soln, c )

            self.calcReducedContacts( soln, c )

## TODO: Prosa will not run when called via conatacSlave, runs as it should when
##        called as c.prosa2003Energy() in the interpreter. What's wroong here?
##        For mow the Prosa calculation is skipped.
##
            self.calcProsa( soln, c )

            self.calcPairScore( soln, c )  ## uses res-contacts

            self.calcFoldX( soln, c ) ##uses rec/lig.info['foldX'] if available

            for method in ['cons_ent', 'cons_max', 'cons_abs']:
                self.calcConservation( soln, c, method )

            c['__version_contacter'] = self.version()

            result[ soln ] = c

            c.slim()

##             if not os.path.exists(T.absfile('~/debug_afterslave.dic') ):
##                 T.dump( cmplxDic,  T.absfile('~/debug_afterslave.dic') )


        print "\navg time for last %i complexes: %f s" %\
              ( len(cmplxDic), (time.time()-startTime)/len(cmplxDic))

        return result


import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test ContactSlave locally without running the master.

    This allows to test the master and slave without using
    external nodes. The master still requires a running PVM
    though.
    """

    TAGS = [ BT.PVM ]

    def prepare(self):
        import tempfile
        self.cl_out = tempfile.mktemp('_test.cl')

    def test_ContactSlave(self):
        """Dock.ContactSlave test (local)"""
        import os
        from Biskit.Dock.ContactMaster import ContactMaster

        ## load complex list (docking result) and reference complex
        lst = T.load( T.testRoot() + "/dock/hex/complexes.cl")
        lst = lst[:3]
        refcom = T.load( T.testRoot() + "/com/ref.complex")

        ## let ContactMaster prepare everything but don't run it
        self.master = ContactMaster( lst, verbose = self.local,
                                     log=self.log,
                                     refComplex = refcom,
                                     outFile = self.cl_out )

        jobs = self.master.data

        self.slave = ContactSlave()
        self.slave.initialize( self.master.getInitParameters(1) )

        if self.local or self.VERBOSITY > 2:
            self.log.writeln("Currently available info records (from hex):")
            self.log.writeln( repr(jobs[0].info.keys()) )
            self.log.writeln( "Calculating all scores for %i complexes..." \
                              % len(jobs) )

        self.result = self.slave.go( jobs )

        if self.local or self.VERBOSITY > 2:
            self.log.writeln("info records after contacting: ")
            self.log.writeln( repr(jobs[0].info.keys()) )
        if self.local:
            print "new scores are available in 'result[0-2].info'"

        ## verify fraction of native atom contacts for third complex
        self.assertAlmostEqual( self.result[2]['fnac_10'],
                                0.11533600168527491, 7 )

    def cleanUp(self):
        T.tryRemove(  self.cl_out )

## ## PROFILING:
## in slave window:
##     slave.stop()
##     import profile
##     profile.run( 'slave._go( slave.dic )', 'report.out' )

## ## Analyzing
##     import pstats
##     p = pstats.Stats('report.out')

##     ## long steps and methods calling them
##     p.sort_stats('cumulative').print_stats(20)
##     p.print_callers(0.1)

##     ## time consuming methods 
##     p.sort_stats('time').print_stats(10)


if __name__ == '__main__':

##     BT.localTest()

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = ContactSlave()
    slave.start()
