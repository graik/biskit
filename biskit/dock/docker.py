## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
Prepare and run a HEX docking.
"""

from biskit import PDBModel, PDBDope, molUtils
import biskit.tools as t

from biskit.dock.hexparser import HexParser
from biskit.dock import ComplexList
from biskit.dock import hextools 

import re
import os.path
import subprocess
from time import localtime, sleep
from threading import Thread, RLock, Condition
import biskit.core.oldnumeric as N0


class DockerError( Exception ):
    pass

class Docker:
    """
    Prepare and run a hex docking for one or several models
    of a receptor and ligand pair. Collect docking results
    into a growing ComplexList.

    The docking runs are started in seperate Threads. External objects
    can register a method with call_when_done() that is called whenever
    a docking has finished. The waitForLastHex() method can be called
    to wait until all currently running dockings are finished.

    The models in the 2 given dictionaries might get different chainIds but
    Docker doesn't save them.

    @todo: implement that as JobMaster / JobSlave instead, so that also the
           Hex-file parsing is distributed.
    """

    def __init__( self, recDic, ligDic, recPdb=None, ligPdb=None,
                  comPdb=None, out='hex_%s', soln=512,
                  recChainId=None, ligChainId=None, macDock=None,
                  bin='hex', verbose=1 ):  ## use ExeConfig instead of bin=...
        """
        @param recDic: receptor dictionary
        @type  recDic: dict
        @param ligDic: ligand dictionary
        @type  ligDic: dict
        @param recPdb: hex-formatted PDB with rec models
        @type  recPdb: str
        @param ligPdb: hex-formatted PDB with lig models
        @type  ligPdb: str
        @param comPdb: hex-formatted PDB with com models
        @type  comPdb: str
        @param soln: number of solutions to keep from HEX (default 512)
        @type  soln: int
        @param recChainId: force chain IDs only for HEX
        @type  recChainId: [str]
        @param ligChainId: force chain IDs only for HEX
        @type  ligChainId: [str]
        @param out: out folder name, will be formatted against date
        @type  out: str
        @param macDock: force macro docking ON or OFF (default: size decides)
        @type  macDock: 1|0
        @param bin: path to HEX binary
        @type  bin: str
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        """
        self.lock = RLock()
        self.lockMsg = Condition( self.lock )

        t.ensure( macDock, int, [None] )
        t.ensure( recDic, dict, )
        t.ensure( ligDic, dict, )

        self.verbose = verbose

        self.recDic = recDic
        self.ligDic = ligDic

        ## force chain IDs
        self.recChainId = recChainId
        self.ligChainId = ligChainId

        ## create folder if necessary
        self.out = self.prepareOutFolder( out )

        self.comPdb = None
        if comPdb:
            self.comPdb = t.absfile( comPdb )

        ## remember which HEX PDB was created for which model
        self.recHexPdbs = {}
        self.ligHexPdbs = {}

        ## remember which macrofile uses which models and which out file
        self.runDic = {}

        if recPdb:
            ## HEX PDB provided externally
            for k in self.recDic:
                self.recHexPdbs[ k ] = t.absfile( recPdb )
        else:
            ## HEX PDB has to be created
            self.recHexPdbs = self.prepareHexPdbs( self.recDic,
                                                   self.recChainId, 'rec')

        if ligPdb:
            for k in self.ligDic:
                self.ligHexPdbs[ k ] = t.absfile( recPdb )
        else:
            self.ligHexPdbs = self.prepareHexPdbs( self.ligDic,
                                                   self.ligChainId, 'lig')

        self.bin = bin

        self.call_when_done = None
        self.call_when_failed=None

        ## remember whether macro dock was used
        self.macroDock = macDock

        self.running = []
        self.hexDone = 0
        self.hexFailed = 0

        self.result = ComplexList()

        self.soln = soln


    def prepareOutFolder( self, fout ):
        """
        Setup an output folder

        @param fout: outfile name
        @type  fout: str

        @return: out path
        @rtype: str
        """
        ## date stamp in folder name
        try:
            datestr = "%02i%02i" % (localtime()[1], localtime()[2])
            fout = fout % ( datestr )
            if not os.path.exists( t.absfile( fout ) ):
                os.mkdir( t.absfile( fout ) )

        ## user provided folder name -> formatting fails
        except:
            if not os.path.exists( t.absfile( fout ) ):
                os.mkdir( t.absfile( fout ) )

        ## create folders for rec / lig HEX pdbs
        if not os.path.exists( t.absfile( fout ) + '/rec' ):
            os.mkdir( t.absfile( fout ) + '/rec' )

        if not os.path.exists( t.absfile( fout ) + '/lig' ):
            os.mkdir( t.absfile( fout ) + '/lig' )

        return t.absfile( fout )


    def __setChainID( self, m, ids ):
        """
        set chaiID for Hex pdb files

        @param m: model
        @type  m: PDBModel
        @param ids: chain id, len(ids) == m.lenChains
        @type  ids: [str]

        @return: m is changed directly
        @rtype: PDBModel
        """
        if ids:
            ids = t.toList( ids )
            cMap = m.chainMap()

            for chain in range( m.lenChains() ):
                idx = N0.nonzero( cMap == chain )
                for i in idx:
                    m.atoms['chain_id'][i] = ids[chain]


    def prepareHexPdbs( self, modelDic, idList, subDir ):
        """
        create HEX-formatted PDB for each model

        @param modelDic: model dictionary to be prepared for HEX
        @type  modelDic: dict
        @param idList: force chain IDs for HEX pdb files
        @type  idList: [str]
        @param subDir: 'rec' or 'lig'
        @type  subDir: str

        @return: model dictionary, { modelNumber : file_name_created, .. }
        @rtype: dict
        """
        result = {}
        for k in modelDic:

            m = modelDic[k].clone()

            fhex = self.out + '/%s/%s_%03d_hex.pdb' % (subDir, m.pdbCode, k)
            if os.path.exists( fhex ) and self.verbose:
                print("using old ", os.path.split( fhex )[1])

            else:
                m.remove( m.maskH() )
                m = molUtils.sortAtomsOfModel(m)
                self.__setChainID( m, idList )

                hextools.createHexPdb_single( m, fhex )
                if self.verbose: print("created ", fhex)

            result[ k ] = fhex

        return result


    def createHexInp( self, nRec, nLig ):
        """
        Create a HEX macro file for docking a rec lig pair.

        @param nRec: model number rec in model dictionaries
        @type  nRec: int
        @param nLig: model number lig in model dictionaries
        @type  nLig: int

        @return: macro file name, future out file name
        @rtype: (str, str)

        @raise DockerError: if macro dock option is different from previous
                            call
        """
        ## fetch PDB file names for given model numbers
        recPdb, ligPdb = self.recHexPdbs[nRec], self.ligHexPdbs[nLig]

        ## fetch the according PDBModels
        rec, lig = self.recDic[ nRec ], self.ligDic[ nLig ]

        fout_base = self.out+'/%s_%i-%s_%i' % (rec.pdbCode, nRec,
                                               lig.pdbCode, nLig)

        if os.path.exists( t.absfile( fout_base + '_hex.mac' ) ) and self.verbose:
            fmac = t.absfile( fout_base + '_hex.mac' )
            fout = t.absfile( fout_base + '_hex.out' )
            macro = self.macroDock
            print("Dock setup: using old ", os.path.split( fmac )[1])
        else:
            silent=1
            if self.verbose: silent=0

            fmac, fout, macro = hextools.createHexInp(recPdb,rec, ligPdb,lig,
                                                      self.comPdb, outFile=fout_base,
                                                      macDock=self.macroDock, sol=self.soln,
                                                      silent=silent)

        if self.macroDock is None:
            self.macroDock = macro
        else:
            if self.macroDock != macro:
                raise DockerError('MacroDock setting changed to %i: %i : %i' %\
                                  ( macro, nRec, nLig) )

        ## remember which models and which out file are used
        self.runDic[ fmac ] = ( nRec, nLig, fout )

        return fmac, fout


    def runHex( self, finp, log=0, ncpu=2, nice=0, host=os.uname()[1] ):
        """
        @param finp: hex macro file
        @type  finp: str
        @param log: write log file
        @type  log: 1|0
        @param ncpu: number of cpus to use
        @type  ncpu: int
        @param nice: nice value for HEX job (default: 0)
        @type  nice: int
        @param host: host to run jon on
        @type  host: str
        """
        flog = finp + '.log'

        cmd = "%s -ncpu %i -batch -nice %i -noexec -e %s -log %s"
        cmd = cmd % (self.bin, ncpu, nice, finp, flog )
        
        if not host in (os.uname()[1], 'localhost'):
            cmd = "ssh %s %s" % (host, cmd)

        runner = RunThread( cmd, self, finp=finp, host=host,
                            log=log, verbose=self.verbose )

        self.lock.acquire()

        self.lastCmd = cmd
        self.running += [runner]

        runner.start()

        self.lock.release()


    def countDifferentModels( self, lst ):
        """
        """
        mrec = []
        mlig = []
        for c in lst:
            if c.rec_model not in mrec:
                mrec += [ c.rec_model ]
            if c.lig_model not in mlig:
                mlig += [ c.lig_model ]

        return len( mrec ), len( mlig )


    def doneHex( self, runner ):
        """
        Do something after hex has finished. Notify all threads waiting
        on self.lockMsg
        """
        self.lock.acquire()

        ## append ComplexList
        self.result += runner.result

        if self.call_when_done:
            self.call_when_done( runner )

        self.running.remove( runner )

        self.hexDone += 1
        
        self.lockMsg.notifyAll()
        self.lock.release()


    def failedHex( self, runner ):
        """
        notify of failed hex run.
        """
        self.lock.acquire()

        self.running.remove( runner )
        self.hexFailed += 1

        if self.call_when_failed:
            self.call_when_failed( runner )

        self.lockMsg.notifyAll()
        self.lock.release()


    def waitForLastHex( self ):
        """
        Return after the last hex thread has finished.
        """
        self.lock.acquire()

        while len( self.running ) > 0:
            self.lockMsg.wait()

        self.lock.release()
        return


    def set_call_when_done( self, funct ):
        """
        funct( RunThread )
        """
        self.lock.acquire()
        self.call_when_done = funct
        self.lock.release()

    def set_call_when_failed( self, funct ):
        """
        funct( RunThread )
        """
        self.lock.acquire()
        self.call_when_failed = funct
        self.lock.release()


class RunThread( Thread ):
    """
    @todo: implement that as JobSlave instead
    """

    def __init__( self, cmd, owner, finp=None, log=0,
                  host=None, verbose=1, **kw ):
        """
        @param cmd: command to execute
        @type  cmd: str
        @param owner: Docker job to run to run
        @type  owner: Docker instance
        @param finp: hex macro file (default: None)
        @type  finp: str
        @param log: write log file (default: 0)
        @type  log: 1|0
        @param host: host to run jon on (default: None, localhost)
        @type  host: str
        @param kw: optional key/value pairs
        @type  kw: key=value
        """
        Thread.__init__( self, name=(host or 'local'), **kw )
        self.cmd = cmd
        self.owner = owner
        self.finp = finp
        self.log = log
        self.host = host or os.uname()[1]
        self.nRec, self.nLig, self.fout = self.owner.runDic[ self.finp ]
        self.output = None
        self.status = 0
        self.result = None
        self.verbose = verbose


    def run( self ):
        """
        Run HEX job.

        @raise DockerError: if HEX exists with error
        """
        try:
            if not os.path.exists( self.fout ):

                if self.verbose:
                    print("Executing on ", self.host, ' with ', \
                          t.stripFilename(self.finp))
                    print("Command: ", self.cmd)

                cmd_lst = self.cmd.split()

                ## self.status = os.spawnvp(os.P_WAIT, cmd_lst[0], cmd_lst )
                
                p = subprocess.Popen( cmd_lst, executable=cmd_lst[0],
                                      universal_newlines=True, 
                                      stdout=subprocess.DEVNULL, ## see flog
                                      )
                self.pid = p.pid
    
                output, error = p.communicate()
                self.status = p.returncode

                if self.status != 0:
                    raise DockerError('Hex returned exit status %i' % self.status)

                waited = 0
                while waited < 25 and not os.path.exists( self.fout ):
                    sleep( 5 )
                    waited += 5


            ## replace model numbers in HEX output file
            self.__hackHexOut( self.nRec, self.nLig, self.fout )

            parser = HexParser(self.fout, self.owner.recDic,
                               self.owner.ligDic)

            ## generate ComplexList from hex output
            self.result = parser.parseHex()

            self.done()

        except:
            self.failed()


    def done( self ):
        """
        HEX job done.
        """
        self.owner.doneHex( self )


    def failed( self ):
        """
        If HEX job fails
        """
        print("FAILED: ", self.host, ' ', t.stripFilename(self.finp))
        print("\tJob details:")
        print("\tCommand: ", self.cmd)
        print("\tinput:   ", self.finp)
        print("\tHex log: ", self.log)
        print("\tHex out: ", self.fout)
        print()
        print("\t", t.lastError())

        self.owner.failedHex( self )


    def __hackHexOut( self, nRec, nLig, fout ):
        """
        Replace current model numbers in HEX out file by given nRec and nLig.
        Only to be used with one-against-one dockings because all
        ReceptorModel and LigandModel entries are replaced by nRec, nLig.

        @param nRec: receptor number
        @type  nRec: int
        @param nLig: ligand number
        @type  nLig: int
        @param fout: name of existing  out file, is overwritten
        @type  fout: str
        """
        old_out = open( fout, 'r' )
        contents = old_out.readlines()

        oldRec = r'ReceptorModel: \d+$'
        oldLig = r'LigandModel: \d+$'
        correctRec = 'ReceptorModel: %i' % nRec
        correctLig = 'LigandModel: %i'   % nLig
        contents = [re.sub(oldRec, correctRec, i ) for i in contents]
        contents = [re.sub(oldLig, correctLig, i ) for i in contents]
        old_out.close()

        new_out = open( fout, 'w')
        new_out.writelines(contents )
        new_out.close()



#############
##  TESTING        
#############
import biskit.test as BT

class TestCore(BT.BiskitTest):
    """Base class for short and long Test case"""

    def prepare(self):
        import tempfile
        self.out_folder = tempfile.mktemp('_test_docker_%s')

    def cleanUp(self):
        t.tryRemove( self.out_folder, tree=1 )

    def dry_or_wet_run( self, run=True ):
        """ """
        import time, os
        import os.path

        ligDic = t.load(t.testRoot('/dock/lig/1A19_model.dic'))
        recDic = t.load(t.testRoot('/dock/rec/1A2P_model.dic'))

        self.d = Docker( recDic, ligDic, out=self.out_folder,
                         verbose=self.local  )

        # dock rec 1 vs. lig 2 on localhost
        fmac1, fout = self.d.createHexInp( 1, 1 )
        if run:
            self.d.runHex( fmac1, log=1, ncpu=6 )
            if self.local: print("ALL jobs submitted.")

            self.d.waitForLastHex()
            
            self.assertEqual(self.d.hexFailed, 0, 'Hex exited with error')
            self.assertTrue(len(self.d.result) > 1, 'no results collected')

            ## dock receptor 1 vs. ligand one on remote host
            # fmac2, fout2= d.createHexInp( 1, 1 )
            # d.runHex( fmac2, log=0, ncpu=2, host='remote_host_name' )
    


class TestLong(TestCore):
    """Test case running a complete hex docking (ca. 30 min/CPU)"""

    TAGS = [ BT.EXE, BT.LONG ]

    def test_DockerLong(self):
        """Dock.Docker real test (20-30min/CPU)"""
        self.dry_or_wet_run( run=True )

class TestShort(TestCore):

    TAGS = [ BT.EXE ]

    def test_DockerShort(self):
        """Dock.Docker dry run test"""
        self.dry_or_wet_run( run=False )

if __name__ == '__main__':

    BT.localTest(debug=False)
