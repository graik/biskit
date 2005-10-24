#!/usr/bin/env python
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
##
##

## last $Author$
## last $Date$
## $Revision$

from Biskit import PDBModel, PDBDope, molUtils
import Biskit.tools as t

from HexParser import HexParser
from ComplexList import ComplexList
import hexTools 
import settings

import sys, os, re, commands
import os.path
from time import localtime, sleep
from threading import Thread, RLock, Condition
import Numeric as N

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

    Todo: implement that as JobMaster / JobSlave instead, so that also the
    Hex-file parsing is distributed.
    """

    def __init__( self, recDic, ligDic, recPdb=None, ligPdb=None,
                  comPdb=None, out='hex_%s', soln=512,
                  recChainId=None, ligChainId=None, macDock=None,
                  bin=settings.hex_bin ):
        """
        recDic, ligDic - model dictionary
        recPdb, ligPdb - str, hex-formatted PDB with rec/lig models
        recNumber, ligNumber - int, model number rec/lig in dict and hex PDB
        recChainId, ligChainId - list of str, force chain IDs only for HEX
        out   - str, out folder name, will be formatted against date
        macDock - force macro docking ON or OFF (default: size decides)
        soln - int, number of solutions to keep from HEX (default 512)
        """
        self.lock = RLock()
        self.lockMsg = Condition( self.lock )

        t.ensure( macDock, int, [None] )
        t.ensure( recDic, dict, )
        t.ensure( ligDic, dict, )
        
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

        self.result = ComplexList()

        self.soln = soln


    def prepareOutFolder( self, fout ):
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
        ids - list of str, len(ids) == m.lenChains
        set chaiID for Hex pdb files
        -> m is changed directly
        """
        if ids:
            ids = t.toList( ids )
            cMap = m.chainMap()
            
            for chain in range( m.lenChains() ):
                idx = N.nonzero( cMap == chain )
                for i in idx:
                    m.atoms[i]['chain_id'] = ids[chain]


    def prepareHexPdbs( self, modelDic, idList, subDir ):
        """
        create HEX-formatted PDB for each model
        subDir - str, 'rec' or 'lig'
        -> { modelNumber : file_name_created, .. }
        """
        result = {}
        for k in modelDic:
            
            m = modelDic[k].clone()

            fhex = self.out + '/%s/%s_%03d_hex.pdb' % (subDir, m.pdbCode, k)
            if os.path.exists( fhex ):
                print "using old ", os.path.split( fhex )[1]

            else:
                m.remove( m.maskH() )
                m = molUtils.sortAtomsOfModel(m)
                self.__setChainID( m, idList )

                hexTools.createHexPdb_single( m, fhex )
                print "created ", fhex

            result[ k ] = fhex

        return result


    def createHexInp( self, nRec, nLig ):
        """
        Create a HEX macro file for docking a rec lig pair.
        nRec, nLig - int, model number rec / lig in model dictionaries
        -> (str, str), macro file name, future out file name
        !! raises DockerError, if macro dock option is different from previous
           call
        """
        ## fetch PDB file names for given model numbers
        recPdb, ligPdb = self.recHexPdbs[nRec], self.ligHexPdbs[nLig]

        ## fetch the according PDBModels
        rec, lig = self.recDic[ nRec ], self.ligDic[ nLig ]

        fout_base = self.out+'/%s_%i-%s_%i' % (rec.pdbCode, nRec,
                                               lig.pdbCode, nLig)

        if os.path.exists( t.absfile( fout_base + '_hex.mac' ) ):
            fmac = t.absfile( fout_base + '_hex.mac' )
            fout = t.absfile( fout_base + '_hex.out' )
            macro = self.macroDock
            print "Dock setup: using old ", os.path.split( fmac )[1]
        else:
            fmac, fout, macro = hexTools.createHexInp(recPdb,rec, ligPdb,lig,
                                   self.comPdb, outFile=fout_base,
                                   macDock=self.macroDock, sol=self.soln)

        if self.macroDock == None:
            self.macroDock = macro
        else:
            if self.macroDock != macro:
                raise DockerError('MacroDock setting changed to %i: %i : %i' %\
                                  ( macro, nRec, nLig) )

        ## remember which models and which out file are used
        self.runDic[ fmac ] = ( nRec, nLig, fout )

        return fmac, fout


    def runHex( self, finp, log=0, ncpu=2, nice=0, host=os.uname()[1] ):

        flog = finp + '.log'

        cmd = "%s -ncpu %i -nice %i -noexec < %s > %s"
        cmd = cmd % (self.bin, ncpu, nice, finp, flog )
        
        cmd = "ssh %s %s" % (host, cmd)

        runner = RunThread( cmd, self, finp=finp, host=host, log=log )

        self.lock.acquire()

        self.lastCmd = cmd
        self.running += [runner]

        runner.start()

        self.lock.release()


    def countDifferentModels( self, lst ):

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
    Todo: implement that as JobSlave instead
    """

    def __init__( self, cmd, owner, finp=None, log=0,
                  host=None, **kw ):
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


    def run( self ):

        try:
            if not os.path.exists( self.fout ):

                print "Executing: ", self.host, ' ', t.stripFilename(self.finp)

                cmd_lst = self.cmd.split()

                self.status = os.spawnvp(os.P_WAIT, cmd_lst[0], cmd_lst )

                waited = 0
                while waited < 25 and not os.path.exists( self.fout ):
                    sleep( 5 )
                    waited += 5

            if self.status != 0:
                raise DockerError, 'Hex returned exit status %i' % self.status

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
        self.owner.doneHex( self )

    def failed( self ):
        print "FAILED: ", self.host, ' ', t.stripFilename(self.finp)
        print "\tJob details:"
        print "\tCommand: ", self.cmd
        print "\tinput:   ", self.finp
        print "\tHex log: ", self.log
        print "\tHex out: ", self.fout
        print
        print "\t", t.lastError()
        
        errorRunner = self
        self.owner.failedHex( self )


    def __hackHexOut( self, nRec, nLig, fout ):
        """
        Replace current model numbers in HEX out file by given nRec and nLig.
        Only to be used with one-against-one dockings because all
        ReceptorModel and LigandModel entries are replaced by nRec, nLig.
        nRec, nLig - int
        fout - str, name of existing  out file, is overwritten
        """
        old_out = open( fout, 'r' )
        contents = old_out.readlines()

        oldRec = 'ReceptorModel: \d+$'
        oldLig = 'LigandModel: \d+$'
        correctRec = 'ReceptorModel: %i' % nRec
        correctLig = 'LigandModel: %i'   % nLig
        contents = [re.sub(oldRec, correctRec, i ) for i in contents]
        contents = [re.sub(oldLig, correctLig, i ) for i in contents]
        old_out.close()

        new_out = open( fout, 'w')
        new_out.writelines(contents )
        new_out.close()


#### TEST ######
        
if __name__ == '__main__':

##     rec = testRoot() + '/dock/pcr_rec/1A2P_1_0.pdb.model'
##     lig = testRoot() + '/dock/pcr_lig/1A19_1_3.pdb.model'
##     com = testRoot() + '/dock/com/1BGS.model'
    testf = t.absfile('~/interfaces/c15/dock_multi_0919')

    recDic = t.Load( testf + '/pcr_rec/1AVV_models.dic' )
    ligDic = t.Load( testf + '/pcr_lig/1SHF_models.dic')
    
    d = Docker( recDic, ligDic,                
                out = '~/data/tmp/hex_%s' )

    fmac1, fout = d.createHexInp( 2, 2 )
    fmac2, fout2= d.createHexInp( 1, 1 )
    
    d.runHex( fmac1, log=1, ncpu=2 )
    d.runHex( fmac2, log=0, ncpu=2, host='riddler' )
    
    print "ALL jobs submitted."
