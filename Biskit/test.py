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


##
## MODULE                 INPUT TEST DATA
##

## AmberCrdEntropist.py
## AmberCrdParser.py
## AmberEntropist.py
## AmberEntropyMaster.py
## AmberEntropySlave.py
## AmberParmBuilder.py
## AmberRstParser.py
## BisList.py
## Blast2Seq.py
## ChainCleaner.py        /rec/1A2P_rec_original.pdb
## ChainSeparator.py      /com/1BGS_original.pdb
## ChainWriter.py
## ColorSpectrum.py       None
## decorators.py
## default_hosts.py
## DictList.py
## difflib_old.py
## DSSP.py                com/1BGS.pdb
## EnsembleTraj.py
## ErrorHandler.py
## Errors.py
## ExeConfigCache.py
## ExeConfig.py
## Executor.py
## Fold_X.py
## FuzzyCluster.py        None
## gnuplot.py
## hist.py
## Hmmer.py               lig/1A19.pdb
## hosts.py
## IcmCad.py
## __init__.py
## LocalPath.py
## LogFile.py
## match2seq.py
## mathUtils.py
## MatrixPlot.py
## ModelList.py
## molTools.py
## molUtils.py
## msms.py
## PCRModel.py
## PDBCleaner.py
## PDBDope.py
## PDBModel.py
## plotUtils.py
## ProfileCollection.py
## Prosa2003.py		lig/1A19.pdb + rec/2A2P.pdb
## Prosa.py
## Pymoler.py
## QualSlave.py
## Ramachandran.py      /lig_pcr_00/traj.dat
## ReduceCoordinates.py
## rmsFit.py
## settings_default.py
## SettingsManager.py
## settings.py
## SparseArray.py
## StructureSlave.py
## SurfaceRacer.py
## surfaceRacerTools.py
## Table.py
## tools.py
## TrajCluster.py
## Trajectory.py
## TrajFlexMaster.py
## TrajFlexSlave.py
## WhatIf.py
## Xplorer.py
## XplorInput.py


## DOCK

## Dock/Analyzer.py
## Dock/ComplexEvolvingList.py
## Dock/ComplexEvolving.py
## Dock/ComplexList.py
## Dock/ComplexModelRegistry.py
## Dock/Complex.py
## Dock/ComplexRandomizer.py
## Dock/ComplexTraj.py
## Dock/ContactMaster.py
## Dock/ContactSlave.py
## Dock/Docker.py
## Dock/FixedList.py
## Dock/HexParser.py
## Dock/hexTools.py
## Dock/__init__.py
## Dock/settings_default.py
## Dock/settings.py

## MOD

## Mod/AlignerMaster.py
## Mod/__Aligner.py
## Mod/Aligner.py
## Mod/AlignerSlave.py
## Mod/Analyse.py
## Mod/Benchmark.py
## Mod/CheckIdentities.py
## Mod/__init__.py
## Mod/Modeller.py
## Mod/ModelMaster.py
## Mod/ModelSlave.py
## Mod/modUtils.py
## Mod/SequenceSearcher.py
## Mod/settings_default.py
## Mod/settings.py
## Mod/TemplateCleaner.py
## Mod/TemplateSearcher.py
## Mod/ValidationSetup.py

## PVM

## PVM/dispatcher.py
## PVM/ExampleMaster.py
## PVM/ExampleSlave.py
## PVM/__init__.py
## PVM/pvm.py
## PVM/PVMThread.py
## PVM/pvmTools.py
## PVM/Status.py
## PVM/TrackingJobMaster.py

## STATISTICS

## Statistics/Density.py
## Statistics/__init__.py
## Statistics/lognormal.py
## Statistics/pstat.py
## Statistics/stats.py

import unittest
import Numeric as N

class Biskit_Test( unittest.TestCase ):

    def test_FuzzyCluster( self ):
        from Biskit.FuzzyCluster import Test
        T = Test()
        self.assertEquals( N.shape(T.run()), T.expected_result() )

        
    def test_Ramachandran( self ):
        from Biskit.Ramachandran import Test
        T = Test()
        self.assertAlmostEquals( T.run(), T.expected_result(), 2 )


    def test_ColorSpectrum( self ):
        from Biskit.ColorSpectrum import Test
        T = Test()
        self.assertEquals( T.run(), T.expected_result() )


    def test_ChainCleaner( self ):
        from Biskit.ChainCleaner import Test
        T = Test()
        self.assertEquals( T.run(), T.expected_result() )

        
    def test_ChainSeparator( self ):
        from Biskit.ChainSeparator import Test
        T = Test()
        self.assertEquals( T.run(), T.expected_result() )


    def test_ChainWriter( self ):
        from Biskit.ChainWriter import Test
        T = Test()
        self.assertEquals( T.run(), T.expected_result() )    

        
class Biskit_Test_Applications( unittest.TestCase ):

    def test_Dssp( self ):
        from Biskit.DSSP import Test
        T = Test()
        self.assertEquals( T.run(), T.expected_result() )


    def test_Prosa2003( self ):
        from Biskit.Prosa2003 import Test
        T = Test()
        self.assertAlmostEquals( N.sum(T.run()),
                                 N.sum(T.expected_result()),
                                 2 )

    def test_Hmmer( self ):
        from Biskit.Hmmer import Test
        T = Test()
        self.assertAlmostEquals( N.sum(T.run()),
                                 N.sum(T.expected_result()),
                                 2 )

def suite_biskit():
    suite = unittest.TestSuite()
    suite.addTest( Biskit_Test( "test_FuzzyCluster" ) )
    suite.addTest( Biskit_Test( "test_Ramachandran" ) )
    suite.addTest( Biskit_Test( "test_ColorSpectrum" ) )
    suite.addTest( Biskit_Test( "test_ChainCleaner" ) )
    suite.addTest( Biskit_Test( "test_ChainSeparator" ) )
    suite.addTest( Biskit_Test( "test_ChainWriter" ) )
    return suite


        
def suite_biskit_applications():
    suite = unittest.TestSuite()
    suite.addTest( Biskit_Test_Applications( "test_Dssp" ) )
    suite.addTest( Biskit_Test_Applications( "test_Prosa2003" ) )
    suite.addTest( Biskit_Test_Applications( "test_Hmmer" ) )
    return suite


if __name__ == '__main__':

    f=open('/shared_bin/test.txt','w')
    
    runner = unittest.TextTestRunner(f, verbosity=2)
    runner.run( suite_biskit_applications() )
    runner.run( suite_biskit() )

    f.close()
    
    r = open('/shared_bin/test.txt','r')
    test_result = r.readlines()
    r.close()

    for line in test_result:
        if line != '\n':
            print line
