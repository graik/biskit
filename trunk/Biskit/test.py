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
##   N --BisList.py
## Blast2Seq.py           None
## ChainCleaner.py        rec/1A2P_rec_original.pdb
## ChainSeparator.py      com/1BGS_original.pdb
## ChainWriter.py         rec/1A2P_rec_original.pdb
## ColorSpectrum.py       None
## decorators.py          None
##   N -- default_hosts.py
##     -- DictList.py
## difflib_old.py
## DSSP.py                com/1BGS.pdb
## EnsembleTraj.py        lig_pcr_00/traj.dat
##   N -- ErrorHandler.py
##   N -- Errors.py
## ExeConfigCache.py
## ExeConfig.py
##     -- Executor.py
## Fold_X.py              /rec/1A2P.pdb
## FuzzyCluster.py        None
## gnuplot.py             None
##   N -- hist.py
## Hmmer.py               lig/1A19.pdb
##   N -- hosts.py
## IcmCad.py              lig/1A19.pdb, lig_pcr_00/traj.dat
##   N -- __init__.py
##     -- LocalPath.py
##   N -- LogFile.py
## match2seq.py           lig_pcr_00/traj.dat
##   N -- mathUtils.py
## MatrixPlot.py          None
## ModelList.py           rec/1A2P.pdb, lig/1A19.pdb, com/1BGS.pdb
## molTools.py            lig/1A19.pdb
##   N -- molUtils.py
## msms.py                lig/1A19.pdb
## PCRModel.py            com/1BGS.psf com/1BGS.pdb rec/1A2P.psf rec/1A2P.pdb
## PDBCleaner.py          rec/1A2P_rec_original.pdb
## PDBDope.py             com/1BGS.pdb
## PDBModel.py            rec/1A2P.pdb
## plotUtils.py           None
## ProfileCollection.py   None
## Prosa2003.py           lig/1A19.pdb + rec/2A2P.pdb
##    -- Prosa.py -- outdated --
## Pymoler.py             lig_pcr_00/traj.dat
## QualSlave.py
## Ramachandran.py        lig_pcr_00/traj.dat
## ReduceCoordinates.py   com/1BGS.pdb
## rmsFit.py              lig_pcr_00/traj.dat
## settings_default.py
## SettingsManager.py
## settings.py
## SparseArray.py         None
##   N -- StructureSlave.py 
## SurfaceRacer.py        lig/1A19.pdb
##   N -- surfaceRacerTools.py
##     -- Table.py -- outdated --
## tools.py               rec/1A2P.pdb
## TrajCluster.py         lig_pcr_00/traj.dat
## Trajectory.py          lig_pcr_00/traj.dat
## TrajFlexMaster.py      lig_pcr_00/traj.dat
##   --  TrajFlexSlave.py
## WhatIf.py              com/1BGS.pdb
##   N -- Xplorer.py
## XplorInput.py          None


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


class Pvm_Depending( unittest.TestCase ):
    
    def test_TrajFlexMaster( self ):
        """
        No testing done, will just display a matrix plot
        """
        from Biskit.TrajFlexMaster import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


                           
class Biskit_Test( unittest.TestCase ):

##     def test_FuzzyCluster( self ):
##         from Biskit.FuzzyCluster import Test
##         t = Test()
##         self.assertEquals( N.shape(t.run()), t.expected_result() )

        
##     def test_Ramachandran( self ):
##         from Biskit.Ramachandran import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 2 )


##     def test_ColorSpectrum( self ):
##         from Biskit.ColorSpectrum import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_ChainCleaner( self ):
##         from Biskit.ChainCleaner import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )

        
##     def test_ChainSeparator( self ):
##         from Biskit.ChainSeparator import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_ChainWriter( self ):
##         from Biskit.ChainWriter import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_Blast2Seq( self ):
##         from Biskit.Blast2Seq import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )        


##     def test_decorators( self ):
##         from Biskit.decorators import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_ReduceCoordinates( self ):
##         from Biskit.ReduceCoordinates import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )
            

##     def test_rmsFit( self ):
##         from Biskit.rmsFit import Test
##         t = Test()
##         self.assertAlmostEquals( N.sum(N.ravel(t.run())),
##                                  N.sum(N.ravel(t.expected_result())),
##                                  4 )


##     def test_SparseArray( self ):
##         from Biskit.SparseArray import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_match2seq( self ):
##         from Biskit.match2seq import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_Trajectory( self ):
##         from Biskit.Trajectory import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


##     def test_EnsembleTraj( self ):
##         from Biskit.EnsembleTraj import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


##     def test_MatrixPlot( self ):
##         """
##         No testing done, will just display a matrix plot
##         """
##         from Biskit.MatrixPlot import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_ModelList( self ):
##         from Biskit.ModelList import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_molTools( self ):
##         from Biskit.molTools import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 8 )
        

##     def test_PCRModel( self ):
##         from Biskit.PCRModel import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


##     def test_PDBCleaner( self ):
##         from Biskit.PDBCleaner import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


##     def test_PDBModel( self ):
##         from Biskit.PDBModel import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


##     def test_plotUtils( self ):
##         """
##         No testing done, will just display a mock plot
##         """
##         from Biskit.plotUtils import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_ProfileCollection( self ):
##         from Biskit.ProfileCollection import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )
        

##     def test_gnuplot( self ):
##         """
##         No testing done, will just display a mock plot
##         """
##         from Biskit.gnuplot import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_tools( self ):
##         """
##         No testing done. Just let run trough
##         """
##         from Biskit.tools import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_TrajCluster( self ):
##         """
##         No testing done. Just let run trough
##         """
##         from Biskit.TrajCluster import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


    def test_XplorInput( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.XplorInput import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )        



class Biskit_Test_Applications( unittest.TestCase ):

    def test_Dssp( self ):
        from Biskit.DSSP import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


##     def test_Prosa2003( self ):
##         from Biskit.Prosa2003 import Test
##         t = Test()
##         self.assertAlmostEquals( N.sum(t.run()),
##                                  N.sum(t.expected_result()),
##                                  2 )

##     def test_Hmmer( self ):
##         from Biskit.Hmmer import Test
##         t = Test()
##         self.assertAlmostEquals( N.sum(t.run()),
##                                  N.sum(t.expected_result()),
##                                  2 )
        
##     def test_Fold_X( self ):
##         from Biskit.Fold_X import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


##     def test_IcmCad( self ):
##         from Biskit.IcmCad import Test
##         t = Test()
##         self.assertAlmostEquals( N.sum(t.run()),
##                                  N.sum(t.expected_result()),
##                                  2  )


##     def test_SurfaceRacer( self ):
##         from Biskit.SurfaceRacer import Test
##         t = Test()
##         self.assertAlmostEquals( N.sum(t.run()),
##                                  N.sum(t.expected_result()),
##                                  4 )

##     def test_Pymoler( self ):
##         """
##         A PyMol wondow should be shown briefly .. not a real test.
##         """
##         from Biskit.Pymoler import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )        


##     def test_msms( self ):
##         from Biskit.msms import Test
##         t = Test()
##         self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


##     def test_PDBDope( self ):
##         """
##         No testing done. Just let run trough.
##         """
##         from Biskit.PDBDope import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


    def test_WhatIf( self ):
        from Biskit.WhatIf import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


#######
## Create test suits
        
def suite_biskit():
    suite = unittest.makeSuite( Biskit_Test )
    return suite


        
def suite_biskit_applications():
    suite = unittest.makeSuite( Biskit_Test_Applications )
    return suite


def pvm_depending():
    suite = unittest.makeSuite( Pvm_Depending )
    return suite


def report( file ):
    r = open('test.log','r')
    test_result = r.readlines()
    r.close()

    print  '\n' + '='*60+'\n======= TEST RESULTS\n' + '='*60 +'\n'
    for line in test_result:
        if line != '\n':
            print line[:-1]
            

if __name__ == '__main__':


    f=open('test.log','w')
    
    runner = unittest.TextTestRunner(f, verbosity=2)

    f.write('\nTESTS FOR MODULES CALLING EXTERNAL APPLICATIONS:\n')
    runner.run( suite_biskit_applications() )
    
    f.write('\nTESTS FOR BISKIT MODULES:\n')
    runner.run( suite_biskit() )
    
#    f.write('\nTESTS FOR MODULES USING PVM:\n')
#    runner.run( pvm_depending() )

    f.close()

    
    report('test.log')
