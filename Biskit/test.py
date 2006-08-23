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

## !! AmberCrdEntropist.py
## !! AmberCrdParser.py
## !! AmberEntropist.py
## !! AmberEntropyMaster.py
## !! AmberEntropySlave.py
## !! AmberParmBuilder.py
## !! AmberRstParser.py
##   N --BisList.py
## Blast2Seq.py           None
## ChainCleaner.py        rec/1A2P_rec_original.pdb
## ChainSeparator.py      com/1BGS_original.pdb
## ChainWriter.py         rec/1A2P_rec_original.pdb
## ColorSpectrum.py       None
## decorators.py          None
##   N -- default_hosts.py
## DictList.py            None
##   -- difflib_old.py
## DSSP.py                com/1BGS.pdb
## EnsembleTraj.py        lig_pcr_00/traj.dat
##   N -- ErrorHandler.py
##   N -- Errors.py
## ExeConfigCache.py      None
## ExeConfig.py           None
## Executor.py            None
## Fold_X.py              /rec/1A2P.pdb
## FuzzyCluster.py        None
## gnuplot.py             None
## hist.py                None
## Hmmer.py               lig/1A19.pdb
##   N -- hosts.py
## IcmCad.py              lig/1A19.pdb, lig_pcr_00/traj.dat
##   N -- __init__.py
## LocalPath.py           None
##   N -- LogFile.py
## match2seq.py           lig_pcr_00/traj.dat
## mathUtils.py           None
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
## !! QualSlave.py
## Ramachandran.py        lig_pcr_00/traj.dat
## ReduceCoordinates.py   com/1BGS.pdb
## rmsFit.py              lig_pcr_00/traj.dat
## !! settings_default.py
## !! SettingsManager.py
## !! settings.py
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

##  N -- Dock/Analyzer.py
## Dock/ComplexEvolvingList.py    dock/hex/complexes.cl"
## Dock/ComplexEvolving.py        com/ref.complex
## Dock/ComplexList.py            dock/hex/complexes.cl
## --not running -- Dock/ComplexModelRegistry.py   dock/hex/complexes.cl
## Dock/Complex.py                com/1BGS.psf com/lig.model com/rec.model
## Dock/ComplexRandomizer.py      rec/1A2P.pdb lig/1A19.pdb rec/1A2P.psf lig/1A19.psf
## Dock/ComplexTraj.py            com/1BGS.pdb
## Dock/ContactMaster.py          dock/hex/complexes.cl  com/ref.complex
## Dock/ContactSlave.py
## Dock/Docker.py                 multidock/lig/1A19_models.dic, multidock/rec/1A2P_model.dic
##  N -- Dock/FixedList.py              
## Dock/HexParser.py              dock/rec/1A2P_model.dic, dock/lig/1A19_model.dic, dock/hex/1A2P-1A19_hex.out
##  Dock/hexTools.py              com/1BGS.pdb
##  N -- Dock/__init__.py
##  N -- Dock/settings_default.py
##  N -- Dock/settings.py

## MOD

## Mod/Aligner.py              /Mod/project/templates/t_coffe/* (1DT7_A.alpha, 1K8U_A.alpha, 1MHO_.alpha, 1NSH_A.alpha, 1J55_A.alpha, 1KSO_A.alpha, 1MQ1_A.alpha, 1K2H_A.alpha, 1M31_A.alpha, 1MWN_A.alpha),  Mod/project/templates/templates.fasta, Mod/project/sequences/nr.fast, Mod/project/target.fasta

## Mod/AlignerMaster.py
## Mod/AlignerSlave.py
## Mod/Analyse.py
## Mod/Benchmark.py
## Mod/CheckIdentities.py
## Mod/__init__.py

## Mod/Modeller.py             Mod/project/templates/modeller/* (1DT7_A.pdb, 1K2H_A.pdb, 1KSO_A.pdb, 1MHO_.pdb, 1MWN_A.pdb, 1J55_A.pdb, 1K8U_A.pdb, 1M31_A.pdb, 1MQ1_A.pdb, 1NSH_A.pdb), Mod/project/t_coffee/final.pir_aln, Mod/project/target.fasta

## Mod/ModelMaster.py
## Mod/ModelSlave.py
## Mod/modUtils.py               Mod/project/target.fasta
## Mod/SequenceSearcher.py       Mod/project/target.fasta
## Mod/settings_default.py
## Mod/settings.py

## Mod/TemplateCleaner.py        Mod/project/templates/nr/* (1DT7.pdb, 1KSO.pdb, 1MWN.pdb, 1J55.pdb, 1M31.pdb, 1NSH.pdb, 1K2H.pdb, 1MHO.pdb, chain_index.txt, 1K8U.pdb, 1MQ1.pdb )

## Mod/TemplateSearcher.py       Mod/project/target.fasta
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


class Mod_Test( unittest.TestCase ):


    def test_SequenceSearcher( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.Mod.SequenceSearcher import Test
        t = Test()
        self.assertEquals( t.run(flavour='blastp'), t.expected_result() )
        self.assertEquals( t.run(flavour='blastpgp'), t.expected_result() )


    def test_TemplateSearcher( self ):
        """
        No testing done. Just let run trough
        """       
        from Biskit.Mod.TemplateSearcher import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_TemplateCleaner( self ):
        """
        No testing done. Just let run trough
        """       
        from Biskit.Mod.TemplateCleaner import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


##     def test_Aligner( self ):
##         """
##         No testing done. Just let run trough
##         """       
##         from Biskit.Mod.Aligner import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )
        

##     def test_Modeller( self ):
##         """
##         No testing done. Just let run trough
##         """       
##         from Biskit.Mod.Modeller import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


    def test_modUtils( self ):
        from Biskit.Mod.modUtils import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

                


class Dock_Test( unittest.TestCase ):


    def test_Complex( self ):
        from Biskit.Dock.Complex import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ComplexList( self ):
        """
        No testing done, will just display a mock plot
        """
        from Biskit.Dock.ComplexList import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_ComplexEvolving( self ):
        from Biskit.Dock.ComplexEvolving import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ComplexEvolvingList( self ):
        from Biskit.Dock.ComplexEvolvingList import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ComplexTraj( self ):
        from Biskit.Dock.ComplexTraj import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

                


#############################################
## Didn't get this running ...
##
##     def test_ComplexModelRegistry( self ):
##         """
##         No testing done. Just let run trough
##         """
##         from Biskit.Dock.ComplexModelRegistry import Test
##         t = Test()
##         self.assertEquals( t.run(inheret='t'), t.expected_result() )


    def test_HexParser( self ):
        from Biskit.Dock.HexParser import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_hexTools( self ):
        from Biskit.Dock.hexTools import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )
        
        


class Dock_Test_Applications( unittest.TestCase ):


##     def test_ComplexRandomizer( self ):
##         from Biskit.Dock.ComplexRandomizer import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


    def test_Docker( self ):
        """
        No testing done. Just let run trough
        NOTE: This test takes 15 minutes if run=1 (a HEX run)!!
        """
        from Biskit.Dock.Docker import Test
        t = Test()
        self.assertEquals( t.run( run=0 ), t.expected_result() )
        

        
class Pvm_Depending( unittest.TestCase ):
    
    
##     def test_TrajFlexMaster( self ):
##         """
##         No testing done, will just display a matrix plot
##         """
##         from Biskit.TrajFlexMaster import Test
##         t = Test()
##         self.assertEquals( t.run(), t.expected_result() )


    def test_ContactMaster( self ):
        from Biskit.Dock.ContactMaster import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )

                           



class Biskit_Test( unittest.TestCase ):

    def test_FuzzyCluster( self ):
        from Biskit.FuzzyCluster import Test
        t = Test()
        self.assertEquals( N.shape(t.run()), t.expected_result() )

        
    def test_Ramachandran( self ):
        from Biskit.Ramachandran import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 2 )


    def test_ColorSpectrum( self ):
        from Biskit.ColorSpectrum import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ChainCleaner( self ):
        from Biskit.ChainCleaner import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ChainSeparator( self ):
        from Biskit.ChainSeparator import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ChainWriter( self ):
        from Biskit.ChainWriter import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_Blast2Seq( self ):
        from Biskit.Blast2Seq import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )        


    def test_decorators( self ):
        from Biskit.decorators import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ReduceCoordinates( self ):
        from Biskit.ReduceCoordinates import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )
            

    def test_rmsFit( self ):
        from Biskit.rmsFit import Test
        t = Test()
        self.assertAlmostEquals( N.sum(N.ravel(t.run())),
                                 N.sum(N.ravel(t.expected_result())),
                                 4 )


    def test_SparseArray( self ):
        from Biskit.SparseArray import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_match2seq( self ):
        from Biskit.match2seq import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_Trajectory( self ):
        from Biskit.Trajectory import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_EnsembleTraj( self ):
        from Biskit.EnsembleTraj import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_MatrixPlot( self ):
        """
        No testing done, will just display a matrix plot
        """
        from Biskit.MatrixPlot import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ModelList( self ):
        from Biskit.ModelList import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_molTools( self ):
        from Biskit.molTools import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )
        

    def test_PCRModel( self ):
        from Biskit.PCRModel import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBCleaner( self ):
        from Biskit.PDBCleaner import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBModel( self ):
        from Biskit.PDBModel import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_plotUtils( self ):
        """
        No testing done, will just display a mock plot
        """
        from Biskit.plotUtils import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ProfileCollection( self ):
        from Biskit.ProfileCollection import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_gnuplot( self ):
        """
        No testing done, will just display a mock plot
        """
        from Biskit.gnuplot import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_tools( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.tools import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_TrajCluster( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.TrajCluster import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_XplorInput( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.XplorInput import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )        


    def test_DictList( self ):
        """
        No testing done. will just display a mock plot
        """
        from Biskit.DictList import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_LocalPath( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.LocalPath import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_ExeConfigCache( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.ExeConfigCache import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ExeConfig( self ):
        """
        No testing done. Just let run trough
        """
        from Biskit.ExeConfig import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_hist( self ):
        from Biskit.hist import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_mathUtils( self ):
        from Biskit.mathUtils import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )



class Biskit_Test_Applications( unittest.TestCase ):

    def test_Dssp( self ):
        from Biskit.DSSP import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_Prosa2003( self ):
        from Biskit.Prosa2003 import Test
        t = Test()
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2 )

    def test_Hmmer( self ):
        from Biskit.Hmmer import Test
        t = Test()
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2 )
        
    def test_Fold_X( self ):
        from Biskit.Fold_X import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_IcmCad( self ):
        from Biskit.IcmCad import Test
        t = Test()
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2  )


    def test_SurfaceRacer( self ):
        from Biskit.SurfaceRacer import Test
        t = Test()
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 4 )

    def test_Pymoler( self ):
        """
        A PyMol wondow should be shown briefly .. not a real test.
        """
        from Biskit.Pymoler import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )        


    def test_msms( self ):
        from Biskit.msms import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBDope( self ):
        """
        No testing done. Just let run trough.
        """
        from Biskit.PDBDope import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


    def test_WhatIf( self ):
        from Biskit.WhatIf import Test
        t = Test()
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_Executor( self ):
        """
        No testing done. Will open emacs breifly.
        """
        from Biskit.Executor import Test
        t = Test()
        self.assertEquals( t.run(), t.expected_result() )


#######
## Create test suits
        
def suite_biskit():
    suite = unittest.makeSuite( Biskit_Test )
    return suite

        
def suite_biskit_applications():
    suite = unittest.makeSuite( Biskit_Test_Applications )
    return suite


def suite_dock():
    suite = unittest.makeSuite( Dock_Test )
    return suite

        
def suite_dock_applications():
    suite = unittest.makeSuite( Dock_Test_Applications )
    return suite


def suite_mod():
    suite = unittest.makeSuite( Mod_Test )
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

##     f.write('\nTESTS FOR BISKIT MODULES CALLING EXTERNAL APPLICATIONS:\n')
##     runner.run( suite_biskit_applications() )
    
##     f.write('\nTESTS FOR BISKIT MODULES:\n')
##     runner.run( suite_biskit() )

##     f.write('\nTESTS FOR BISKIT.DOCK MODULES:\n')
##     runner.run( suite_dock() )

##     f.write('\nTESTS FOR BISKIT.DOCK MODULES CALLING EXTERNAL APPLICATIONS:\n')
##     runner.run( suite_dock_applications() )
    
##     f.write('\nTESTS FOR MODULES USING PVM:\n')
##     runner.run( pvm_depending() )


    f.write('\nTESTS FOR BISKIT.MOD MODULES:\n')
    runner.run( suite_mod() )
    

    f.close()

    
    report('test.log')
