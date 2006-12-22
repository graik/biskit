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
## last $Date  $
## $Revision$


#############################################################
##                                                         ##
##            MODULES THAT DON'T HAVE ANT TEST             ##
##                                                         ##
#############################################################

######################## BISKIT #############################
## !! AmberCrdEntropist.py   |
## !! AmberCrdParser.py      |
## !! AmberEntropist.py      |
## !! AmberEntropyMaster.py  |- I can't test these. No Amber.
## !! AmberEntropySlave.py   |
## !! AmberParmBuilder.py    |
## !! AmberRstParser.py      |

##  settings.py          |
##  BisList.py           |
##  default_hosts.py     |- Nothing to test
##  Errors.py            |
##  hosts.py             |
##  __init__.py          |

## difflib_old.py          |
## Table.py -- outdated -- |- Not our modules or should be removed
## Prosa.py -- outdated -- |

############################# DOCK ############################
## Dock/__init__.py
## Dock/settings.py

############################# MOD #############################
## Mod/settings.py
## Mod/__init__.py

############################# PVM #############################
## PVM/__init__.py
## PVM/dispatcher.py
## PVM/pvm.py
## PVM/PVMThread.py
## PVM/pvmTools.py
## PVM/Status.py
## PVM/TrackingJobMaster.py

########################### STATISTICS #########################
## Statistics/__init__.py
## Statistics/pstat.py
## Statistics/stats.py



#############################################################
##                                                         ##
##   MODULES WITH TESTS (FOLLOWED BY TEST INPUT DATA)      ##
##                                                         ##
#############################################################

## (V) - verbose module, try to fix
## (P) - displays plot also when not executed localy
## (D) - depends on another module calling an external application
## (E) - is not using the Executor

###### BISKIT ########
## SettingsManager.py     None
## Blast2Seq.py           None
## ChainCleaner.py        rec/1A2P_rec_original.pdb
## ChainSeparator.py      com/1BGS_original.pdb
## ChainWriter.py         rec/1A2P_rec_original.pdb
## ColorSpectrum.py       None
## decorators.py          None
## DictList.py            None
## DSSP.py                com/1BGS.pdb
## EnsembleTraj.py (P)    lig_pcr_00/traj.dat
## ExeConfigCache.py      None
## ErrorHandler.py        None
## ExeConfig.py           None
## Executor.py            None
## Fold_X.py              /rec/1A2P.pdb
## FuzzyCluster.py        None
## gnuplot.py (P)         None
## hist.py                None
## Hmmer.py               lig/1A19.pdb
## IcmCad.py              lig/1A19.pdb, lig_pcr_00/traj.dat
## LocalPath.py           None
## LogFile.py             None
## match2seq.py           lig_pcr_00/traj.dat
## molUtils.py            lig/1A19.pdb
## mathUtils.py           None
## MatrixPlot.py (P)      None
## ModelList.py           rec/1A2P.pdb, lig/1A19.pdb, com/1BGS.pdb
## molTools.py            lig/1A19.pdb
## msms.py                lig/1A19.pdb
## PCRModel.py            com/1BGS.psf com/1BGS.pdb rec/1A2P.psf rec/1A2P.pdb
## PDBCleaner.py          rec/1A2P_rec_original.pdb
## PDBDope.py (D)         com/1BGS.pdb
## PDBModel.py            rec/1A2P.pdb
## plotUtils.py (P)       None
## ProfileCollection.py   None
## Prosa2003.py           lig/1A19.pdb + rec/2A2P.pdb
## Pymoler.py             lig_pcr_00/traj.dat
## Ramachandran.py (D)    lig_pcr_00/traj.dat
## ReduceCoordinates.py   com/1BGS.pdb
## rmsFit.py              lig_pcr_00/traj.dat
## SparseArray.py         None
## SurfaceRacer.py        lig/1A19.pdb
## surfaceRacerTools.py   lig/1A19.pdb
## tools.py               rec/1A2P.pdb
## TrajCluster.py         lig_pcr_00/traj.dat
## Trajectory.py          lig_pcr_00/traj.dat
## TrajFlexMaster.py      lig_pcr_00/traj.dat
## TrajFlexSlave.py       -- tested via the master --
## WhatIf.py              com/1BGS.pdb
## XplorInput.py          None
## Xplorer.py (E)         lig/1A19.pdb, lig/1A19.psf
## StructureMaster.py     lig/1A19.pdb, rec/1A2P.pdb, com/1BGS.pdb
## StructureSlave.py      -- tested via the master --
## QualMaster.py          lig_pcr_00/traj.dat
## QualSlave.py           -- tested via the master --

############################# DOCK ############################
## Analyzer.py                    rec/1A2P.pdb, dock/hex/complexes.cl,
##                                   lig_pcr_00/traj.dat, com/ref.complex
##                                   dock/lig/1A19.model, lig/1A19_dry.model
##                                   dock/rec/1A2P.model, rec/1A2P_dry.model
## Dock/ComplexModelRegistry.py   dock/hex/complexes.cl
## Dock/ComplexEvolvingList.py    dock/hex/complexes.cl
## Dock/ComplexEvolving.py        com/ref.complex
## Dock/ComplexList.py            dock/hex/complexes.cl
## Dock/Complex.py (D)            com/1BGS.psf com/lig.model com/rec.model
## Dock/ComplexRandomizer.py      rec/1A2P.pdb, lig/1A19.pdb,
##                                  rec/1A2P.psf, lig/1A19.psf
## Dock/ComplexTraj.py            com/1BGS.pdb
## Dock/ContactMaster.py (D)      dock/hex/complexes.cl  com/ref.complex
## Dock/ContactSlave.py           -- tested via the master --
## Dock/Docker.py (E)             multidock/lig/1A19_models.dic,
##                                  multidock/rec/1A2P_model.dic
##                                  multidock/lig/1A19_15_7.model
##                                  multidock/lig/1A19_45_8.mod
## Dock/HexParser.py              dock/rec/1A2P_model.dic,
##                                  dock/lig/1A19_model.dic,
##                                  dock/hex/1A2P-1A19_hex.out
## Dock/hexTools.py               com/1BGS.pdb
## Dock/FixedList.py              None

############################# MOD #############################
## Mod/Aligner.py              Mod/project/templates/t_coffe/*
##                               (1DT7_A.alpha, 1K8U_A.alpha, 1MHO_.alpha,
##                                1NSH_A.alpha, 1J55_A.alpha, 1KSO_A.alpha,
##                                1MQ1_A.alpha, 1K2H_A.alpha, 1M31_A.alpha,
##                                1MWN_A.alpha),
##                                Mod/project/templates/templates.fasta,
##                                Mod/project/sequences/nr.fasta,
##                                Mod/project/target.fasta
## Mod/Modeller.py (E)        Mod/project/templates/modeller/*
##                                (1DT7_A.pdb, 1K2H_A.pdb, 1KSO_A.pdb,
##                                 1MHO_.pdb, 1MWN_A.pdb, 1J55_A.pdb,
##                                 1K8U_A.pdb, 1M31_A.pdb, 1MQ1_A.pdb,
##                                 1NSH_A.pdb),
##                                Mod/project/t_coffee/final.pir_aln,
##                                Mod/project/target.fasta
## Mod/TemplateCleaner.py      Mod/project/templates/nr/*
##                                (1DT7.pdb, 1KSO.pdb, 1MWN.pdb, 1J55.pdb,
##                                 1M31.pdb, 1NSH.pdb, 1K2H.pdb, 1MHO.pdb,
##                                 chain_index.txt, 1K8U.pdb, 1MQ1.pdb )
## Mod/TemplateSearcher.py     Mod/project/target.fasta
## Mod/Analyse.py              Mod/project/modeller/PDBModels.list,
##                                Mod/project/modeller/target.B99990001.pdb
##                                Mod/project/t_coffee/final.pir_aln,
##                                Mod/project/validation/,
##                                Mod/project/validation/1DT7/modeller/Modeller_Score.out,
##                                Mod/project/validation/1DT7/identities_cov.out,
##                                Mod/project/validation/1DT7/benchmark/rmsd_aa.out,
##                                Mod/project/validation/1DT7/benchmark/rmsd_ca.out,
##                                Mod/project/validation/1DT7/benchmark/PDBModels.list,
##                                Mod/project/validation/1DT7/modeller/target.B99990001.pdb
##                                Mod/project/validation/1J55/modeller/Modeller_Score.out,
##                                Mod/project/validation/1J55/identities_cov.out,
##                                Mod/project/validation/1J55/benchmark/rmsd_aa.out,
##                                Mod/project/validation/1J55/benchmark/rmsd_ca.out,
##                                Mod/project/validation/1J55/benchmark/PDBModels.list
##                                Mod/project/validation/1J55/modeller/target.B99990001.pdb
## Mod/Benchmark.py            Mod/project/validation/1DT7/modeller/PDBModels.list,
##                                Mod/project/validation/1DT7/reference.pdb
##                                Mod/project/validation/1DT7/modeller/target.B999900??.pdb
##                                  (target.B99990001.pdb, target.B99990002.pdb
##                                   target.B99990003.pdb, target.B99990004.pdb
##                                   target.B99990005.pdb, target.B99990006.pdb
##                                   target.B99990007.pdb, target.B99990008.pdb
##                                   target.B99990009.pdb, target.B99990010.pdb)
## Mod/CheckIdentities.py      Mod/project/t_coffee/final.pir_aln
## Mod/AlignerMaster.py        same as Aligner.py
## Mod/AlignerSlave.py         -- tested via the master --
## Mod/ModelMaster.py          same as Modell.py
## Mod/ModelSlave.py           -- tested via the master --
## Mod/modUtils.py             Mod/project/target.fasta
## Mod/SequenceSearcher.py     Mod/project/target.fasta
## Mod/ValidationSetup.py      Mod/project/templates/modeller/*
##                              (1DT7_A.pdb, 1K2H_A.pdb, 1KSO_A.pdb,
##                               1MHO_.pdb, 1MWN_A.pdb, 1J55_A.pdb,
##                               1K8U_A.pdb, 1M31_A.pdb, 1MQ1_A.pdb,
##                               1NSH_A.pdb),
##                             Mod/project/templates/nr/* (1DT7.pdb,
##                               1KSO.pdb, 1MWN.pdb, 1J55.pdb, 1M31.pdb,
##                               1NSH.pdb, 1K2H.pdb, 1MHO.pdb, chain_index.txt,
##                               1K8U.pdb, 1MQ1.pdb )
     

############################# PVM ############################
## PVM/ExampleMaster.py   None
## PVM/ExampleSlave.py    -- tested via the master --

############################ STATISTICS #####################
## Statistics/Density.py       None
## Statistics/lognormal.py     None


#############################################################
##                                                         ##
##       CHECKED IN TEST DATA  ( in T.testRoot() )         ##
##                                                         ##
#############################################################

## rec/1A2P_rec_original.pdb
## rec/1A2P.pdb
## rec/1A2P.psf
## rec/1A2P_dry.model
## lig/1A19.pdb
## lig/1A19.psf
## lig/1A19_dry.model
## lig_pcr_00/traj.dat
## com/1BGS_original.pdb
## com/1BGS.pdb
## com/1BGS.psf 
## com/lig.model
## com/rec.model
## com/ref.complex

## dock/hex/complexes.cl
## dock/rec/1A2P_model.dic
## dock/rec/1A2P.model 
## dock/lig/1A19_model.dic
## dock/lig/1A19.model 


## dock/hex/1A2P-1A19_hex.out
## multidock/lig/1A19_models.dic
## multidock/rec/1A2P_model.dic
## multidock/lig/1A19_15_7.model
## multidock/lig/1A19_45_8.model

## Mod/project/templates/t_coffe/* (1DT7_A.alpha, 1K8U_A.alpha, 1MHO_.alpha,
##                                   1NSH_A.alpha, 1J55_A.alpha, 1KSO_A.alpha,
##                                   1MQ1_A.alpha, 1K2H_A.alpha, 1M31_A.alpha,
##                                   1MWN_A.alpha),
## Mod/project/templates/templates.fasta
## Mod/project/sequences/nr.fasta
## Mod/project/target.fasta
## Mod/project/templates/modeller/* (1DT7_A.pdb, 1K2H_A.pdb, 1KSO_A.pdb,
##                                   1MHO_.pdb, 1MWN_A.pdb, 1J55_A.pdb,
##                                   1K8U_A.pdb, 1M31_A.pdb, 1MQ1_A.pdb,
##                                   1NSH_A.pdb),
## Mod/project/t_coffee/final.pir_aln
## Mod/project/templates/nr/* (1DT7.pdb, 1KSO.pdb, 1MWN.pdb, 1J55.pdb,
##                             1M31.pdb, 1NSH.pdb, 1K2H.pdb, 1MHO.pdb,
##                             chain_index.txt, 1K8U.pdb, 1MQ1.pdb )
## Mod/project/modeller/PDBModels.list
## Mod/project/modeller/target.B99990001.pdb

## Mod/project/validation/1DT7/modeller/Modeller_Score.out
## Mod/project/validation/1DT7/modeller/PDBModels.list
## Mod/project/validation/1DT7/reference.pdb
## Mod/project/validation/1DT7/identities_cov.out
## Mod/project/validation/1DT7/benchmark/rmsd_aa.out
## Mod/project/validation/1DT7/benchmark/rmsd_ca.out
## Mod/project/validation/1DT7/benchmark/PDBModels.list
## Mod/project/validation/1DT7/modeller/target.B999900??.pdb
##                                  (target.B99990001.pdb, target.B99990002.pdb
##                                   target.B99990003.pdb, target.B99990004.pdb
##                                   target.B99990005.pdb, target.B99990006.pdb
##                                   target.B99990007.pdb, target.B99990008.pdb
##                                   target.B99990009.pdb, target.B99990010.pdb)
## Mod/project/validation/1J55/modeller/Modeller_Score.out
## Mod/project/validation/1J55/identities_cov.out
## Mod/project/validation/1J55/benchmark/rmsd_aa.out
## Mod/project/validation/1J55/benchmark/rmsd_ca.out
## Mod/project/validation/1J55/benchmark/PDBModels.list
## Mod/project/validation/1J55/modeller/target.B99990001.pdb

 



import unittest
import Numeric as N
import Biskit.tools as T

import string

#############################################################

class Biskit( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit that don't
    call external applications.

    Modules that are dependent of PVM have been excluded and
    can be found in L{Biskit_Pvm} and modules that call external
    applications can be found in L{Biskit_Applications}.    
    """

    def test_SettingsManager( self ):
        """
        Testing L{Biskit.SettingsManager}
        
        No testing done. Just let run trough
        """
        from Biskit.SettingsManager import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_FuzzyCluster( self ):
        """
        Testing L{Biskit.FuzzyCluster}
        """         
        from Biskit.FuzzyCluster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( N.shape(t.run()), t.expected_result() )

        
    def test_Ramachandran( self ):
        """
        Testing L{Biskit.Ramachandran}
        """           
        from Biskit.Ramachandran import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 2 )


    def test_ColorSpectrum( self ):
        """
        Testing L{Biskit.ColorSpectrum}
        """    
        from Biskit.ColorSpectrum import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ChainCleaner( self ):
        """
        Testing L{Biskit.ChainCleaner}
        """    
        from Biskit.ChainCleaner import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ChainSeparator( self ):
        """
        Testing L{Biskit.ChainSeparator}
        """    
        from Biskit.ChainSeparator import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ChainWriter( self ):
        """
        Testing L{Biskit.ChainWriter}
        """
        from Biskit.ChainWriter import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Blast2Seq( self ):
        """
        Testing L{Biskit.Blast2Seq}
        """
        from Biskit.Blast2Seq import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )        


    def test_decorators( self ):
        """
        Testing L{Biskit.decorators}
        """
        from Biskit.decorators import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ReduceCoordinates( self ):
        """
        Testing L{Biskit.ReduceCoordinates}
        """
        from Biskit.ReduceCoordinates import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
            

    def test_rmsFit( self ):
        """
        Testing L{Biskit.rmsFit}
        """
        from Biskit.rmsFit import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( N.sum(N.ravel(t.run())),
                                 N.sum(N.ravel(t.expected_result())),
                                 4 )

    def test_SparseArray( self ):
        """
        Testing L{Biskit.SparseArray}
        """
        from Biskit.SparseArray import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_match2seq( self ):
        """
        Testing L{Biskit.match2seq}
        """
        from Biskit.match2seq import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Trajectory( self ):
        """
        Testing L{Biskit.Trajectory}
        """
        from Biskit.Trajectory import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_EnsembleTraj( self ):
        """
        Testing L{Biskit.EnsembleTraj}
        """
        from Biskit.EnsembleTraj import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_ErrorHandler( self ):
        """
        Testing L{Biskit.ErrorHandler}
        No testing done. Just let run trough
        """      
        from Biskit.ErrorHandler import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

            
    def test_MatrixPlot( self ):
        """
        Testing L{Biskit.MatrixPlot}
        No testing done, will just display a matrix plot
        """
        from Biskit.MatrixPlot import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ModelList( self ):
        """
        Testing L{Biskit.ModelList}
        """
        from Biskit.ModelList import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_molTools( self ):
        """
        Testing L{Biskit.molTools}
        """
        from Biskit.molTools import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 3 )


    def test_molUtils( self ):
        """
        Testing L{Biskit.molUtils}
        """
        from Biskit.molUtils import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )      


    def test_PCRModel( self ):
        """
        Testing L{Biskit.PCRModel}
        """
        from Biskit.PCRModel import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBCleaner( self ):
        """
        Testing L{Biskit.PDBCleaner}
        """
        from Biskit.PDBCleaner import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBModel( self ):
        """
        Testing L{Biskit.PDBModel}
        """
        from Biskit.PDBModel import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


    def test_plotUtils( self ):
        """
        Testing L{Biskit.plotUtils}
        No testing done, will just display a mock plot
        """
        from Biskit.plotUtils import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ProfileCollection( self ):
        """
        Testing L{Biskit.ProfileCollection}
        """
        from Biskit.ProfileCollection import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_gnuplot( self ):
        """
        Testing L{Biskit.gnuplot}        
        No testing done, will just display a mock plot
        """
        from Biskit.gnuplot import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_tools( self ):
        """
        Testing L{Biskit.tools}        
        No testing done. Just let run trough
        """
        from Biskit.tools import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_surfaceRacerTools( self ):
        """
        Testing L{Biskit.surfaceRacerTools}        
        No testing done. Just let run trough
        """
        from Biskit.surfaceRacerTools import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )
        

    def test_TrajCluster( self ):
        """
        Testing L{Biskit.TrajCluster}  
        No testing done. Just let run trough
        """
        from Biskit.TrajCluster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_XplorInput( self ):
        """
        Testing L{Biskit.XplorInput}  
        No testing done. Just let run trough
        """
        from Biskit.XplorInput import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )        


    def test_DictList( self ):
        """
        Testing L{Biskit.DictList}  
        No testing done. will just display a mock plot
        """
        from Biskit.DictList import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_LocalPath( self ):
        """
        Testing L{Biskit.LocalPath}  
        No testing done. Just let run trough
        """
        from Biskit.LocalPath import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_LogFile( self ):
        """
        Testing L{Biskit.LogFile} 
        No testing done. Just let run trough
        """
        from Biskit.LogFile import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ExeConfigCache( self ):
        """
        Testing L{Biskit.ExeConfigCache} 
        No testing done. Just let run trough
        """
        from Biskit.ExeConfigCache import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ExeConfig( self ):
        """
        Testing L{Biskit.ExeConfig} 
        No testing done. Just let run trough
        """
        from Biskit.ExeConfig import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_hist( self ):
        """
        Testing L{Biskit.hist} 
        """
        from Biskit.hist import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_mathUtils( self ):
        """
        Testing L{Biskit.mathUtils} 
        """
        from Biskit.mathUtils import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


#############################################################

class Biskit_Applications( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit that
    call external applications.

    Modules that are dependent of PVM has been excluded
    and can be found in L{Biskit_Pvm}.    
    """
    
##     def test_Xplorer( self ):
##         """
##         Testing L{Biskit.Xplorer}
##         No testing done. Just let run trough
##         """
##         from Biskit.Xplorer import Test
##         t = Test()
##         print 'Testing: ', t.__module__
##         self.assertEquals( t.run(), t.expected_result() )


    def test_Dssp( self ):
        """
        Testing L{Biskit.DSSP} 
        """
        from Biskit.DSSP import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Prosa2003( self ):
        """
        Testing L{Biskit.Prosa2003} 
        """
        from Biskit.Prosa2003 import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2 )

    def test_Hmmer( self ):
        """
        Testing L{Biskit.Hmmer} 
        """
        from Biskit.Hmmer import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2 )
        
    def test_Fold_X( self ):
        """
        Testing L{Biskit.Fold_X} 
        """
        from Biskit.Fold_X import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_IcmCad( self ):
        """
        Testing L{Biskit.IcmCad} 
        """
        from Biskit.IcmCad import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 2  )


    def test_SurfaceRacer( self ):
        """
        Testing L{Biskit.SurfaceRacer} 
        """
        from Biskit.SurfaceRacer import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( N.sum(t.run()),
                                 N.sum(t.expected_result()),
                                 4 )

    def test_Pymoler( self ):
        """
        Testing L{Biskit.Pymoler} 
        A PyMol wondow should be shown briefly .. not a real test.
        """
        from Biskit.Pymoler import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )        


    def test_msms( self ):
        """
        Testing L{Biskit.msms} 
        """
        from Biskit.msms import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_PDBDope( self ):
        """
        Testing L{Biskit.PDBDope} 
        No testing done. Just let run trough.
        """
        from Biskit.PDBDope import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_WhatIf( self ):
        """
        Testing L{Biskit.WhatIf} 
        """        
        from Biskit.WhatIf import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )


    def test_Executor( self ):
        """
        Testing L{Biskit.Executor} 
        No testing done. Will open emacs breifly.
        """
        from Biskit.Executor import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


#############################################################
        
class Biskit_Pvm( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit that depend on PVM.
    To run these tests the PVM deamon must be running.

    Modules that that call external  applications has
    been excluded and can be found in L{Biskit_Applications}.
    The remaining Biskit test are found in L{Biskit}.
    """

    
    def test_TrajFlexMaster( self ):
        """
        Testing L{Biskit.TrajFlexMaster}
        No testing done, will just display a matrix plot
        """
        from Biskit.TrajFlexMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_StructMaster( self ):
        """
        Testing L{Biskit.StructureMaster}
        No testing done. Just let run trough
        """       
        from Biskit.StructureMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_QualMaster( self ):
        """
        Testing L{Biskit.QualMaster}
        No testing done. Just let run trough
        """       
        from Biskit.QualMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

    
############################################################# 

class Dock( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Dock that don't
    call external applications.

    Modules that are dependent of PVM has been excluded and
    can be found in L{Dock_Pvm} and modules that call external
    applications can be found in L{Dock_Applications}.       
    """
    
    def test_Complex( self ):
        """
        Testing L{Biskit.Dock.Complex} 
        """
        from Biskit.Dock.Complex import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ComplexList( self ):
        """
        Testing L{Biskit.Dock.ComplexList}
        No testing done, will just display a mock plot
        """
        from Biskit.Dock.ComplexList import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_ComplexEvolving( self ):
        """
        Testing L{Biskit.Dock.ComplexEvolving} 
        """
        from Biskit.Dock.ComplexEvolving import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ComplexEvolvingList( self ):
        """
        Testing L{Biskit.Dock.ComplexEvolvingList} 
        """
        from Biskit.Dock.ComplexEvolvingList import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ComplexTraj( self ):
        """
        Testing L{Biskit.Dock.ComplexTraj} 
        """
        from Biskit.Dock.ComplexTraj import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_ComplexModelRegistry( self ):
        """
        Testing L{Biskit.Dock.ComplexModelRegistry}
        No testing done. Just let run trough
        """      
        from Biskit.Dock.ComplexModelRegistry import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_HexParser( self ):
        """
        Testing L{Biskit.Dock.HexParser} 
        """
        from Biskit.Dock.HexParser import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_hexTools( self ):
        """
        Testing L{Biskit.Dock.hexTools} 
        """
        from Biskit.Dock.hexTools import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )

        
    def test_FixedList( self ):
        """
        Testing L{Biskit.Dock.FixedList}
        No testing done. Just let run trough
        """    
        from Biskit.Dock.FixedList import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )   


    def test_Analyzer( self ):
        """
        Testing L{Biskit.Dock.Analyzer} 
        """
        from Biskit.Dock.Analyzer import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

#############################################################

class Dock_Applications( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Dock that 
    call external applications.

    Modules that are dependent of PVM has been excluded and
    can be found in L{Dock_Pvm}.    
    """
    
    def test_ComplexRandomizer( self ):
        """
        Testing L{Biskit.Dock.ComplexRandomizer} 
        """        
        from Biskit.Dock.ComplexRandomizer import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Docker( self ):
        """
        Testing L{Biskit.Dock.Docker} 
        No testing done. Just let run trough
        @note: This test takes 15 minutes if run=1 (a HEX run)!!
        """
        from Biskit.Dock.Docker import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run( run=0 ), t.expected_result() )


#############################################################

class Dock_Applications_Long( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Dock that 
    call external applications.

    Modules that are dependent of PVM has been excluded and
    can be found in L{Dock_Pvm}.    
    """

    def test_Docker( self ):
        """
        Testing L{Biskit.Dock.Docker} 
        No testing done. Just let run trough
        @note: This test takes 15 minutes if run=1 (a HEX run)!!
        """
        from Biskit.Dock.Docker import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run( run=1 ), t.expected_result() )

        
#############################################################
        
class Dock_Pvm( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Dock that depend on PVM.
    To run these tests the PVM deamon must be running.

    Modules that that call external  applications has
    been excluded and can be found in L{Dock_Applications}.
    The remaining Biskit.Dock test are found in L{Dock}.
    """

    def test_ContactMaster( self ):
        """
        Testing L{Biskit.Dock.ContactMaster}
        """           
        from Biskit.Dock.ContactMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 6 )


#############################################################

class Mod( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Mod. Most of
    these modules call external applications.

    Modules that are dependent of PVM has been excluded and
    can be found in L{Pvm_Depending}.
    """

    def test_SequenceSearcher( self ):
        """
        Testing L{Biskit.Mod.SequenceSearcher}
        No testing done. Just let run trough
        """
        from Biskit.Mod.SequenceSearcher import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(flavour='blastp'), t.expected_result() )
        self.assertEquals( t.run(flavour='blastpgp'), t.expected_result() )


    def test_TemplateSearcher( self ):
        """
        Testing L{Biskit.Mod.TemplateSearcher}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.TemplateSearcher import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_TemplateCleaner( self ):
        """
        Testing L{Biskit.Mod.TemplateCleaner}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.TemplateCleaner import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Aligner( self ):
        """
        Testing L{Biskit.Mod.Aligner}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.Aligner import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

    def test_Modeller( self ):
        """
        Testing L{Biskit.Mod.Modeller}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.Modeller import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_modUtils( self ):
        """
        Testing L{Biskit.Mod.modUtils}
        """           
        from Biskit.Mod.modUtils import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_ValidationSetup( self ):
        """
        Testing L{Biskit.Mod.ValidationSetup}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.ValidationSetup import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_CheckIdentities( self ):
        """
        Testing L{Biskit.Mod.CheckIdentities}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.CheckIdentities import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_Benchmark( self ):
        """
        Testing L{Biskit.Mod.Benchmark}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.Benchmark import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )       


    def test_Analyse( self ):
        """
        Testing L{Biskit.Mod.Analyse}
        No testing done. Just let run trough
        """       
        from Biskit.Mod.Analyse import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        

#############################################################


class Mod_Long( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Mod. Most of
    these modules call external applications.

    Modules that are dependent of PVM has been excluded and
    can be found in L{Pvm_Depending}.
    """

    def test_Aligner( self ):
        """
        Testing L{Biskit.Mod.Aligner}
        No testing done. Just let run trough
        @note: this test will take a minute or two
        """       
        from Biskit.Mod.Aligner import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(run=1), t.expected_result() )
        

    def test_Modeller( self ):
        """
        Testing L{Biskit.Mod.Modeller}
        No testing done. Just let run trough
        @note: this test will take a few minutes to finish
        """       
        from Biskit.Mod.Modeller import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(run=1), t.expected_result() )

        
#############################################################

        
class Mod_Pvm( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Mod that depend on PVM.
    To run these tests the PVM deamon must be running.

    The remaining Biskit.Mod test are found in L{Mod}.
    """

    def test_AlignerMaster( self ):
        """
        Testing L{Biskit.Mod.AlignerMaster}
        No testing done. Just let run trough
        @note: This test takes 4 minutes
        """       
        from Biskit.Mod.AlignerMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )

        
    def test_ModelMaster( self ):
        """
        Testing L{Biskit.Mod.ModelMaster}
        No testing done. Just let run trough
        @note: This test takes 10 minutes
        """       
        from Biskit.Mod.ModelMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


#############################################################

        
class Mod_Pvm_Long( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Mod that depend on PVM.
    To run these tests the PVM deamon must be running.

    The remaining Biskit.Mod test are found in L{Mod}.
    """

    def test_AlignerMaster( self ):
        """
        Testing L{Biskit.Mod.AlignerMaster}
        No testing done. Just let run trough
        @note: This test takes 4 minutes
        """       
        from Biskit.Mod.AlignerMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run( run=1 ), t.expected_result() )

        
    def test_ModelMaster( self ):
        """
        Testing L{Biskit.Mod.ModelMaster}
        No testing done. Just let run trough
        @note: This test takes 10 minutes
        """       
        from Biskit.Mod.ModelMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run( run=1 ), t.expected_result() )


#############################################################
        
class Statistics( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/Statistics.
    """

    def test_Density( self ):
        """
        Testing L{Biskit.Statistics.Density}
        No testing done. Just let run trough
        """
        from Biskit.Statistics.Density import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )


    def test_lognormal( self ):
        """
        Testing L{Biskit.Statistics.lognormal}
        No testing done. Just let run trough
        """
        from Biskit.Statistics.lognormal import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertAlmostEquals( t.run(), t.expected_result(), 8 )

        
#############################################################

class Pvm( unittest.TestCase ):
    """
    Tests for modules in ~biskit/Biskit/PVM
    
    Tests for modules that depend on PVM. To run these
    tests the PVM deamon must be running.
    """    

    def test_ExampleMaster( self ):
        """
        Testing L{Biskit.PVM.ExampleMaster}
        No testing done. Just let run trough
        """       
        from Biskit.PVM.ExampleMaster import Test
        t = Test()
        print 'Testing: ', t.__module__
        self.assertEquals( t.run(), t.expected_result() )
        


