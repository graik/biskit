##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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

## default error handler
from ErrorHandler import ErrorHandler
EHandler = ErrorHandler()

## public classes
try:
    from BisList import BisList, BisListError, ConditionError, AmbiguousMatch,\
         ItemNotFound
    from DictList import DictList

    from LogFile import LogFile, StdLog, ErrLog
    from Errors import BiskitError

    from ExeConfig import ExeConfig, ExeConfigError
    from ExeConfigCache import ExeConfigCache
    from Executor import Executor, TemplateError

    from AmberCrdParser import AmberCrdParser, ParseError
    from AmberRstParser import AmberRstParser
    from AmberParmBuilder import AmberParmBuilder
    from PDBCleaner import PDBCleaner, CleanerError
    from Blast2Seq import Blast2Seq
    from ChainCleaner import ChainCleaner
    from ChainSeparator import ChainSeparator
    from ChainWriter import ChainWriter

    from EnsembleTraj import EnsembleTraj
    from LocalPath import LocalPath, LocalPathError

    from PCRModel import PCRModel
    from PDBModel import PDBModel, PDBProfiles, PDBError

    from ProfileCollection import ProfileCollection, ProfileError
    from Prosa import ProsaII
    from Pymoler import Pymoler

    from TrajCluster import TrajCluster
    from Trajectory import Trajectory, TrajError, TrajProfiles
    from XplorInput import XplorInput, XplorInputError
    from Xplorer import Xplorer, XplorerError, RunError
    from ColorSpectrum import ColorSpectrum
    from MatrixPlot import MatrixPlot

    from Hmmer import Hmmer
    from Fold_X import Fold_X, Fold_XError
    from WhatIf import WhatIf, WhatIf_Error
    from SurfaceRacer import SurfaceRacer, SurfaceRacer_Error
    from DSSP import Dssp, Dssp_Error
    from FuzzyCluster import FuzzyCluster
    from PDBDope import PDBDope
    
    from ModelList import ModelList
    from Ramachandran import Ramachandran

## PVM-dependent modules

    from QualMaster import QualMaster
    from StructureMaster import StructMaster
    from StructureSlave import StructureSlave
    from TrajFlexMaster import TrajFlexMaster, FlexError

except Exception, why:
    EHandler.warning('Could not import all biskit modules:', trace=1 )

try:
    import biggles
    del biggles
except:
    EHandler.warning('Could not import biggles module -- plotting is not available.')
