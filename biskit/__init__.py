##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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

__version__ = '3.0.0.a'

import logging

## public classes
try:
    ## default error handler
    from biskit.errorHandler import ErrorHandler
    EHandler = ErrorHandler()

    from biskit.logFile import LogFile, StdLog, ErrLog
    from biskit.errors import BiskitError

##     from Blast2Seq import Blast2Seq

##     from EDParser import EZDParser

    from biskit.pdbModel import PDBModel, PDBProfiles, PDBError
    from biskit.xplorModel import XplorModel

    from biskit.profileCollection import ProfileCollection, ProfileError
##     from ProfileMirror import ProfileMirror
    
    from biskit.pdbCleaner import PDBCleaner, CleanerError

##     from ModelList import ModelList
##     from CommandLine import CommandLine
    
    from .amberResidues import AmberResidueType, AmberPrepParser
    from .amberResidueLibrary import AmberResidueLibrary,\
                                           AmberResidueLibraryError
    
    from .atomCharger import AtomCharger
    

    from .pdbDope import PDBDope
##     from Ramachandran import Ramachandran
    
    from .colorspectrum import ColorSpectrum, ColorError, colorRange
    from .matrixPlot import MatrixPlot
    
    from biskit.core.localpath import LocalPath, LocalPathError
    from biskit.core.dictlist import DictList

    
## ## PVM-dependent modules

##     from QualMaster import QualMaster
##     from StructureMaster import StructMaster
##     from StructureSlave import StructureSlave
##     from TrajFlexMaster import TrajFlexMaster, FlexError

except Exception as why:
    logging.warning('Could not import all biskit modules: ' + repr(why))
    raise

## clean up namespace
del logging
