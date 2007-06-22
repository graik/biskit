##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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
"""
High-level parallelisation with PVM.

Python jobs can be distributed from a TrackingJobMaster to many
JobSlaves running on different machines. The general mechanism is shown in
ExampleMaster.py and ExampleSlave.py 

The compilation of PVM/pypvm can be tricky on some architectures. In order
to support installations that don't need parallelisation, only a warning is
issued if pypvm is missing and the public classes are exported as
Pseudo-classes. Pseudo classes are empty and raise an ImportError when you
try to initialize them. See also L{Biskit.tools.tryImport}.
"""

##
## error-tolerant export of public classes
##
from Biskit import EHandler
import Biskit.tools as T

r = True

r = T.tryImport( 'TrackingJobMaster', 'TrackingJobMaster', namespace=globals())
r = T.tryImport( 'dispatcher', 'JobSlave', namespace=globals() ) and r

if not r:
    EHandler.warning('Could not import PVM (Parallel Virtual Machine) modules.'+
        ' Please check that PVM and pypvm are installed!\n'+
        '\tParallelisation is not available.')

##
## clean up
##
del r, EHandler, T
