##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
"""
High-level parallelisation with PVM.

Python jobs can be distributed from a TrackingJobMaster to many
JobSlaves running on different machines. The general mechanism is shown in
ExampleMaster.py and ExampleSlave.py 

The package depends on the pypvm_core c-library which has to be compiled
in the package directory. With Python versions >= 2.3, platform specific
compilations can be placed in subfolders. The name of this folder should
be assembled from: 'py' + version + '_' + platform.machine()
Where platform is a python module available since Python 2.3.
See Biskit.tools.platformFolder!
"""

##
## platform-dependent import of PVM core module
## Look for a folder ala py2.3_i686 
##
try:
    import Biskit.tools as T

    __path__.insert( 0, T.platformFolder( __path__[0] ) )

except ImportError, why:
    pass

##
## export public classes
##
from TrackingJobMaster import TrackingJobMaster
from dispatcher import JobSlave
