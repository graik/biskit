##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$

"""
This module provides global settings as fields. Throughout Biskit a
(environment-dependent) parameter such as, e.g., ssh_bin can be addressed as:
>>> import Biskit.settings as S
>>> bin = S.ssh_bin

However, since a user should not be required to hack python modules, ssh_bin
is not actually defined in settings.py. Instead, settings_default.py defines
a parameter ssh_bin (with its default value and a comment), which
tells settings.py that it should export a parameter of that name. The value
of this parameter is taken from the Biskit configuration file
.biskit/settings.dat -- which should have an entry like
'ssh_bin=/bin/ssh  # comment'. If this entry (or the config file) is not
found, settings.py uses the default value from settings_default.py,
and informs the user that parameter ssh_bin is missing from the configuration
file.

The configuration file .biskit/settings.dat is created with the
script scripts/Biskit/setup_env.py and then modified by the user.

Summary for Biskit users:
If you want to change a biskit parameter, do so in .biskit/settings.dat

Summary for Biskit developpers:
If you want to create a new user-adjustable parameter, do so in
settings_default.py.

Summary for all:
!DON'T TOUCH settings.py!
(Unless you came up with a better way of organising the whole parameter
system or have bug-fixes to make.)
"""
import sys

import Biskit as B
import Biskit.tools as T
import Biskit.SettingsManager as M
import Biskit.settings_default as D

##
## Parse settings_default.py and create one static field for each of its
## fields.
##
try:

    this = locals()   ## get pointer to name space of this module

    m = M.SettingsManager( defaults_module=D )
    m.updateNamespace( this )
    
except Exception, why:
    B.EHandler.warning( 'Error importing %s' % D, trace=1 )

##
## Create some settings on the fly
##
python_bin = sys.executable
xterm_bin  = T.absbinary('xterm')


###################################
## required environment variables.
## format: ENV_NAME : path_example

env = {}

pymol_env = {'PYMOL_PATH':'/home/Bis/shared/rh73/pymol-cvs-20021115',
             'PYMOL_EXTLIBPATH':'/home/Bis/shared/rh73/lib/'}
	     
pvm_env = {'PVM_ROOT':'~/pvm3',
           'PVM_RSH':'/usr/bin/ssh'}

blast_env = {'BLASTDB':'/home/Bis/raik/data/prog/blast/db',
             'BLASTMA':'/home/Bis/johan/APPLICATIONS/blast'}
	     
amber_env = {'AMBERHOME_8':'/Bis/shared/rh73/amber8_intel-7.1'}

prosaII_env = {'PROSA_BASE':'/home/Bis/shared/rh73/prosa/prosabase/'}

env.update(pvm_env)
env.update(blast_env)
env.update(amber_env)
env.update(prosaII_env)

