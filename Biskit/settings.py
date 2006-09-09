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

"""
Settings
========

This module provides global settings as fields. Throughout Biskit a
(environment-dependent) parameter such as, e.g., ssh_bin can be addressed as:

  >>> import Biskit.settings as S
  >>> bin = S.ssh_bin

However, since a user should not be required to hack python modules,
ssh_bin is not actually defined in settings.py. Instead, the value is
taken from C{~/.biskit/settings.cfg} -- which should have an entry
like C{ssh_bin=/bin/ssh # comment}. If this entry (or the config file)
is not found, settings.py uses the default value from
C{biskit/external/defaults/settings.cfg}.

If missing, the user configuration file C{~/.biskit/settings.cfg} is
created automatically during the startup of Biskit (i.e. for any
import). The auto-generated file only contains parameters for which
the default values don't seem to work (invalid paths or binaries).

See L{Biskit.SettingsManager}

Summary for Biskit users
------------------------
  If you want to change a biskit parameter, do so in C{~/.biskit/settings.cfg}

Summary for Biskit developpers
------------------------------
  If you want to create a new user-adjustable parameter, do so in
  C{biskit/external/defaults/settings.cfg}.

Summary for all
---------------
  !Dont't touch C{settings.py}!
"""
import Biskit as B
import Biskit.tools as T
import Biskit.SettingsManager as M

import user, sys

__CFG_DEFAULT = T.projectRoot() + '/external/defaults/settings.cfg'
__CFG_USER    = user.home + '/.biskit/settings.cfg'

try:
    m = M.SettingsManager(__CFG_DEFAULT, __CFG_USER, createmissing=True  )

    m.updateNamespace( locals() )

except Exception, why:
    B.EHandler.fatal( 'Error importing Biskit settings')

##
## Create some settings on the fly
##
python_bin = sys.executable
xterm_bin  = T.absbinary('xterm')


###################################
## required environment variables.
## format: ENV_NAME : path_example

## Todo: These need to go to their Exe_*.dat files

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

######################
## clean up name space

del B, T, M, user, sys
