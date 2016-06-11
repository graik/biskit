##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
##
##

"""
Settings
========

This module provides Mod-global settings as fields. Throughout
Biskit.Mod a (environment-dependent) parameter such as, e.g., ssh_bin
can be addressed as:

  >>> import Biskit.Mod.settings as S
  >>> bin = S.ssh_bin

However, since a user should not be required to hack python modules,
ssh_bin is not actually defined in settings.py. Instead, the value is
taken from C{~/.biskit/settings_Mod.cfg} -- which should have an entry
like C{ssh_bin=/bin/ssh # comment}. If this entry (or the config file)
is not found, settings.py uses the default value from
C{biskit/Biskit/data/defaults/settings_Mod.cfg}.

If missing, the user configuration file C{~/.biskit/settings_Mod.cfg} is
created automatically during the startup of Biskit.Mod (i.e. for any
import). The auto-generated file only contains parameters for which
the default values don't seem to work (invalid paths or binaries).

See L{Biskit.SettingsManager}

Summary for Biskit users
------------------------
    If you want to change a biskit parameter, do so in
    C{~/.biskit/settings_Mod.cfg}

Summary for Biskit developpers
------------------------------
  If you want to create a new user-adjustable parameter, do so in
  C{biskit/Biskit/data/defaults/settings_Mod.cfg}.

Summary for all
---------------
  !Dont't touch C{settings.py}!
"""
import Biskit as B
import Biskit.tools as T
import Biskit.SettingsManager as M

import user, sys

__CFG_DEFAULT = T.dataRoot() + '/defaults/settings_Mod.cfg'
__CFG_USER    = user.home + '/.biskit/settings_Mod.cfg'

try:
    m = M.SettingsManager(__CFG_DEFAULT, __CFG_USER, createmissing=True  )

    m.updateNamespace( locals() )

except Exception, why:
    B.EHandler.fatal( 'Error importing Biskit.Mod settings' )



##############################
## Check environment variables
env = {}

net_env = {'http_proxy':'http://cache.pasteur.fr:8080'}

##either this set or .ncbi
net_ncbi = { 'BLASTDB': '/Bis/db/blastdb',
             'BLASTMAT' : '/Bis/shared/rh73/ncbi' }

env.update(net_env)
env.update(net_ncbi)


## .ncbi example:
## [NCBI]
## DATA=/home/Bis/shared/rh73/ncbi/data/

## [NET_SERV]
## SRV_CONN_MODE=FIREWALL
## SRV_HTTP_PROXY_HOST=cache.pasteur.fr
## SRV_HTTP_PROXY_PORT=8080


######################
## clean up name space

del B, T, M, user, sys
del __CFG_DEFAULT, __CFG_USER, m

################
## empty test ##
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, settings is always executed anyway."""
    pass
