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

This module provides global settings for Biskit.Mod as fields. Not very
elegantly, it is simply a copy of Biskit.settings. The only difference
is that it takes its default parameters from L{Biskit.Mod.settings_default}
instead of L{Biskit.settings_default}. What follows is the description of
Biskit.settings:

Throughout Biskit a (environment-dependent) parameter such
as, e.g., ssh_bin can be addressed as:

  >>> import Biskit.settings as S
  >>> bin = S.ssh_bin

However, since a user should not be required to hack python modules, ssh_bin
is not actually defined in settings.py. Instead, settings_default.py defines
a parameter ssh_bin (with its default value and a comment), which
tells settings.py that it should export a parameter of that name. The value
of this parameter is taken from the Biskit configuration file
C{~/.biskit/settings.dat} -- which should have an entry like
C{ssh_bin=/bin/ssh  # comment}. If this entry (or the config file) is not
found, settings.py uses the default value from settings_default.py,
and informs the user that parameter ssh_bin is missing from the configuration
file.

The configuration file C{~/.biskit/settings.dat} is created with the
script C{scripts/Biskit/setup_env.py} and then modified by the user.

Summary for Biskit users
------------------------
  If you want to change a biskit parameter, do so in C{~/.biskit/settings.dat}

Summary for Biskit developpers
------------------------------
  If you want to create a new user-adjustable parameter, do so in
  C{settings_default.py}.

Summary for all
---------------
  !DON'T TOUCH C{settings.py}!
  (Unless you came up with a better way of organising the whole parameter
  system or have bug-fixes to make.)
"""
import Biskit as B
import Biskit.SettingsManager as M
import settings_default as D

#
## Parse settings_default.py and create one static field for each of its
## fields.
##
try:
    m = M.SettingsManager( defaults_module=D )

    m.updateNamespace( locals())

except Exception, why:
    B.EHandler.warning( 'Error importing %s' % D, trace=1 )


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
