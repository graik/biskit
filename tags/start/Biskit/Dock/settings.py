##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$

"""
This module provides global settings for Biskit.Dock as fields. Not very
elegantly, it is simply a copy of Biskit.settings. The only difference
is that it takes its default parameters from Biskit.Dock.settings_default
instead of Biskit.settings_default. What follows is the description of
Biskit.settings:

Throughout Biskit a
(environment-dependent) parameter such as, e.g., ssh_bin can be addressed as:
>>> import Biskit.settings as S
>>> bin = S.ssh_bin

However, since a user should not be required to hack python modules, ssh_bin
is not actually defined in settings.py. Instead, settings_default.py defines
a parameter ssh_bin (with its default value and comment), which
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

hex_env = {'HEX_ROOT':'/home/Bis/johan/APPLICATIONS/HEX',
           'HEX_CACHE':'/home/Bis/johan/APPLICATIONS/HEX/hex_cache',
           'HEX_VERSION':'4b'}

prosaII_env = {'PROSA_BASE':'/home/Bis/shared/rh73/prosa/prosabase/'}

env.update(hex_env)
env.update(prosaII_env)



