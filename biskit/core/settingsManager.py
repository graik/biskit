##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
Manage Biskit settings.
"""

import biskit as B

## allow relative imports when calling module as script for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.core; __package__ = "biskit.core"

from .. import tools as T
from . import settingsParser as P

import os

class WriteCfgError( P.SettingsError ):
    pass

class SettingsManager:
    """
    SettingsManager merges the parameters from a default and a user
    configuration file into a python module where they are published as
    normal fields. The general flow is like this::

      default.cfg ---[SettingsParser]---\
                                         [SettingsManager]--->[settings]
                                               /
                user.cfg---[SettingsParser]---/

    See :class:`P.SettingsParser`
    See :class:`B.settings`

    The default configurations should be located in:

    * C{biskit/Biskit/data/defaults/settings.cfg}      --> :class:`B.settings`
    * C{biskit/Biskit/data/defaults/settings_Mod.cfg}  --> :class:`B.Mod.settings`
    * C{biskit/Biskit/data/defaults/settings_Dock.cfg} --> :class:`B.Dock.settings`

    The user configurations are expected in files of the same name in
    C{~/.biskit/}.
    """

    USER_HEADER = """
##     This is a Biskit user configuration file. The parameters in
##     this file are overriding the default parameters given in
##     %(fdefault)s.
##     If missing, Biskit creates a new user configuration file with
##     those parameters for which the default value seems
##     invalid. The remaining parameters are commented out.

##     Parameters in this file will be accessible from within python as
##     fields of Biskit.settings. For example::
##
##       leaprc = some/path/to/leaprc  # some comment
##
##     will lead to a variable in Biskit.settings::
##    
##     >>> import Biskit.setting as S
##     >>> S.leaprc
##     >>> 'some/path/to/leaprc'

##     ...If, and only if, leaprc also exists in the default settings
##     file.  Parameters that are not listed in the default settings file
##     are ignored.

##     The default type of parameters is str. A prefix to the name like
##     'int-', 'float-', 'bool-', etc. will be interpreted as
##     type-casting. For example::
##
##       float-nice_value = 10  # some comment
##
##     will lead to a variable in Biskit.settings::
##
##     >>> S.nice_value
##     >>> 10.0

"""

    def __init__( self, fdefault, fuser, createmissing=False, verbose=1 ):
        """
        :param fdefault: default configuration file
        :type  fdedault: str
        :param fuser: user configuration file
        :type  fuser: str
        :param createmissing: create user config file if missing
        :type  createmissing: bool
        :param verbose: verbosity level (default: 1)
        :type  verbose: 1|0
        """
        self.verbose = verbose
        self.fdefault = fdefault
        self.fuser = fuser
        self.createmissing = createmissing
        self.fusermissing = not os.path.exists( T.absfile(fuser) )

        self.settings = []  #: will hold extracted Setting's 

    def __update( self, cfg_default, cfg_user ):
        """
        Override default settings by valid (or equally invalid) user settings.

        :param cfg_default: settings read in from default file
        :type  cfg_default: dict {'str':SettingsParser.Setting}
        :param cfg_user   : settings read in from user config file
        :type  cfg_user   : dict {'str':SettingsParser.Setting}

        :return: configuration with valid user settings overriding default ones
        :rtype: dict {'str':SettingsParser.Setting}
        """
        r = {}

        for name, default in cfg_default.items():

            next = cfg_user.get( name, default )

            if next.error and not default.error:

                if self.verbose: B.EHandler.warning(\
                    'User setting %s is reset to default (%r),\n\treason: %s'\
                    % (name, default.value, next.error)\
                    + '\n\tPlease check %s!' % self.fuser )

                next = default

            r[name] = next

        return r


    def collectSettings( self ):
        """
        Parse and combine default and user-defined config files.
        """
        try:
            pdefault = P.SettingsParser( self.fdefault )
            cdefault = pdefault.parse()

            try:
                puser = P.SettingsParser( self.fuser )
                cuser = puser.parse()

            except IOError as e:
                if self.verbose: B.EHandler.warning(
                    'Could not find file with user-defined settings in %s' \
                    % self.fuser, trace=0, error=0)

                cuser = {}

            self.settings = self.__update( cdefault, cuser )

        except P.SettingsError as e:
            B.EHandler.fatal( str(e) )


    def writeUserSettings( self, errorsonly=False ):
        """
        Create a settings file with all options that are invalid with their
        default value.
        """
        try:
            T.backup( self.fuser )  ## create backup if file already exists

            fpath = os.path.dirname(self.fuser)
            if not os.path.exists( fpath ):
                if self.verbose:
                    B.EHandler.warning('Creating folder %s for Biskit settings.'\
                                       %fpath )
                os.mkdir( fpath )

            sections = [P.Setting.NORMAL, P.Setting.PATH, P.Setting.BIN]
            r = {}

            for section in sections:

                r[ section ] = [ s for s in self.settings.values() \
                                 if s.section == section]
                r[ section ].sort()

            f = open( self.fuser, 'w' )

            f.write( SettingsManager.USER_HEADER % self.__dict__ )

            for section in sections:

                f.write( '[%s]\n' % section )
                f.write('\n')

                for param in r[section]:

                    if (not errorsonly) or param.error:
                        f.write( param.formatted() + '\n')
                    else:
                        f.write( '## ' + param.formatted() + '\n') 

                f.write('\n')

            f.close()

        except OSError as e:
            raise WriteCfgError(e)

    def settings2dict( self ):
        """
        Create dictionary from settings.
        :return: dictionary of parameter names (keys) and values
        :rtype: dict {str : any}
        """
        return dict( [ (s.name, s.value) for s in self.settings.values() ] )


    def updateNamespace( self, ns ):
        """
        1. Parse in default configuration and user configuration file
        2. Merge the two, preferring valid user settings
        3. Create missing user configuration file if createmissing=True
        4. Insert parameters into the given namespace

        :param ns: namespace of a module ( obtained with locals() )
        :type  ns: dict {str:any}
        """
        self.collectSettings()

        if self.fusermissing and self.createmissing:
            if self.verbose:
                B.EHandler.warning('Creating new user configuration file %s.' \
                                   % self.fuser, trace=0, error=0)
            self.writeUserSettings( errorsonly=True )

        d = self.settings2dict()

        ns.update( d )


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_SettingsManager(self):
        """SettingsManager test"""

        f_in = T.dataRoot() + '/defaults/settings.cfg'
        self.f_out =  T.tempDir() + '/settings.cfg'

        self.m = SettingsManager( f_in, self.f_out,
                                  createmissing=True,
                                  verbose=self.local )

        ns = locals()             ## fetch local namespace

        self.m.updateNamespace( ns ) ## parse and insert options into namespace

        if self.local:
            globals().update( locals() ) ## publish namespace for debugging

        r = self.m.settings2dict()['testparam']

        self.assertEqual( r, 42) ## from 'int-testparam = 42' in settings.cfg

    def cleanUp(self):
        T.tryRemove( self.f_out )


if __name__ == '__main__':

    BT.localTest()
