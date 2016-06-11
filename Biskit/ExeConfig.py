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


"""
Collect settings for an external program from a configuration file.
"""

import ConfigParser
import user, os.path

from Biskit.Errors import BiskitError
from Biskit import EHandler

import Biskit.tools as T


class ExeConfigError( BiskitError ):
    pass

class CaseSensitiveConfigParser( ConfigParser.SafeConfigParser ):
    """
    Change ConfigParser so that it doesn't convert option names to lower case.
    """
    def optionxform(self, optionstr):
        return optionstr


class ExeConfig( object ):

    """
    ExeConfig
    =========
    
    Manage the settings that Executor needs for calling an external
    program.

    ExeConfig is initialised with a program name and expects to find
    the configuration for this program in
    C{ ~/.biskit/exe_|name|.dat }. Only if nothing is found there, it looks
    for C{ Biskit/data/defaults/exe_|name|.dat } in the biskit installation
    folder.  If neither of the two files are found, an error is raised
    (strict=1) or the binary is assumed to be |name| and must be
    accessible in the search path (strict=0).


    Example
    -------
    
    The configuration file (exe_name.dat) should look like this::

      ---- start example configuration file ----
      [BINARY]

      comment=the emacs editor
      bin=emacs
      cwd=
      shell=0
      shellexe=
      pipes=0
      ## Use new environment containing only variables given below
      replaceEnv=0

      [ENVIRONMENT]

      HOME=
      EMACS_CONFIG=~/.emacs/config.dat
      ---- end of example file ----


    This example config would ask Executor to look for an executable
    called 'emacs' in the local search path. Before running it,
    executor should check that a variable $HOME exists in the local
    shell environment (or raise an error otherwise) and set the
    variable $EMACS_CONFIG to the given path.

    The other settings specify how the program call is done (see also
    Python 2.4 subprocess.Popen() ):

      - cwd   ...    working directory (empty -- current working directory)
      - shell ...    wrap process in separate shell
      - shellexe ... which shell (empty -- sh)
      - pipes ...    paste input via STDIN, collect output at STDOUT

    Missing options are reset to their default value; See
    L{ ExeConfig.reset() }.  All entries in section BINARY are put into the
    name space of the ExeConfig object. That means an ExeConfig object x
    created from the above file can be used as follows:

      >>> x = ExeConfig( 'emacs' )
      >>> x.cwd is None
      >>> True
      >>> print x.comment
      >>> the emacs editor
    """

    ## static fields
    PATH_CONF   = user.home + '/.biskit'
    PATH_CONF_DEFAULT = os.path.join( T.dataRoot(), 'defaults' )
    SECTION_BIN = 'BINARY'
    SECTION_ENV = 'ENVIRONMENT'

    def __init__( self, name, strict=1 ):
        """
        @param name: unique name of the program
        @type  name: str
        @param strict: insist on a config file exe_name.dat
                       and do not tolerate missing environment variables
                       (default: 1)
        @type  strict: 0|1
        
        @raise ExeConfigError: if strict==1 and config file incomplete/missing
        """
        self.name = name    #: identifier
        #: path to configuration file
        self.dat  = os.path.join( self.PATH_CONF, 'exe_%s.dat' % name )

        if not os.path.exists( self.dat ):
            self.dat = os.path.join( self.PATH_CONF_DEFAULT,'exe_%s.dat'%name )

        #: True if a configuration file was found
        self.dat_found =  os.path.exists( self.dat )

        self.strict = strict

        self.env_checked = 0 ## environment was verified

        if strict and not self.dat_found:

            raise ExeConfigError,\
                  'Could not find configuration file %s for program %s.'\
                  % (self.dat, self.name)

        self.conf = CaseSensitiveConfigParser()
        self.conf.read( self.dat )

        self.reset()
        self.update()


    def reset( self ):
        """
        Reset all required parameters. Called at creation
        """
        ## default values
        self.comment = 'no comment or missing configuration file'
        self.bin = self.name
        self.shell = 0
        self.shellexe = None
        self.pipes = 0
        self.cwd = None   #'./'

        self.replaceEnv = 0
        self.env = None


    def update( self ):
        """
        Load settings from associated configuration file (if available).
        Is automatically called at creation.
        
        @raise ExeConfigError: if section [BINARY] was not found in the file
        """
        ## get parameters from config file if available; type-cast values 
        try:
            dconf = self.conf.items( self.SECTION_BIN )

            for key, value in dconf:

                ## default type is string
                t = type( self.__dict__.get( key, '' ) )
                if t is type(None):
                    t = str

                ## leave default value if None is given
                if value is not '':
                    self.__dict__[ key ] = t( value )

        except ConfigParser.NoSectionError:
            if self.strict:
                raise ExeConfigError,\
                      'Could not find BINARY section in %s.' % self.dat

        try:
            self.env = dict( self.conf.items( self.SECTION_ENV ) )
        except:
            pass


    def validate( self ):
        """
        Validate the path to the binary.
        
        @raise ExeConfigError: if environment is not fit for running
                               the program
        """
        try:
            self.bin = T.absbinary( self.bin ) ## raises IOError if not found

            missing = self.update_environment()
            report = '%s is missing environment variables: %r'\
                     % (self.name, missing )

            if missing and self.strict:
                raise ExeConfigError, report

            if missing:
                EHandler.warning( report )

        except IOError, why:
            raise ExeConfigError, str(why) + ' Check %s!' % self.dat


    def environment( self ):
        """
        Get needed environment variables.
        
        @return: dictionary with environment for subprocess.Popen;
                 empty, if no environment was specified
        @rtype: {str:str} OR None

        @raise ExeConfigError: if env was not yet checked by update_environment
        """
        if not self.env_checked:
            raise ExeConfigError, 'Environment not yet checked, validate()!'

        return self.env


    def update_environment( self ):
        """
        Check for missing environment settings.
        
        @return: names of required but missing environment variables
        @rtype: [str]
        """
        missing = []

        if self.env:

            for key,value in self.env.items():

                if value is '':

                    if os.getenv( key ) is None:
                        missing += [ key ]
                    else:
                        self.env[ key ] = os.getenv( key )

        self.env_checked = 1

        return missing


    def __repr__( self ):
        s = 'ExeConfig for %s' % self.name
        for k,v in self.__dict__.items():
            s += '\n%10s \t%r' % (k,v)
        return s

    def __str__( self ):
        return self.__repr__()



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """ExeConfig test"""
    
    def test_ExeConfig( self ):
        """ExeConfig test (validate xclock)"""

        x = ExeConfig( 'ls', strict=1 )
        x.validate()

        if self.local:
            print x.bin

        self.assertEquals( True, 'ls' in x.bin )
    
        
if __name__ == '__main__':

    BT.localTest()
