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
Parse a Biskit settings file.
"""

import os
import user
import sys
import ConfigParser

import Biskit as B
import Biskit.tools as T

class SettingsError( Exception ):
    pass

class InvalidType( SettingsError ):
    pass

class InvalidValue( SettingsError ):
    pass

class SettingsWarning( SettingsError ):
    pass

class InvalidPath( SettingsWarning ):
    pass

class InvalidBinary( SettingsWarning ):
    pass


class SettingsParser:
    """
    A config file parser on steroids -- performs the following tasks:

    1. read a ini-style settings file
    2. type-cast options (e.g. of the form int-some_name into int(some_name) )
    3. validate that all entries of section [PATH] point to existing paths
    4. absolutize all valid paths
    5. validate that all entries of section [BINARIES] point to binaries
    """

    NOVALIDATION = 0
    VALID_PATH = 1
    VALID_BIN = 2

    def __init__(self, ini):

        self.f_ini = ini
        self.result = {}


    def __validPath( self, v ):
        """
        @param v: potential path name
        @type  v: str

        @return: validated absolute Path
        @rtype : str

        @raise InvalidPath: if path is not found
        """
        try:
            v = T.absfile( v )
            if not os.path.exists( v ):
                raise InvalidPath, 'invalid path %r' % v

            return v

        except InvalidPath, e:
            raise
        except Exception, e:
            raise InvalidPath, 'error during path validation: %r' % str(e)
        

    def __validBinary( self, v ):
        """
        @param v: potential binary path
        @type  v: str

        @return: validated absolute path to existing binary
        @rtype : str

        @raise InvalidBinary: if path is not found
        """
        try:
            return T.absbinary( v )
        except IOError, msg:
            raise InvalidBinary( str(msg) )


    def __type( self, option, default=str ):
        """
        Extract type from option name.
        
        @param option: name of parameter
        @type  option: str
        @param default: default type [str]
        @type  default: type

        @return: type, stripped option name (e.g. 'int_var1' -> int, 'var1')
        @rtype: type, str

        @raise TypeError, if type cannot be interpreted
        """
        t = default
        o = option

        if option.count('-') > 0:

            try:
                
                splt = option.split('-')

                s = splt[0]
                o = ''.join( splt[1:] )
                
                t = eval( s )

                if not type(t) is type:
                    raise TypeError, '%s is not a valid type' % s
                
            except Exception, e:
                raise TypeError, 'Cannot extract type from %s: %r'\
                      % option, e

        return t, o


    def __process( self, option, value, valid=SettingsParser.NOVALIDATION ):
        """
        @param option: option name
        @type  option: str

        @param value: option value
        @type  value: str

        @param path: validate as path name (warning if not existent)
        @type  path: bool
        
        @return: processed option name, processed value
        @rtype: str, any

        @raise SettingsError: InvalidPath/Binary/Type/Value
        """
        try:

            t, name = self.__type( option )  ## target type and target name

            v = value.split('#')[0].strip()  ## strip off comments

            if not v:                        ## don't return empty strings
                v = None

            if valid is SettingsParser.VALID_PATH:
                v = self.__validPath( v )

            if valid is SettingsParser.VALID_BIN:
                v = self.__validBinary( v )

            try:
                v = t( v )
            except ValueError, e:
                raise InvalidValue, '%s: cannot convert "%s" to %r.' %\
                  (option,value,t)

        except SettingsError, e:
            self.error[ option ] = str(e)
            raise e

        return name, v

    def __processSection( self, items, valid=SettingsParser.NOVALIDATION ):
        """
        @param items: section comming from ConfigParser
        @type  items: [ ( str, str ) ]

        @param path: validate path names
        @type  path: bool

        @return: validated path values
        @rtype : dict
        """
        r = {}

        for name, value in items:

            try:

                n, v = self.__process( name, value, valid )

                r[ n ] = v

            except SettingsWarning, e:
                B.EHandler.warning( str(e), trace=0, error=0 )

        return r


    def parse( self ):
        """
        @return dict of type-cast params contained in fini
        @rtype: dict of dict
        """
        try:

            ## read from file
            c = ConfigParser.ConfigParser()
            c.read( self.f_ini  )

            r = {}
            self.error = {}  ## options with errors

            for section in c.sections():

                valid = SettingsParser.NOVALIDATION or \
                        SettingsParser.VALID_PATH * (section == 'PATHS') or \
                        SettingsParser.VALID_BIN  * (section == 'BINARIES')

                r.update( self.__processSection(c.items(section), valid) )

            return r

        except Exception, e:
            B.EHandler.fatal('Biskit could not read its settings file %s'\
                           % self.f_ini )
            

#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self ):
        """
        run function test

        @return: 1
        @rtype: int
        """
        p = SettingsParser( T.projectRoot()+'/external/defaults/settings.cfg')

        globals().update( locals() )

        r = p.parse()

        globals().update( locals() )

        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    
        

if __name__ == '__main__':

    test = Test()

    assert test.run( ) == test.expected_result()
