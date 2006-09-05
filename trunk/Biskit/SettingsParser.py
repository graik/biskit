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

class InvalidPath( SettingsError):
    pass

class SettingsParser:
    """
    A config file parser on steroids -- performs the following tasks:

    1. read a ini-style settings file
    2. type-cast options (e.g. of the form int-some_name into int(some_name) )
    3. validate that all entries of section [PATH] point to existing paths
    4. absolutize all valid paths
    5. validate that all entries of section [EXE] point to existing binaries
    """

    def __init__(self, ini):

        self.f_ini = ini
        self.result = {}


    def __validPath( self, v ):
        """
        @param v: potential path name
        @type  v: str

        @return: validated Path
        @rtype : str

        @raise InvalidPath: if path is not found
        """
        try:
            v = T.absfile( v )
            if not os.path.exists( v ):
                raise InvalidPath, 'invalid path' % v

            return v

        except Exception:
            raise InvalidPath, 'error during path validation: %r'
        

    def __validBinary( self, v ):
        try:
            v = T.absbinary( v )
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
                    raise SettingsError, '%s is not a valid type' % s
                
            except:
                EHandler.warning('Ignoring type %s of option %s'% s,option)

        return t, o


    def __process( self, option, value, path=False ):
        """
        @param option: option name
        @type  option: str

        @param value: option value
        @type  value: str

        @param path: validate as path name (warning if not existent)
        @type  path: bool
        
        @return: processed option name, processed value
        @rtype: str, any

        @raise ValueError: if target type and value are incompatible
        """

        t, name = self.__type( option )  ## target type and target name

        v = value.split('#')[0].strip()  ## strip off comments

        if not v:   ## don't return empty strings
            v = None

        try:

            if path:
                v = self.__validPath( v )

            v = t( v )

        except InvalidPath, e:
            self.error[ name ] = str(e)

        except ValueError:
            msg =  '%s: cannot convert value %s to %r.\n' %\
                  (option,value,t)
            self.error[ option ] = msg
            msg += 'Please verify %s.' % self.f_ini
            raise ValueError( msg )

        return name, v

    def __processSection( self, items, path=False ):
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

                n, v = self.__process( name, value, path )

                r[ n ] = v

            except ValueError, e:
                B.EHandler.warning( str(e), trace=0, error=0 )

        return r


    def parse( self ):
        """
        @return dict of type-cast params contained in fini
        @rtype: dict of dict
        """
        ## read from file
        c = ConfigParser.ConfigParser()
        c.read( self.f_ini  )

        r = {}
        self.error = {}  ## options with errors

        for section in c.sections():

            paths = (section == 'PATHS')

            r.update( self.__processSection(c.items(section), paths) )

        return r


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
        p = SettingsParser( T.projectRoot()+'/external/defaults/settings.dat')
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
