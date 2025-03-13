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
Parse a Biskit settings file.
"""

import os
import configparser

import biskit as B
import biskit.tools as T

class SettingsError( Exception ):
    pass

class InvalidType( SettingsError ):
    pass

class InvalidValue( SettingsError ):
    pass

class InvalidFile( SettingsError ):
    pass

class SettingsWarning( SettingsError ):
    pass

class InvalidPath( SettingsWarning ):
    pass

class InvalidBinary( SettingsWarning ):
    pass


class CaseSensitiveConfigParser( configparser.ConfigParser ):
    """
    Change ConfigParser so that it doesn't convert option names to lower case.
    """
    def optionxform(self, optionstr):
        return optionstr

class Setting:
    """
    Simple container for a single parameter
    """
    ## types of settings // section name in cfg file
    NORMAL = 'NORMAL'
    PATH = 'PATHS'
    BIN = 'BINARIES'

    def __init__( self, name=None, value=None, vtype=str, comment=None,
                  error=None, section=NORMAL ):
        self.name = name
        self.value = value
        self.vtype = vtype
        self.comment = comment
        self.error = error
        self.section = section

    def typeCast( self, vtype ):
        """
        Recast value to a new type. Value None remains unchanged.
        :param vtype: new type for value
        :type  vtype: type
        :raise InvalidValue: if current value is incompatible with vtype
        """
        try:
            if not self.value is None:
                self.value = vtype( self.value )
            self.vtype = vtype
        except ValueError as e:
            raise InvalidValue('%s: cannot convert "%s" to %r.' %\
                  (self.name,self.value,vtype))

    def __repr__( self, tab='' ):
        error = ''
        if self.error: error = '(!)'

        comment = self.comment or ''
        if comment: comment = ' # '+comment

        return 'SettingsParser.Setting: %s%s = %s%s(%s)%s%s' %\
               (error, self.name, tab, self.vtype.__name__, str(self.value),\
                tab, comment)

    def __str__( self ):
        return self.__repr__( tab='\t' )

    def __lt__(self, other):
        """Compare settings instances by their name"""
        if isinstance( other, self.__class__ ):
            return self.name < other.name
        return self < other
    
    def __eq__( self, other ):
        """Compare settings instances by their name"""
        if isinstance( other, self.__class__ ):
            return self.name == other.name
        return self == other


    def formatted( self ):
        """
        :return: parameter formatted for setting file
        :rtype: str
        """
        comment = ''
        error = ''
        name = self.name
        value = self.value or ''

        if self.vtype != str:
            name = self.vtype.__name__ + '-' + name

        if self.comment:
            comment = '\t## ' + self.comment

        if self.error:
            error = '\t#! ' + self.error + ' !'

        return '%s = %s%s%s' % (name, str(value), comment, error)


class SettingsParser(object):
    """
    A config file parser on steroids -- performs the following tasks:

      1. read a ini-style settings file
      2. type-cast options (e.g. of the form int-some_name into int(some_name))
      3. validate that all entries of section [PATHS] point to existing paths
      4. absolutize all valid paths
      5. validate that all entries of section [BINARIES] point to binaries
    """

    def __init__(self, ini):

        self.f_ini = T.absfile( ini )
        self.result = {}


    def __validPath( self, v ):
        """
        :param v: potential path name
        :type  v: str

        :return: validated absolute Path
        :rtype : str

        :raise InvalidPath: if path is not found
        """
        try:
            v = T.absfile( v )
            if not v or not os.path.exists( v ):
                raise InvalidPath('invalid path %r' % v)

            return v

        except InvalidPath as e:
            raise
        except Exception as e:
            raise InvalidPath('error during path validation: %r' % str(e))


    def __validBinary( self, v ):
        """
        :param v: potential binary path
        :type  v: str

        :return: validated absolute path to existing binary
        :rtype : str

        :raise InvalidBinary: if path is not found
        """
        try:
            if not v:
                raise IOError('empty path')

            return T.absbinary( v )

        except IOError as msg:
            raise InvalidBinary( str(msg) )


    def __type( self, option, default=str ):
        """
        Extract type from option name.

        :param option: name of parameter
        :type  option: str
        :param default: default type [str]
        :type  default: type

        :return: type, stripped option name (e.g. 'int_var1' -> int, 'var1')
        :rtype: type, str

        :raise TypeError: if type cannot be interpreted
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
                    raise TypeError('%s is not a valid type' % s)

            except Exception as e:
                raise TypeError('Cannot extract type from %s: %r'\
                      % option).with_traceback(e)

        return t, o


    def __process( self, option, value, section=Setting.NORMAL ):
        """
        :param option: option name
        :type  option: str

        :param value: option value
        :type  value: str

        :param section: which section are we working on
        :type  section: str

        :return: new setting
        :rtype: Setting

        :raise SettingsError: InvalidType or Value
        """
        r = Setting( section=section )

        try:

            x = value.split('#')             ## split off comments
            r.value = x[0].strip() or None   ## don't return empty strings

            if len(x) > 1:
                r.comment = ''.join( x[1:] )

            vtype, r.name = self.__type( option )
            r.typeCast( vtype )

            if section == Setting.PATH:
                r.value = self.__validPath( r.value )

            if section == Setting.BIN:
                r.value = self.__validBinary( r.value )

        except SettingsWarning as e:           ## catch and record warnings
            r.error = str(e)

        return r


    def __processSection( self, items, section=Setting.NORMAL, verbose=False ):
        """
        :param items: section comming from ConfigParser
        :type  items: [ ( str, str ) ]

        :param section: which config section are we working on?
        :type  section: str

        :return: validated path values
        :rtype : dict, {str: Setting}
        """
        r = {}

        for name, value in items:

            s = self.__process( name, value, section )

            r[ s.name ] = s

            if verbose and s.error:
                B.EHandler.warning( s.error, trace=0, error=0 )

        return r


    def parse( self ):
        """
        :return: dict of type-cast params contained in fini
        :rtype: dict, {str: Setting}

        :raise IOError: if the settings file does not exist
        :raise SettingsError: (InvalidFile, InvalidValue, InvalidType)
        """
        try:
            ## read from file
            c = CaseSensitiveConfigParser()

            if c.read( self.f_ini ) != [ self.f_ini ]:
                raise IOError('Settings file %s not found.' % self.f_ini)

            for section in c.sections():

                self.result.update(
                    self.__processSection( c.items(section), section) )

        except configparser.Error as e:
            raise InvalidFile('Error parsing settings file %s: ' %\
                  self.f_ini + str(e))

        return self.result

    def __repr__( self):
        r = super( SettingsParser, self).__repr__()
        err = len( [ s for s in self.result.values() if s.error ] )
        r += ' -- %i entries, (!) %i errors' % (len( self.result ), err)

        values = list(self.result.values())
        values.sort()

        for v in values:
            r+= T.clipStr( '\n - %s' % str(v), 75)

        return r

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_SettingsParser(self):
        """SettingsManager test"""
        p = SettingsParser( T.dataRoot() + '/defaults/settings.cfg')

        p.parse()

        t = p.result.get('testparam', Setting())

        self.assertEqual( (t.name, t.value), ('testparam', 42) )

        return t.name, t.value


if __name__ == '__main__':

    BT.localTest()
