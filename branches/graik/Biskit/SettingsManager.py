##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$

import os
import user
import sys
import ConfigParser

import Biskit
import Biskit.tools as T

class SettingsError( Exception ):
    pass

class SettingsManager:
    """
    Merge and export information from settings_default and
    .biskit/settings.dat into a module settings. SettingsManager is
    automatically called when importing Biskit.

    See Biskit.settings for details.
    """

    ## Static fields
    conf_file = user.home + '/.biskit/settings.dat'

    ## report missing settings at the end
    warning_0="Biskit is missing %i parameters which were set to defaults:\n\n"
    warning_1= "%s/ %s: \t%s\n"
    warning_2= "\nPlease create or verify %s with setup_env.py!\n"

    def __init__( self, defaults_module ):
        """
        defaults_module - module, source of parameter names and default values
        """
        self.defaults_module = defaults_module
        
        ## Read configuration file
        self.conf     = ConfigParser.ConfigParser()
        self.conf.read( SettingsManager.conf_file )

        self.missing_settings = []


    def getFromConfig( self, section, option, value_type=str ):
        """
        section, option - str
        value_type      - type, try to convert value to this type [str]
        -> any, the value extracted from the config file
        !! ValueError, if type conversion fails
        !! SettingsError, if there is no parameter of that section / name
        """
        if not self.conf.has_option( section, option):
            raise SettingsError, 'cannot find %s / %s in %s' % \
                  (section, option, self.conf_file )
        
        setting = self.conf.get( section, option )

        r = setting.split('#')[0].strip()  ## strip off comments
        if r:                      ## don't return empty strings

            return value_type( r )
       

    def __getSetting( self, section, option, default ):
        """
        Try to fetch a setting value from file and type cast it.
        section, option - str, name of section and parameter in settings.dat
        default - any, default value if parameter is not found in settings.dat
        -> any, value taken from settings.dat or default, same type as default
        """
        try:

            return self.getFromConfig( section, option, type( default ) )
        
        except ValueError:
            msg =  '%s: cannot convert value %s to %r.\n' %\
                  (option,r,type(default))
            msg += 'Please verify %s.' % self.conf_file
            Biskit.EHandler.warning( msg, trace=0, error=0 )

        except SettingsError:
            pass

        ## didn't find parameter in .biskit/settings.dat
        self.missing_settings += [ (section, option, default) ]

        return default


    def getDefaults( self, module=None ):
        """
        Parse default values defined in module into a dictionary.
        module - python module, Biskit.(Mod/Dock.)settings_default
        -> { str_section: {str_param: (value, str_comment) } }
        """
        module = module or self.defaults_module
        result = {}

        for param, descr in module.__dict__.items():

            if param[0] != '_':   ## ignore private fields of settings_default

                default = descr[0]
                section = descr[1]
                comment = None
                if len( descr ) > 2:
                    comment = descr[2]

                if not section in result:
                    result[ section ] = {}

                result[ section ][ param ] = ( default, comment )

        return result


    def updateNamespace( self, namespace ):
        """
        Parse settings_default.py (and .biskit/settings.dat) and create one
        entry in namespace for each of its fields.
        namespace - dict, name space of module that wants to export parameters
        """
        defaults = self.getDefaults()
            
        for section, sect_dict in defaults.items():

            for param, descr_tuple in sect_dict.items():

                default = descr_tuple[0]
                value   = self.__getSetting( section, param, default ) 

                namespace[ param ] = value  ## add param to module namespace

        self.reportMissing()


    def reportMissing( self ):
        """Report missing parameters (which were recorded by __getSetting )."""
        if len( self.missing_settings ) > 0:

            msg = self.warning_0 % len( self.missing_settings )

            for m in self.missing_settings:
                msg += self.warning_1 % m

            msg += self.warning_2 % SettingsManager.conf_file

            Biskit.EHandler.warning( msg )

            self.missing_settings = []


    def __backup( self, fname ):
        """
        Create backup of file if it already exists.
        """
        if os.path.exists( fname ):
            os.rename( fname, fname + '~' )
            

    def writeConfig( self, settings, fname ):
        """
        Create a config file from a dictionary. This method is used by
        setup_biskit.py.
        settings - { str_section: {str_param: (value, str_comment)}}
        fname    - str, target file name (overwritten)
        !! IOError, if the file cannot be written
        """
        try:

            if not os.path.exists( os.path.dirname(fname) ):
                os.mkdir( os.path.dirname(fname) )

            self.__backup( fname )

            f = open( T.absfile( fname ), 'w' )

            section_names = settings.keys()
            section_names.sort()

            ## create sections in alphabetical order
            for sec_name in section_names:

                f.write( '[%s]\n' % sec_name)

                params = settings[ sec_name ].keys()
                params.sort()

                ## create entries in alphabetical order
                for param in params:

                    descr = settings[ sec_name ][ param ]

                    value   = descr[0]
                    comment = descr[1]

                    f.write( '## %s\n' % comment or 'no comment')

                    f.write( '%s = %s\n\n' % (param, value ) )

                f.write('\n')
        finally:
            try:
                f.close()
            except:
                pass


### TEST ###

if __name__ == '__main__':

    ## importing Biskit calls SettingsManager automatically
    import Biskit.settings_default as D

    m = SettingsManager( defaults_module=D )

    defaults = m.getDefaults( D )

    m.writeConfig( defaults, T.absfile('testconfig.dat') )


