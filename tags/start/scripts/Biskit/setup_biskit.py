#!/usr/bin/env python2.4
## last $Author$
## last $Date$
## $Revision$

import os.path as osp

import Biskit as B
import Biskit.tools as T
from Biskit.SettingsManager import SettingsManager, SettingsError

import Biskit.settings_default as Defaults_Biskit
import Biskit.Mod.settings_default as Defaults_Mod
import Biskit.Dock.settings_default as Defaults_Dock


def validatePath( value, default ):
    """
    value   - str, value taken from settings.dat or default
    default - str, default value
    -> any, valid value or default or ''
    """
    try:

        if osp.exists( value ):
            print ' -- OK.'
            return value

        if osp.exists( default ):
            print ' -- using default.'
            return default

        print ' -- invalid path %s' % value

        return '%s  # invalid path!' % value

    except TypeError:
        print ' -- not checked.'
        return value  ## value is not a string -- don't validate then


def validateBinary( value, default ):
    """
    value   - str, value taken from settings.dat or default
    default - str, default value
    -> any, valid value or default or ''
    """
    try:
        if T.binExists( value ):
            print ' -- OK.'
            return value

        if T.binExists( default ):
            print ' -- using default.'
            return default

        print ' -- invalid binary %s' % value

        return '%s  # invalid binary!' % value

    except TypeError:
        print ' -- not checked.'
        return value  ## value is not a string -- don't validate then


def isPath( value ):
    return type( value ) is str and osp.sep in value and not '//' in value


def validate( param_dic, default_dic ):
    """
    param_dict   - { str : { str : any } }
    default_dict - { str : { str : any } }
    """
    for sec_name, section in default_dic.items():

        for param_name, default in section.items():

            print sec_name, param_name,

            param, comment = param_dic[ sec_name ][ param_name ]
            default = default[0]

            if param_name[:3] == 'bin' or param_name[-3:] == 'bin':

                param_dic[sec_name][param_name]=\
                              validateBinary(param,default), comment

            elif isPath( param ) or isPath( default ):

                param_dic[sec_name][param_name]= validatePath(param,default),\
                                                 comment

            else:
                print ' -- not checked.'

        print

    return param_dic
                                                                   

def load_and_validate_settings():
    """main method"""
    ## prepare validated settings
    result   = {}
    defaults = {}
    
    for module in [Defaults_Biskit, Defaults_Mod, Defaults_Dock]:

        print "Analyzing %s ..." % module.__name__

        m = SettingsManager( module )
        defaults = m.getDefaults()
        
        for sec_name, section in defaults.items():

            result[ sec_name ] = {}

            for param_name, param in section.items():

                default, comment = param
                value = default
                
                try:
                    value = \
                       m.getFromConfig( sec_name, param_name, type( default ))

                except SettingsError, why:   ## value not found in config file
                    B.EHandler.warning( why, error=0, trace=0 )

                except ValueError:           ## error during type conversion
                    B.EHandler.warning( 'cannot convert %s / %s to %r' %\
                                        (sec_name, param_name, type(default)),
                                        error=0, trace=0)

                result[ sec_name ][ param_name ] = value, comment
                    
        ## validate parameters of this module
        print 'Validating paths and binaries ...'

        result = validate( result, defaults )

    ## finally, create a new settings file
    print
    print 'Writing the new settings file to ', m.conf_file
    m.writeConfig( result, m.conf_file )



#########
## MAIN
#########

try:
    load_and_validate_settings()

except Exception, why:
    B.EHandler.error( 'Error importing Biskit.settings.' )

