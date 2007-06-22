#!/usr/bin/env python

## Setup a host list (read from from Biskit/hosts.py)

## last $Author$
## last $Date$
## $Revision$
  
import ConfigParser
import os, re
import string as S
import types
import sys
import socket

import Biskit.tools as T

import Biskit.hosts as hosts


def _use():
    print """
setup_hosts.py: Setup the host list neded for Biskit distributed calculations.


Usage: Run this script once and let it create the empty host list in
       ~/.biskit/hosts.dat.
       Add your avaliable hosts to the list. There are three different sections
       to which you can choose to add a host:
           - own_hosts:    omputers reserved for own use, highest priority
           - shared_hosts: computers shared with others, medium priority
           - others_hosts: computers mainly used by others, lowest priority
       Add your hosts to the corresponding 'dual' or 'single' cpu option.
       Separat the different hosts with a blank space. If you whish to
       temporarily exclude a host from being used, add it to the 'exclude' option

       Optional settings (noce and ram):
       The nice settings can be changed for a specific computer (default values
       are 0 for 'own' and 'shared' and 17 for 'others'). To add nice value
       add the host(s) and the nice value separated by a colon (:) to the
       'option' 'nice'. Separate multiple hosts with a blank space.
       Example: computer1.local.net:12  computer2:8 computer3.local.net:5
       In the same way avaliablbe RAM in GB can be added. The default values
       here are 0.5 for a single cpu machine and 1.0 GB for a dual cpu machine.

       
Syntax:  setup_hosts  -i |list, input file| -w |str, out file|


Options:
     -i  |filename| read variable names from this file
     -w  |filename| write variables to this file, if the file already exists
           it will be updated,
     -d  |Yes| accept default values for -i and -w
     
Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)

## hosts file
hosts_file = os.path.expanduser('~/.biskit/hosts.dat')

def defOptions():
    bispath = os.getenv('BISKIT') +'/Biskit/'
    return {'i':bispath+'hosts.py',
            'w':hosts_file}


def parseHosts( file, pattern, config ):
    """
    Parse hosts.py and collect list of all host options
    -> ComfigParser object
    """
    settings = open( file, "r")

    p = re.compile( pattern )

    for l in settings.readlines():
        
        if not l[0] == '#':
            m = p.search(l)
            if type(m) != types.NoneType:
                s = m.group('section')
                o = m.group('option')
                c = m.group('comment')
                if not config.has_section( s ):
                    config.add_section( s )
                if not o in config.options( s ):
                    config.set( s, o, c )
        
    settings.close()
    return config


def getHosts( conf, section, option ):
    """
    Strip comment from host list (separator: #)
    """
    setting = conf.get( section, option )
    setting = S.split( setting, '#' )
    return S.strip( setting[0] )


def toList( s ):
    """
    Convert space separated string to list. If the substrings
    contain a colon strip everything off after the colon.
    """
    lst = [ i for i in S.split(s)]
    result = []
    for l in lst:
        d = S.split(l,':')
        if len(d) >= 2:
            result += [d[0]]
        else:
            result+= [l]     
    return result

    
def checkHost( conf ):
    """
    Check that the hosts are valid.
    """
    report = {}
    for s in conf.sections():
        option = {}
        for o in conf.options(s):
            host_fail = []
            host_pass = []
            var = toList( getHosts( conf, s, o ) )
            for v in var:
                try:
                    socket.gethostbyaddr(v)
                    host_pass += [v]
                except:
                    host_fail += [v]
            option[o] = [ host_fail, host_pass ]
        report[s] = option
    return report


def mergeConfig( conf_1, conf_2 ):
    """
    merge settings betweem conf_1 and conf_2 by adding settings in
    conf_2 that is not already in conf_1
    """
    for s in conf_2.sections():
        ## check that conf_1 has section 
        if not conf_1.has_section( s ):
            conf_1.add_section( s )
 
        for o in conf_2.options( s ):
            ## update conf_1
            if not conf_1.has_option( s, o ):
                conf_1.set( s, o, conf_2.get(s,o) )
            ## append comment from conf_2 to conf_1 (if not identical)
            else:
                com_2 = S.strip( S.split( conf_2.get(s,o), '#' )[-1] )
                com_1 = S.strip( S.split( conf_1.get(s,o), '#' )[-1] )
                if not com_2 == com_1 and len(com_2)>0:
                    com = '%s # %s'%( conf_1.get(s,o), com_2 )
                    conf_1.set( s, o, com )
            
    return conf_1


def writeConfig( conf, file ):
    """
    write config file
    """
    if not file:
        conf.write(sys.stdout)
    else:
        if not os.path.exists( os.path.dirname(file) ):
            os.mkdir( os.path.dirname(file) )
        config_file = open( file, 'w' )
        conf.write( config_file )
        config_file.close()


def test():
    options = defOptions()
    options['i'] = '/home/Bis/johan/biskit/Biskit/hosts.py'
    options['w'] = '/home/Bis/johan/.biskit/hosts.dat'
    return options

        
########
## Main
if __name__ == '__main__':
    options = T.cmdDict( defOptions() )
    try:
        if options['d'][0] == 'Y':
            pass
    except:
        _use()
        sys.exit()

#    options = test()

in_file = options['i']
out_file= options['w']

##########################################################
## Setup variables

## initiate configParser
conf_settings = ConfigParser.ConfigParser()
conf_read = ConfigParser.ConfigParser()


## Get Biskit settings (bin, paths and databases)
varExp = "[-\w]*\s*=\s*get[HostsDict]{4,5}" # matches: ssh_bin = getSetting
sectExp = "[('\s]*(?P<section>[-\w]*)"      #          ( 'biskit_bin
optExp = "[,'\s]*(?P<option>[-\w]*)[)"      #          ' , 'ssh_bin')
comExp = "[)'\s]*(?P<comment>.*)"           #          # Path to ssh executable

pattern = varExp + sectExp + optExp + comExp

conf_settings = parseHosts( in_file, pattern, conf_settings )

## check ~/.biskit/hosts.dat for settings that are already there
conf_read.read( hosts_file )

## merge settings in  ~/.biskit/hosts.dat with all needed settings 
conf_read = mergeConfig( conf_read, conf_settings )

## check that variables in ~/.biskit/hosts.dat are valid
report = checkHost( conf_read )

## write configuration file
writeConfig( conf_read, out_file )

## print host list status
print '\nSECTION \tOPTION     \tHOST INFORMATION'
print '\t\t\t\tPASSED \tFAILED'
for s in report.keys():
    print '%s'%s
    for o in report[s]:
        fail = report[s][o][0]
        if len(fail) == 0:
            fail = 'None'
        print '\t\t%s     \t%i \t%s'%(o, len(report[s][o][1]), fail)

## Print message
print '\n%s \nHost list have been written to %s. \nPlease edit this file according \
to your environment and \nthen run this script again to verify the list.\n%s'\
%('='*60, out_file, '='*60)
