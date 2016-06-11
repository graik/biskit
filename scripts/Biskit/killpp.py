#!/usr/bin/env python
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
##
## kill all processes containing a given phrase on a list of hosts

import os, sys
import string
import commands
import getpass

## modules = '~/biskit'
## modules = os.path.expanduser(modules)
## sys.path.insert(0, modules)

from Biskit.PVM.hosts import nodes_all
from Biskit.tools import *
import Biskit.mathUtils as MU

def pp( host, name, user=None ):
    cmd = "ssh %s '/usr/bin/env ps -efl | grep %s | grep -v grep'"  % (host, name)
    l = commands.getoutput( cmd )

    if l == '':
        return []
    
    l = l.split('\n')

    l = [ p for p in l if p.find('killpp') == -1 ]

    if user:
        l = [ p for p in l if p.split()[2]==user  ]

    return l



def kill( host, name, user=None, ask=1 ):

    l = pp( host, name, user )

    if l:

        print "\n%s:" % host
        for p in l:
            print p

        msg = "KILL Y/N "

        if (not ask) or force or string.upper( raw_input( msg ) ) == 'Y':

            for p in l:
                pid = p.split()[3]
                cmd = "ssh %s 'kill -9 %s'" % (host, pid)
                r = os.system( cmd )
                if r != 0:
                    print "exit code ",r

############
### MAIN ###

if len( sys.argv ) < 2:
    print "kill all processes containing a given phrase on all hosts"
    print "Syntax: killpp.py -n |part_of_name| [-a -f "
    print "                  -h |host1 host2 host3 ..| -e |host4 host5 ..|]"
    print "        -n ... search phrase"
    print "        -a ... kill all without asking"
    print "        -f ... don't ask anything at all"
    print "        -h ... only look on these hosts  (default: all)"
    print "        -e ... exclude one or more hosts from killing"
    sys.exit(0)

o = cmdDict( { 'h':nodes_all } )
phrase = o['n']
hosts  = toList( o['h'] )
user = getpass.getuser()
ask  = 'a' not in o
force= 'f' in o
if 'e' in o:
    hosts = MU.difference( hosts, toList( o['e'] ) )

msg = "KILL All jobs matching '%s' on all hosts?? Y/N " % phrase

if (not ask) and (force or string.upper( raw_input( msg ) ) != 'Y'):
    print "Aborted."
    sys.exit()

for h in hosts:
    try:
        kill( h, phrase, user, ask )
    except KeyboardInterrupt:
        print "\nAborted."
        sys.exit()
    except Exception, why:
        print "Error killing %s: " % str(why),  lastError()
