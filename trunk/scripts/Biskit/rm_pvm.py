#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## delete all /tmp/pvm* files on all hosts

import os
import sys, string
from Biskit.hosts import nodes_all

def cleanPvm( host ):
    cmd = "ssh %s 'rm /tmp/pvm*'"
##     cmd = "ssh %s 'rm -f /work/@*'"
    os.system( cmd % host )

def df( host ):
    cmd = "ssh %s 'cd /tmp; df -h .'"
    os.system( cmd % host )

s = "removing local /tmp/pvm* on ALL hosts! Continue? (Y/N)+ENTER: "
if string.upper( raw_input(s) ) != 'Y':
    print "Cancelled"
    for h in nodes_all:
        try:
            print h
            df( h )
        except:
            pass
    sys.exit(0)

for h in nodes_all:

    try:
        print h
        cleanPvm( h )
    except Exception, why:
        pass
##         print "Error cleaning %s: ", why
