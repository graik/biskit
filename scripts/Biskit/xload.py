#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

import os, commands
import time
import string
from Biskit.tools import *

def _use():
    print """xload: show osxview for several machines.\n
    Syntax:     xload |machine_file|
    machine_file ... text file with one machine name per line\n"""

def xload(machines):
    """
    Take list of machine names, run xosview on them,
    -> lits of pids
    """
    pids = []
    for machine in machines:
        try:
            print "launching on ", machine,
            cmd = "ssh %s 'xosview &' &" % string.strip(machine)
            r = os.system( cmd)
            time.sleep(0.6)
            if r == 0:
                p="""ps -ef | grep "ssh %s xosview" | grep -v grep >> pids.dat"""\
                         % string.strip(machine)
                os.system(p)

        except Exception:
            print "xload: ", lastError()

    ## extract pids from temp.dat
    try:
        f = open("pids.dat")
        lines = f.readlines()
        for line in lines:
            pid = string.split(line)[1]
            pids += [pid]
        f.close()
        os.remove("pids.dat")
    except Exception, why:
        print "xload Error: ", lastError()

    return pids

def xload_fromFile( fname ):
    """
    xload for all machines in textfile.
    -> list of pids
    """
    f = open( fname )
    machines = f.readlines()
    f.close()
    return xload( machines )
        
def kill_allFromFile( machines ):
    """
    Kill all process ids with -9
    incomplete
    """
    for id in pids:
        try:
            os.kill( id, 9 )
        except Exception, why:
            print "Error killing ",id, ": ", lastError()

##################################################
# main function
##################################################
if __name__ == '__main__':
    if len (sys.argv) < 2:                      # at least 1 argument
        _use()
        sys.exit(0)
    else:
        s = xload_fromFile( sys.argv[1] )
        for pid in s:
            print pid,
        print
