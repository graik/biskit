#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Date$
## last $Author$

from Biskit.tools import *
from Biskit.hosts import hosts_all, nice_dic
from PVM.TrackingJobMaster import *

def _use( options ):
    print """
Syntax: a_trajQuality -i |traj_1 traj_2 .. traj_n| [-a -h |n_hosts| -w]
    pvm must be running on the local machine!

Result: eps with quality plots in folder of traj files
        
Options:
    -h    number of hosts to be used
    -a    first add hosts to pvm
    -w    display a xterm window for each node

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


class QualMaster(TrackingJobMaster):

    def __init__(self, trajFiles, n_hosts=20, **kw):
        """
        dat - data dictionary
        hosts - list of host names
        """
        dat = {}
        i = 0
        for f in trajFiles:
            dat[i] = absfile( f )
            i += 1

        niceness = nice_dic
        hosts = cpus_all[ :n_hosts ]

        project_path = projectRoot() + '/scripts/modules/'
        slave_script = project_path + 'QualSlave.py'

        TrackingJobMaster.__init__(self, dat, 1, hosts, niceness,
                           slave_script, **kw)


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.
        """
        return 1

    def done(self):
        self.exit()
        

if __name__ == '__main__':

    options = cmdDict( {'h':20} )

    if len( sys.argv ) < 2:
        _use( options )

    ## extract file names
    files = toList( options['i'] )

    master = QualMaster( files, int(options['h']),
                         show_output=options.has_key('w'),
                         add_hosts=options.has_key('a') )

    master.start()
