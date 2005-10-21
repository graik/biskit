#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## last $Date$
## $Revision$

import sys
import Biskit.tools as T
from Biskit.PVM.TrackingJobMaster import restart

def use():
    if len( sys.argv ) < 2:
        print \
"""
Restart a distributed calculation.
Syntax:  restartPVM.py -i |rst_file| [-a]
Options:
         i  .. restart file containing result of TrackingJobMaster.getRst()
         a  .. add hosts to PVM
"""
        sys.exit(0)

## MAIN ##

use()

cmd = T.cmdDict()

T.flushPrint('Loading restart data...')
rst = T.Load( cmd['i'] )

hosts = []
if 'a' in cmd:
    hosts = [ h['host'] for h in rst['hosts'] ]

master = restart( rst, hosts=hosts )

T.flushPrint('Master initialized for restart.')

master.start()





