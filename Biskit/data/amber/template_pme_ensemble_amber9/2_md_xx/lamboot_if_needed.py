#!/usr/bin/env python

import commands
import os

lamReturn = commands.getoutput('lamnodes')

if lamReturn[:10] == 10*'-':
    ## LAM not running
    print "no LAM nodes found -> lamboot"
    os.system( 'lamboot' )
else:
    print "found LAM nodes: ", lamReturn
