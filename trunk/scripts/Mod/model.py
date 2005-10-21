#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## $Date$
## $Revision$

import Biskit.Mod.modUtils as modUtils
from Biskit.Mod import *
from Biskit.Mod import Modeller as M
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

import sys, os.path

def _use( o ):
    print """
Build model using Modeller.

Syntax: model.py [ -o |outFolder| -log |logFile| -h |host_computer| ]

Options:
    -o       output folder for results      (default: .)
    -log     log file                       (default: STDOUT)
    -h       host computer for calculation  (default: local computer)
             -> must be accessible w/o password via ssh, check!
    -? or help .. this help screen

input: templates/modeller/*.pdb
       t_coffee/final.pir_aln

output: modeller/modeller.log
                /*.B9999000??   <- models

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
    sys.exit(0)


def defaultOptions():
    return {'o':'.',
            'log': None,
            'h':None
            }

### MAIN ###

options   = tools.cmdDict( defaultOptions() )
outFolder = tools.absfile( options['o'] )
host = options['h']

log = None
if options['log']:
    log = LogFile( outFolder + '/' + options['log'], 'a' ) 

if '?' in options or 'help' in options:
    _use( defaultOptions() )


###################
## Modeller
##
## Build model using Modeller.

## input: templates/modeller/*.pdb
##        t_coffee/final.pir_aln
##
## output: modeller/modeller.log
##                 /*.B9999000??   <- models
try:

    m8 = M(outFolder, log)
   
    r = m8.prepare_modeller( )

    m8.go(host)

    m8.postProcess()
    
except:
    EHandler.error( 'Error while modelling.')
