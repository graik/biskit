#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## $Date$
## $Revision$

import Biskit.Mod.modUtils as modUtils
from Biskit.Mod import *
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

import sys, os.path

def _use( o ):
    print """
Syntax: align.py [ -o |outFolder| -log |logFile| -h |host_computer| ]

Options:
    -o       output folder for results      (default: .)
    -log     log file                       (default: STDOUT)
    -h       host computer for calculation  (default: local computer)
             -> must be accessible w/o password via ssh, check!
    -? or help .. this help screen

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
## Aligner
##
## Create a sequence-structure alignment using T-coffee.
## Convert the alignment into Modeller compatible format 

## input: sequences/nr.fasta
##        templates/templates.fasta
##        templates/t_cofee/*.alpha
##!
## output: t_coffee/fast_pair.lib
##                 /final.score_html
##                 /struct.aln
##                 /t_coffee.log_*
##                 /final.aln
##                 /lalign_id_pair.lib
##                 /struct.aln_original
##                 /final.phylip
##                 /sap_pair.lib
##                 /t_coffee.inp
##                 /final.pir_aln             (input for Modeller)
##                 /sap_pair.lib_original

## note 1: If there are more than approximately 50 sequences overall
##         t_coffe will eat all the memory and the job will not finish

## (note 2:) Now this is taken care of by Aligner.py
##         If there is only one template structure step 2 of T-coffee
##         will not work. Solution, skipp the structural alignment if
##         only one teemplate structure is provided.

## note 3: In quite som cases the sequence retrieved from the nrpdb
##         sequence database is different from the sequence extracted
##         from the coordinates in the pdb-file. This will sometimes cause
##         t-coffee to terminate with an error (2 sequences with the same name
##         but with different sequences). Temporary solution: Choose another
##         structure from the same cluster as the troublemaker.
try:
    a = Aligner( outFolder, log )

    a.align_for_modeller_inp()

    a.go(host)

except:
    EHandler.error( 'Error while building alingnments.')
    print "\nalign.py -? or align.py -help for help screen"
