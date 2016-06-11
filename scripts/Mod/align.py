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
##

import Biskit.Mod.modUtils as modUtils
from Biskit.Mod import *
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

import sys, os.path

def _use( o ):
    print """
Syntax: align.py [ -o |outFolder| -log |logFile| -h |host_computer| -nosap ]

Options:
    -o       output folder for results      (default: .)
    -log     log file                       (default: STDOUT)
    -nosap   skip structural alignment      (default: don't skip)  
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
sap  = not 'nosap' in options

log = None
if options['log']:
    log = LogFile( outFolder + '/' + options['log'], 'a' )
    
if not (os.path.exists( outFolder +'/templates' ) ):
    print 'Current directory is not a valid modeling folder (missing /templates).' 
    _use( defaultOptions() )
          
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
##
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
##         This should be fixed in more recent versions of T-Coffee
##         (v > 3.2) where T-Coffee, according to the manual "switches
##         to a heuristic mode, named DPA, where DPA stands for Double
##         Progressive Alignment."
    
## note 2: If there is only one template structure step 2 of T-coffee
##         will not work. Solution, skip the structural alignment if
##         only one template structure is provided.

## note 3: In quite som cases the sequence retrieved from the nrpdb
##         sequence database is different from the sequence extracted
##         from the coordinates in the pdb-file. This will sometimes
##         cause t-coffee to terminate with an error (two sequences
##         with the same name but with different sequences). Temporary
##         solution: Choose another  structure from the same cluster
##         as the troublemaker.


try:
    a = Aligner( outFolder, log, verbose=1, sap=sap )

    a.align_for_modeller_inp()

    a.go(host)

except:
    EHandler.error( 'Error while building alingnments.')
    print "\nalign.py -? or align.py -help for help screen"

