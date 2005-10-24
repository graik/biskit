#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
Syntax: search_sequences.py [-q |target.fasta| -o |outFolder| -log |logFile|
               -db |database| -e |e-value-cutoff| -aln |n_alignments|
               -psi
               -... additional options for blastall (see SequenceSearcher.py) ]

Result: folder 'sequences' with blast result, all found sequences in all.fasta,
        clustered sequences in nr.fasta
        
Options:
    -q       fasta file with query sequence (default: ./target.fasta)
    -o       output folder for results      (default: .)
    -log     log file                       (default: STDOUT)
    -db      sequence data base
    -e       E-value cutoff for sequence search
    -aln     number of alignments to be returned
    -simcut  similarity threshold for blastclust (score < 3 or % identity)
    -simlen  length threshold for clustering
    -ncpu    number of CPUs for clustering
    -psi     use PSI Blast instead, experimental!!

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value
    sys.exit(0)


def defaultOptions():
    return {'q':None,
            'o':'.',
            'db' : 'sprot.dat',
            'log': None,
            'e':0.01,
            'aln':500,
            'simcut':1.75,
            'simlen':0.9,
            'ncpu':1
            }

def blastOptions( options ):
    """
    identify options that have to be passed on to blastall
    """
    result = {}
    def_keys = defaultOptions().keys() + ['psi']

    for k, v in options.items():
        if not k in def_keys:
            result[ k ] = v

    return result

### MAIN ###

options   = tools.cmdDict( defaultOptions() )
outFolder = tools.absfile( options['o'] )
f_target  = tools.absfile( options['q'] )
f_target = f_target or outFolder + SequenceSearcher.F_FASTA_TARGET

if not (f_target and os.path.exists( f_target ) ):
    _use( defaultOptions() )

if '?' in options or 'help' in options:
    _use( defaultOptions() )

seq_db = options['db']
e = float( options['e'] )
aln = int( options['aln'])
simCut = float( options['simcut'] )
simLen = float( options['simlen'] )
nCpu = int( options['ncpu'] )

log = None
if options['log']:
    log = LogFile( outFolder + '/' + options['log'], 'a' ) 

ext_options = blastOptions( options )

###################
## SequenceSearcher
##
## Find homologues to the target sequence using blast against "seq_db"
## Cluster the seuences and write the result to nr.fasta

## input: target.fasta
## 
## output: sequences/all.fasta
##                  /blast.out
##                  /cluster_result.out
##                  /nr.fasta                 (input for Aligner)

searcher = SequenceSearcher( outFolder=outFolder, verbose=1, log=log )

if 'psi' in options:
    ## local PSIBlast - not fully implemented!!
    searcher.localPSIBlast( f_target, seq_db, e=e, alignments=aln,
                            **ext_options )
else:
    ## local Blast
    searcher.localBlast( f_target, seq_db, 'blastp', alignments=aln, e=e,
                         **ext_options )

## cluster blast results. Defaults: simCut=1.75, lenCut=0.9, ncpu=1
## expects all.fasta
#searcher.clusterFasta( simCut=simCut, lenCut=simLen, ncpu=nCpu )
searcher.clusterFastaIterative( simCut=simCut, lenCut=simLen, ncpu=nCpu )
searcher.writeFastaClustered()
