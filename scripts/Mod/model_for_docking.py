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

from Biskit.Mod import SequenceSearcher
from Biskit.Mod import TemplateSearcher
from Biskit.Mod import TemplateCleaner
from Biskit.Mod import Aligner
from Biskit.Mod import Modeller

import Biskit.Mod.modUtils as modUtils
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile


## def test():
##     options = {}
##     dir = '/home/Bis/raik/data/tb/interfaces/mod_c01/rec_mod/'
##     options['q'] = dir + 'target.fasta'
##     options['o'] = dir
##     return options
  
def testOptions():
    options = {}
    options['q'] = None
    options['o'] = '/home/Bis/raik/data/mod_project_1/'
    options['log'] = 'modelling.log'
    return options

def defaultOptions():
    return {'q':None,
            'o':'.' }

options = tools.cmdDict( defaultOptions() )
options = testOptions()

outFolder = tools.absfile( options['o'] )
f_target = tools.absfile( options['q'] )
f_target = f_target or outFolder + SequenceSearcher.F_FASTA_TARGET

log = None
if 'log' in options:
    log = LogFile( outFolder + '/' + options['log'] ) 

seq_db = 'sprot.dat'
tmp_db = 'pdbaa'


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

try:
    searcher = SequenceSearcher( outFolder=outFolder, verbose=1 )

    ## local PSIBlast - not fully implemented!!
    #searcher.localPSIBlast( target, seq_db, e=0.1, alignments=1000)

    ## local Blast
    searcher.localBlast( f_target, seq_db, 'blastp', alignments=500, e=0.0001 )

    ## cluster blast results. Defaults: simCut=1.75, lenCut=0.9, ncpu=1
    ## expects all.fasta
    searcher.clusterFasta() # simCut=0.45, lenCut=0.9, ncpu=1 )  

    searcher.writeFastaClustered()

except:
    EHandler.error( 'Error while searching for homologues.')

###################
## TemplateSearcher
##
## Find modelling templates, blasting the target sequence against "tmp_db"
## Cluster the sequences and download the pdbs to templates/all

## input: target.fasta
##
## output: templates/blast.out
##         templates/all.fasta
##         templates/cluster_result.out
##         templates/nr.fasta              (input for Aligner)
##         templates/all/*.pdb
##         templates/nr/chain_index.txt    (input for TemplateCleaner)
##                     /*.pdb              (  "    "         "       )

try:
    searcher = TemplateSearcher( outFolder, verbose=1 )

    searcher.localBlast(f_target, tmp_db, 'blastp', alignments=200, e=0.001)

    searcher.retrievePDBs()

    ## expects all.fasta
    searcher.clusterFasta( simCut=1.75, lenCut=0.9, ncpu=1 )  

    searcher.writeFastaClustered()

    fn = searcher.saveClustered()

except:
    EHandler.error( 'Error while searching for templates')

###################
## TemplateCleaner
##
## Prepare pdb files in templates/nr for T-coffee and modeller
## (replace nonstandard residues, remove hydrogens,
#    remove atoms with nultiple configurations, etc.)

## input: templates/nr/*.pdb
##        templates/nr/chain_index.txt
##
## output: templates/t_coffee/*.alpha    (input for Alignar)
##         templates/modeller/*.pdb      (input for Modeller)

try:
    cleaner = TemplateCleaner( outFolder, log )

    inp_dic = modUtils.parse_tabbed_file( outFolder +
                                          TemplateSearcher.F_NR +
                                          TemplateSearcher.F_CHAIN_INDEX )

    cleaner.process_all( inp_dic )

except:
    EHandler.error( 'Error while cleaning templates')


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

## note 2: If there is only one template structure step 2 of T-coffee
##         will not work. Solution, skipp the structural alignment if
##         only one teemplate structure is provided.

## note 3: In quite som cases the sequence retrieved from the nrpdb
##         sequence database is different from the sequence extracted
##         from the coordinates in the pdb-file. This will sometimes cause
##         t-coffee to terminate with an error (two sequences with the same name
##         but with different sequences). Temporary solution: Choose another
##         structure from the same cluster as the troublemaker.

try:
    a = Aligner( outFolder, log )

    a.align_for_modeller_inp()

    a.go('magneto')
except:
    EHandler.error( 'Error while building alingnments.')


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
    m = Modeller( outFolder, log )

    r = m.prepare_modeller( )

    m.go('magneto')
except:
    EHandler.error( 'Error while modelling.')
