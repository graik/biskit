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


from Biskit.Mod import SequenceSearcher
from Biskit.Mod import TemplateSearcher
from Biskit.Mod import TemplateCleaner
from Biskit.Mod import Aligner
from Biskit.Mod import Modeller

import Biskit.Mod.modUtils as modUtils
import Biskit.tools as tools
from Biskit import EHandler
from Biskit import LogFile

import Biskit.Trajectory as Trajectory
import Biskit.Pymoler as Pymoler
import glob
import sys


options = {'q':None,
           'o':'.'}

def _use():
    print """
Biskit.Mod example script that models a structure from a fasta
formated sequence file in 4 steps:

1) Searches for homologe sequences and clusters the found
   sequences to a representative set using NCBI-Tools.
2) Searches for temptale structures for the homology modeling.
   Similar structures are removed by clustering.
3) Build a combined sequence/structure alignment using T-Coffee.
4) Build models using Modeller.

Syntax: modelling_example.py -q |query file| -o |outputFolder|
                            [-h |host| -log  -view ]

Options:
   -q     file; fasta formated sequence file to model
   -o     folder; directory in which all project files will be
            written
   -h     host name; the quite cpu consuming stasks of aligning
            and modeling can be sent to a remote host that also
            has access to the output directory
   -log   write stdOut messages to log file (~project/modelling.log)
   -view  show the superimposed models in PyMol


HINT: If you want to inspect the alignment used for modeling: 
      ~project/t_coffee/final.score_html 


Default options:"""
    for key in options.keys():
        print "\t-",key, "\t",options[key]

    sys.exit(0)


def testOptions():
    options = {}
    options['q'] = tools.testRoot()+ '/Mod/project/target.fasta'
    options['o'] = tools.testRoot()+ '/Mod/project'
    options['log'] = '1'
    options['view'] = '1'
    return options


###########################
# MAIN
###########################

if len(sys.argv) < 3:
    _use()
    
options = tools.cmdDict( options )
#options = testOptions()

outFolder = tools.absfile( options['o'] )
f_target = tools.absfile( options['q'] )
f_target = f_target or outFolder + SequenceSearcher.F_FASTA_TARGET

log = None
if 'log' in options:
    log = LogFile( outFolder + '/modelling.log'  )


## databases used
seq_db = 'swissprot'
tmp_db = 'pdbaa'


###############
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

tools.flushPrint('Searching for homologues ...')

try:
    # initiate 
    searcher = SequenceSearcher( outFolder=outFolder, verbose=1, log=log )

#    ## local PSIBlast - not fully implemented!!
#    searcher.localPSIBlast( target, seq_db, e=0.1, alignments=1000)

    ## local Blast
    searcher.localBlast( f_target, seq_db, 'blastp', alignments=500, e=0.0001 )

    ## cluster blast results. Defaults: simCut=1.75, lenCut=0.9, ncpu=1
    ## expects all.fasta
#    searcher.clusterFastaIterative( )
    searcher.clusterFasta() 

    searcher.writeFastaClustered()
    
    tools.flushPrint('Done.\n')
    
except:
    EHandler.error( 'Error while searching for homologues.')
    

###############
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

tools.flushPrint('Searching for templates ...')

try:
    searcher = TemplateSearcher( outFolder, verbose=1, log=log )

    searcher.localBlast(f_target, tmp_db, 'blastp', alignments=200, e=0.001)

    searcher.retrievePDBs()

    ## expects all.fasta
    searcher.clusterFasta( simCut=1.75, lenCut=0.9, ncpu=1 )  

    searcher.writeFastaClustered()

    fn = searcher.saveClustered()

    tools.flushPrint('Done.\n')
    
except:
    EHandler.error( 'Error while searching for templates')

    
###################
## TemplateCleaner
##
## Prepare pdb files in templates/nr for T-coffee and modeller
## (replace nonstandard residues, remove hydrogens,
##   remove atoms with multiple configurations, etc.)

## input: templates/nr/*.pdb
##        templates/nr/chain_index.txt
##
## output: templates/t_coffee/*.alpha    (input for Alignar)
##         templates/modeller/*.pdb      (input for Modeller)

tools.flushPrint('Cleaning template structures...')

try:
    cleaner = TemplateCleaner( outFolder, log=log )

    inp_dic = modUtils.parse_tabbed_file( outFolder +
                                          TemplateSearcher.F_NR +
                                          TemplateSearcher.F_CHAIN_INDEX )

    cleaner.process_all( inp_dic )

    tools.flushPrint('Done.\n')
                     
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
##         Remedy: Use clusterFastaIterative() when searching for sequences
    
## note 2: If there is only one template structure step 2 of T-coffee
##         will not work. Solution, skipp the structural alignment if
##         only one template structure is provided.

## note 3: In quite som cases the sequence retrieved from the nrpdb
##         sequence database is different from the sequence extracted
##         from the coordinates in the pdb-file. This will sometimes
##         cause t-coffee to terminate with an error (two sequences
##         with the same name but with different sequences). Temporary
##         solution: Choose another  structure from the same cluster
##         as the troublemaker.

tools.flushPrint('Creating sequence/structure alignment...')

try:
    a = Aligner( outFolder, log=log )

    a.align_for_modeller_inp()

    if options.has_key('host'):
        a.go(options['host'])
    else:        
        a.go( )

    tools.flushPrint('Done.\n')
    
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

tools.flushPrint('Building model...')

try:
    m = Modeller( outFolder, log=log )

    r = m.prepare_modeller( )

    if options.has_key('host'):
        m.go(options['host']) 
    else:
        m.go()

    m.postProcess()
    
    tools.flushPrint('Done.\n')
    
except:
    EHandler.error( 'Error while modelling.')



###################
##
## Superimpose the models and look at the result in Pymol
##

def __printMatrix( matrix ):
    """
    Print the part right of the diagonal in a matrix
    """
    nr = len( matrix )
    for i in range(nr): print '%5i'%(i+1),
    for i in range(nr):
        print '\n%2i'%(i+1),
        for k in range(i):
            print ' '*5,
        for j in range(i, nr):
            print '%5.2f'%matrix[i,j],
            
    
## get filenames of all models 
models = glob.glob( '%s/modeller/%s*.pdb'%(outFolder,
                                           tools.stripFilename(f_target)) )

## create a Trajectory object with the models
traj = Trajectory( pdbs=models )

## fit the models against the average structure iteratively
traj.blockFit2ref()

## calculate and print rmsd matrix
rmsHeavy = traj.pairwiseRmsd()
print '\nHEAVY ATOM RMSD BETWEEN MODELS::'
__printMatrix( rmsHeavy )

## same thing for backbone atoms
BBMask = traj[0].maskBB()
traj.blockFit2ref( mask = BBMask )
rmsBB = traj.pairwiseRmsd( aMask = BBMask  )
print '\nBACKBONE RMSD BETWEEN MODELS:'
__printMatrix( rmsBB )

if options.has_key('view'):
    
    ## show backbone superimposed structures in pymol
    pm = Pymoler( )
    for t in traj:
        pm.addPdb( t )

    pm.show()
