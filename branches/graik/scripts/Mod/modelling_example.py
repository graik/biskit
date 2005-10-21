from Biskit.Mod import SequenceSearcher, TemplateSearcher, Aligner, Modeller
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
import Biskit.tools as tools

## def test():
##     options = {}
##     options['q'] = '/home/Bis/johan/mod_test/target.fasta'
##     options['o'] = '/home/Bis/johan/mod_test'
##     return options


def test():
    options = {}
    options['q'] = '/home/Bis/oreste/fichier5'
    options['o'] = '/home/Bis/oreste/fichier5/target.fasta'
    return options


  
options = test()

outFolder = outFolder=tools.absfile( options['o'] )
query_seq = options['q']

seq_db = 'sprot.dat'
tmp_db = 'pdbaa'

###############
# SequenceSearcher

searcher = SequenceSearcher( outFolder=outFolder, verbose=1 )

## ## remote Blast
## r = searcher.remoteBlast( options['q'], 'swissprot',
##                           'blastp', alignments=50 )

## ## remote PSIBlast
## r = searcher.remotePSIBlast( options['q'], 'swissprot',
##                              'blastp', alignments=50 )

## ## local PSIBlast
## rl, el = searcher.localPSIBlast( options['q'], 'pdbaa',
##                                  'blastpgp', npasses=2)


## file name of target sequence
ftarget = searcher.outFolder + searcher.F_FASTA_TARGET

## local Blast
searcher.localBlast( ftarget, seq_db, 'blastp', alignments=500, e=0.01 )

searcher.clusterFasta()  ## expects all.fasta

searcher.writeFastaClustered()


###############
# TemplateSearcher

searcher = TemplateSearcher( outFolder, verbose=1 )

searcher.localBlast( query_seq, tmp_db, 'blastp', alignments=200, e=0.0001)

searcher.retrievePDBs()

searcher.clusterFasta()  ## expects all.fasta

searcher.writeFastaClustered()

fn = searcher.saveClustered()


###############
# Aligner

a = Aligner( outFolder )

a.align_for_modeller_inp()

a.go('scarescrow')


###############
# Aligner

m = Modeller( outFolder )

r = m.prepare_modeller( )
