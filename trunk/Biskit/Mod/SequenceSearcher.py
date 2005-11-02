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
## last $Author$
## last $Date$
"""
Search one or more sequence databases (SwissProt, Tremble) with a sequence
using the method of choice (Blast , FastA)
"""

import Bio.SwissProt.SProt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIStandalone
from Bio import Fasta

import Biskit.tools as tools
from Biskit.Errors import *
from Biskit import StdLog, EHandler

from Biskit.Mod import settings

import re, os
import commands
import copy
from sys import *

class BlastError( BiskitError ):
    pass

class InternalError( BiskitError ):
    pass

class SequenceSearcher:
    """
    Take a sequence and return a list of nonredundant homolog
    sequences as fasta file.

    ToDo: copy blast output

    =============================== localBlast ==============================
    
    localBlast( seqFile, db, method='blastp', resultOut=None, e=0.01, **kw )
        seqFile   - str, file name with search sequence as FASTA
        db        - list of str, e.g. ['swissprot', 'pdb']
        method    - str, e.g. 'blastp', 'fasta'
        resultOut - str, save blast output to this new file
        e         - expectaion value cutoff
        
    You may pass more parameters to **kw to change the behavior of
    the search.  Otherwise, default values will be chosen by blastall:
    
    --- Scoring ---
    matrix          Matrix to use (default BLOSUM62).
    gap_open        Gap open penalty (default 0).
    gap_extend      Gap extension penalty (default 0).

    --- Algorithm ---
    gapped          Whether to do a gapped alignment. T/F (default T)
    wordsize        Word size (blastp default 11).
    keep_hits       Number of best hits from a region to keep (default off).
    xdrop           Dropoff value (bits) for gapped alignments (blastp default 25).
    hit_extend      Threshold for extending hits (blastp default 11).

    --- Processing ---
    filter          Filter query sequence? (T/F, default F)
    restrict_gi     Restrict search to these GI's.
    believe_query    Believe the query defline? (T/F, default F)
    nprocessors     Number of processors to use (default 1).
    
    --- Formatting ---
    alignments      Number of alignments. (default 250)


    ============================== localPSIBlast =============================
    
    localPSIBlast(seqFile, db, method='blastpgp', resultOut=None, e=0.01, **kw)
        seqFile   - str, file name with search sequence as FASTA
        db        - list of str, e.g. ['swissprot', 'pdb']
        method    - str, e.g. 'blastp', 'fasta'
        e         - float, expectation value cutoff
        resultOut - str, save blast output to this new file

    You may pass more parameters to **kw to change the behavior of
    the search.  Otherwise, default values will be chosen by blastpgp:

    --- Scoring --- 
    matrix           Matrix to use (default BLOSUM62).
    gap_open         Gap open penalty (default 11).
    gap_extend       Gap extension penalty (default 1).
    window_size      Multiple hits window size (default 40).
    npasses          Number of passes (default 1).
    passes           Hits/passes (Integer 0-2, default 1).

    --- Algorithm --- 
    gapped           Whether to do a gapped alignment (T/F, default T).
    wordsize         Word size (default 3).
    keep_hits        Number of beset hits from a region to keep (default 0).
    xdrop            Dropoff value (bits) for gapped alignments (default 15).
    hit_extend       Threshold for extending hits (default 11).
    nbits_gapping    Number of bits to trigger gapping (default 22).
    pseudocounts     Pseudocounts constants for multiple passes (default 9).
    xdrop_final      X dropoff for final gapped alignment (default 25).
    xdrop_extension  Dropoff for blast extensions (default 7).
    model_threshold  E-value threshold to include in multipass model (default 0.005).
    required_start   Start of required region in query (default 1).
    required_end     End of required region in query (default -1).

    --- Processing --- 
    filter           Filter query sequence with SEG? (T/F, default F)
    believe_query    Believe the query defline? (T/F, default F)
    nprocessors      Number of processors to use (default 1).

    --- Formatting --- 
    alignments         Number of alignments (default 250).


    ==================== remoteBlast, remotePSIBlast =========================

    NCBIWWW.blast(program, database, query[, query_from][, query_to]
        [, entrez_query][, filter][, expect]
        [, word_size][, other_advanced][, cdd_search]
        [, composition_based_statistics][, matrix_name][, run_psiblast]
        [, i_thresh][, genetic_code][, show_overview][, ncbi_gi]
        [, format_object][, format_type][, descriptions][, alignments]
        [, alignment_view][, auto_format][, cgi][, timeout]) -> handle

    Blast against the NCBI Blast web page.  This uses the NCBI web
    page cgi script to BLAST, and returns a handle to the
    results. See:
    
    http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.html
    
    for more descriptions about the options.

    Required Inputs:
      program  - The name of the blast program to run (ie. blastn, blastx...)
      database - The database to search against (ie. nr, dbest...)
      query    - The input type for the search, which NCBI tries to guess.
                 Ideally, this would be a sequence in FASTA format.

    General Options:
      filter, expect, word_size, other_advanced

    Formatting Options:
      show_overview, ncbi_gi, format_object, format_type, descriptions,
      alignments, alignment_view, auto_format

    Protein specific options:
      cdd_search, composition_based_statistics, matrix_name, run_psiblast,
      i_thresh

    Translated specific options:
      genetic code

    """

    F_RESULT_FOLDER = '/sequences'

    ## standard file name for unclusterd sequences
    F_FASTA_ALL = F_RESULT_FOLDER + '/all.fasta'
    
    ## standard file name for non-redundant seqs
    F_FASTA_NR =  F_RESULT_FOLDER + '/nr.fasta'

    ## clusters
    F_CLUSTER_RAW = F_RESULT_FOLDER + '/cluster_raw.out'
    F_CLUSTER_LOG = F_RESULT_FOLDER + '/cluster_result.out'

    ## pseudo blast output
    F_BLAST_OUT = F_RESULT_FOLDER + '/blast.out'

    ## default location of existing fasta file with target sequence
    F_FASTA_TARGET = '/target.fasta'

    def __init__(self, outFolder='.', clusterLimit=50, verbose=0, log=None ):
        """
        clusterLimit - int, maximal number of returned sequence clusters [50]
        verbose - 1||0, keep temporary files [0]
        log     - LogFile, log file instance, if None, STDOUT is used [None]
        """
        self.outFolder = tools.absfile( outFolder )

        ##
        ## NOTE: If you get errors of the type "Couldn't find ID in"
        ##       check these regexps (see getSequenceIDs function)!
        ##
        self.ex_gi    = re.compile( '^>ref|gb|ngb|emb|dbj|prf\|{1,2}([A-Z_0-9.]+)\|' )
        self.ex_swiss = re.compile(  '^>sp\|([A-Z0-9_]{5,7})\|' )
        self.ex_pdb   = re.compile( '^>pdb\|([A-Z0-9]{4})\|' )

        self.verbose = verbose
        self.log = log or StdLog()

        self.record_dic = None ## dict ID - Bio.Fasta.Record
        self.clusters = None ## list of lists of ids
        self.bestOfCluster = None ## list of non-redundant ids

        self.clusterLimit = 50 ## the maximal number of clusters to return
        self.clustersCurrent = None ## current number of clusters

        self.prepareFolders()


    def prepareFolders( self ):
        if not os.path.exists(self.outFolder + self.F_RESULT_FOLDER):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def getRecords( self ):
        """
        -> [ Bio.Fasta.Record ], all found protein homologues
        """
        if not self.record_dic:
            raise BlastError( "No sequences found (yet)." )
        
        return self.record_dic.values()

    def getClusteredRecords( self ):
        """
        -> [ Bio.Fasta.Record ], best record of each cluster
        !! BlastError, if called before clustering
        """
        if not self.clusters:
            raise BlastError( "Sequences are not yet clustered." )

        return [ self.record_dic[id] for id in self.bestOfCluster ]


    def __blast2dict( self, parsed_blast, db ):
        """
        blast2dic( Bio.Blast.Record.Blast ) -> { str:Bio.Fasta.Record }
        Convert parsed blast result into dictionary of FastaRecords indexed
        by sequence ID.
        """
        ids = self.getSequenceIDs( parsed_blast )
        self.record_dic = self.fastaFromIds( db, ids )

        if self.verbose:
            self.writeFasta( self.record_dic.values(),
                             self.outFolder + self.F_FASTA_ALL )

            self.__writeBlastResult( parsed_blast,
                                   self.outFolder + self.F_BLAST_OUT)


    def remoteBlast( self, seqFile, db, method, e=0.01, **kw ):
        """
        seqFile- str, file name with search sequence as FASTA
        db     - list of str, e.g. ['swissprot', 'pdb']
        method - str, e.g. 'blastp', 'fasta'
        e      - float, expectation value cutoff
        **kw   - blast options
        """
        fasta = Fasta.Iterator( open(seqFile) )
        query = fasta.next()
        
        blast_result = NCBIWWW.qblast( method, db, query, expect=e, **kw)

        p = NCBIWWW.BlastParser()
        parsed = p.parse( blast_result )

        self.__blast2dict( parsed, db )


    def remotePSIBlast( self, seqFile, db, method, e=0.01, **kw ):
        """
        seqFile- str, file name with search sequence as FASTA
        db     - list of str, e.g. ['swissprot', 'pdb']
        method - str, e.g. 'blastp', 'fasta'
        e      - float, expectation value cutoff
        **kw   - blast options
        """
        fasta = Fasta.Iterator( open(seqFile) )
        query = fasta.next()
        
        blast_result = NCBIWWW.blast( method, db, query, run_psiblast=1,
                                      expect=e, **kw)

        ## DOES THAT WORK???
        p = NCBIWWW.BlastParser()
        parsed = p.parse( blast_result )

        self.__blast2dict( parsed, db )
    

    def localBlast( self, seqFile, db, method='blastp',
                    resultOut=None, e=0.01, **kw ):
        """
        seqFile- str, file name with search sequence as FASTA
        db     - list of str, e.g. ['swissprot', 'pdb']
        method - str, e.g. 'blastp', 'fasta'
        e - float, expectation value cutoff
        resultOut - str, save blast output to this new file
        """
        results = err = p = None
        try:
            self.l = {"blast.bin":settings.blast_bin,
                      "m":method,
                      "db":db,
                      "s":seqFile,
                      "e":e,
                      "k":kw}

            results, err = NCBIStandalone.blastall( settings.blast_bin,
                                                    method, db, seqFile,
                                                    expectation=e, **kw)
                                                                
            p = NCBIStandalone.BlastParser()

            parsed = p.parse( results )
            
            self.__blast2dict( parsed, db )

        except Exception, why:
            self.error = { 'results':results,'err':err.read(), 'parser':p,
                           'blast_bin':settings.blast_bin, 'seqFile':seqFile,
                           'method':method, 'db':db,'kw':kw,
                           'exception':str(why) }
            print tools.lastErrorTrace()
            raise BlastError( str(why) + "\n" +
                              str( self.error ) ) 


    def localPSIBlast( self, seqFile, db, method='blastpgp',
                       resultOut=None, e=0.01, **kw ):
        """
        seqFile- str, file name with search sequence as FASTA
        db     - list of str, e.g. ['swissprot', 'pdb']
        method - str, e.g. 'blastp', 'fasta'
        e - float, expectation value cutoff
        resultOut - str, save blast output to this new file
        """
        try:
            results, err = NCBIStandalone.blastpgp( settings.psi_blast_bin,
                                                        method, db, seqFile,
                                                        expectation=e, **kw)

            p = NCBIStandalone.PSIBlastParser()
            parsed = p.parse( results )

            self.__blast2dict( parsed, db )
            
        except Exception, why:
           self.error = { 'results':results,'err':err, 'parser':p }
           tools.errWriteln( "BlastError " + tools.lastErrorTrace() )
           raise BlastError( why ) 


    def getSequenceIDs( self, blast_records ):
        """
        getSequenceIDs(  Bio.Blast.Record.Blast ) -> [ str ]
        extract sequence ids from BlastParser result.
        """
        result = []
        ids = []
        for a in blast_records.alignments:

            for pattern in [self.ex_gi, self.ex_swiss, self.ex_pdb]:
                ids = pattern.findall( a.title )
                if ids:
                    break

            if not ids:
                raise BlastError( "Couldn't find ID in " + a.title)
            
            result += [ ids[0] ]

        return result


    def fastaRecordFromId( self, db, id ):
        """
        fastaRecordFromId( db, id ) -> Bio.Fasta.Record
        !! BlastError
        """
        cmd = settings.fastacmd_bin + ' -d %s -s %s' % (db, id)

        err, o = commands.getstatusoutput( cmd )
        if err:
            raise BlastError( err )
        
        frecord = Fasta.Record()
        frecord.title = id
        
        try:
            for line in o.split('\n'):
                if line[0] == '>':
                    frecord.annotation = line[1:]
                else:
                    frecord.sequence += line.strip()

        except IndexError:
            raise InternalError, "Couldn't fetch fasta record %s from database %s" \
                  % (id,db)

        return frecord
        

    def fastaFromIds( self, db, id_lst, fastaOut=None ):
        """
        fastaFromIds( id_lst, fastaOut )
        -> { str: Bio.Fasta.Record }
        """
        result = {}
        for i in id_lst:
            try:
                r = self.fastaRecordFromId( db, i )
                result[i] = r
            except BlastError, why:
                tools.errWriteln("ERROR (ignored): couldn't fetch "+ str(i) )

        return result


    def copyClusterOut( self, raw=None):
        if self.verbose and raw:
            f = open( self.outFolder + self.F_CLUSTER_RAW, 'w', 1)

            f.write(raw)
            f.close()


    def reportClustering( self, raw=None ):
        """
        """
        try:
            if self.verbose:
                f = open( self.outFolder + self.F_CLUSTER_LOG, 'w', 1)

                for cluster in self.clusters:
                    f.write( "%i\t%s\n" % ( len( cluster ), str( cluster )))

                f.close()

                self.copyClusterOut( raw=raw )
                
        except IOError, why:
            EHandler.warning( "Can't write cluster report." + str(why) )


    def clusterFasta( self, fastaIn=None, simCut=1.75, lenCut=0.9, ncpu=1 ):
        """
        fastaClust( fastaIn [, simCut, lenCut, ncpu] )

        Cluster sequences. The input fasta titles must be the IDs.
        fastaIn - str, file name
        simCut  - double, similarity threshold (score < 3 or %identity) [1.75]
        lenCut  - double, length threshold [ 0.9 ]
        ncpu    - int, number of CPUs
        !! BlastError
        """
        fastaIn = fastaIn or self.outFolder + self.F_FASTA_ALL

        if tools.fileLength( fastaIn ) < 1:
            raise IOError( "File %s empty. Nothing to cluster"%fastaIn )
        
        self.log.add( "\nClustering sequences\n"+20*"-")
        self.log.add( "sequences: \n%s\n ..." % fastaIn ) 

        cmd = settings.blastclust_bin + ' -i %s -S %f -L %f -a %i' %\
              (fastaIn, simCut, lenCut, ncpu)

        if self.verbose:
            self.log.add("\nclustering command:")
            self.log.add(cmd)

        ## bugfix: at all cost prevent blastclust from using shared temp folder
        tmp = os.environ.get( 'TMPDIR', None )
        if tmp:
            del os.environ['TMPDIR']

        err, o = commands.getstatusoutput( cmd )
        if err:
            raise BlastError( err )

        if tmp:
            os.environ['TMPDIR'] = tmp

        lines = o.split('\n')[2:]
        self.clusters = [ l.split() for l in lines ]

        self.reportClustering( raw=o )

        self.log.add( "returning %i clusters.\n" % ( len( self.clusters)))

        self.bestOfCluster = [ self.selectFasta( ids )
                                for ids in self.clusters ]


    def clusterFastaIterative(self, fastaIn=None, simCut=1.75, lenCut=0.9,
                              ncpu=1 ):
        """
        Run cluterFasta iteratively, with tighter clustering settings, until
        the number of clusters are less than self.clusterLimit.
        """
        iter = 1
        while self.clustersCurrent > self.clusterLimit \
                  or self.clustersCurrent == None:

            self.clusterFasta( fastaIn, simCut, lenCut, ncpu )
            self.clustersCurrent = len( self.clusters )
        
            self.log.add( "\nClustering iteration %i produced %i clusters." \
                          %( iter, len( self.clusters) ))
            self.log.add( "\tusing a simCut of %.2f and a lenCut of %.3f.\n" \
                          %( simCut, lenCut ) )
            
            simCut -= 0.15
            lenCut -= 0.015
            iter += 1

            
    def selectFasta( self, ids_from_cluster ):
        """
        select one member of cluster of sequences.
        ids - [ str ], list of sequence ids defining the cluster
        -> str, id of best in cluster
        """
        return ids_from_cluster[0]


    def writeFasta( self, frecords, fastaOut ):
        """
        Create fasta file for given set of records.
        frecords - [ Bio.Blast.Record ]
        fastaOut - str, file name
        """
        f = open( tools.absfile(fastaOut), 'w' )
        for r in frecords:
            f.write( str( r ) + '\n' )
        f.close()

    def writeFastaAll( self, fastaOut=None ):
        """
        Write all found template sequences to fasta file.
        fastaOut - str, [ 'sequences/all.fasta' ]
        """
        fastaOut = fastaOut or self.outFolder + self.F_FASTA_ALL
        self.writeFasta( self.frecords, fastaOut )


    def writeFastaClustered( self, fastaOut=None ):
        """
        Write non-redundant set of template sequences to fasta file.
        fastaOut - str, [ 'sequences/nr.fasta' ]
        """
        fastaOut = fastaOut or self.outFolder + self.F_FASTA_NR
        self.writeFasta( self.getClusteredRecords(), fastaOut )


    def __writeBlastResult( self, parsed_blast, outFile ):
        """
        writeBlastResult( parsed_blast, outFile )
        Simulate blast output from parsed Blast result.
        parsed_blast -  Bio.Blast.Record.Blast
        outFile  - str
        """
        f = open( tools.absfile( outFile ), 'w' )

        for (description, alignment) in \
                zip(parsed_blast.descriptions,parsed_blast.alignments):

            for hsp in alignment.hsps:
                f.write('\n****Alignment****\n')
                f.write('sequence:' + alignment.title + '\n' )
                f.write('length:' + str(alignment.length) )
                f.write('\te value: '+str(hsp.expect) + '\n')
                f.write( hsp.query[0:75] + '...\n' )
                f.write( hsp.match[0:75] + '...\n' )
                f.write( hsp.sbjct[0:75] + '...\n' )

        f.close()

      
def test():
    options = {}
    options['q'] = '/home/Bis/raik/homopipe/test/nmr_hit.fasta'
    options['o'] = '/home/Bis/johan/blats.test'
    return options
  

##########
## TEST ##
##########

if __name__ == '__main__':
    
    options = test()
    #    options = cmdDict( _defOptions() )

    searcher = SequenceSearcher(
        outFolder=tools.projectRoot()+'/test/Mod/project',
        verbose=1 )

    ## remote Blast
#    r = searcher.remoteBlast( options['q'], 'swissprot',
#                              'blastp', alignments=50 )

    ## remote PSIBlast
#    r = searcher.remotePSIBlast( options['q'], 'swissprot',
#                              'blastp', alignments=50 )

##    rl, el = searcher.localBlast( options['q'], 'pdbaa', 'blastp')

    ## local PSIBlast
##    rl, el = searcher.localPSIBlast( options['q'], 'pdbaa',
##                                     'blastpgp', npasses=2)

    db = 'swissprot'
    db = 'nr'
    f_target = searcher.outFolder + searcher.F_FASTA_TARGET
    searcher.localBlast( f_target, db, 'blastp', alignments=500, e=0.01 )

    searcher.clusterFasta()  ## expects all.fasta

    searcher.writeFastaClustered()


# fastacmd -d pdbaa -s 2098544
# blast db has to be built with -o potion
#    cat db1.dat db2.dat | formatdb -i stdin -o T -n indexed_db



