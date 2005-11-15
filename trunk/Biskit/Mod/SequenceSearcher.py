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
## $Revision$
"""
Search one or more sequence databases (SwissProt, Tremble) with a sequence
using the method of choice (Blast , FastA)
"""

import Bio.SwissProt.SProt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
from Bio import Fasta

import Biskit.tools as tools
from Biskit.Errors import *
from Biskit import StdLog, EHandler

from Biskit.Mod import settings

import re, os
import string, copy
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
    
    localPSIBlast(seqFile, db, resultOut=None, e=0.01, **kw)
        seqFile   - str, file name with search sequence as FASTA
        db        - list of str, e.g. ['swissprot', 'pdb']
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
    NCBIWWW.qblast

    Do a BLAST search using the QBLAST server at NCBI.
      program        BLASTP, BLASTN, BLASTX, TBLASTN, or TBLASTX.
      database       Which database to search against.
      sequence       The sequence to search.
      ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.  Def FALSE.
      descriptions   Number of descriptions to show.  Def 500.
      alignments     Number of alignments to show.  Def 500.
      expect         An expect value cutoff.  Def 10.0.
      matrix         Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
      filter         'none' turns off filtering.  Default uses 'seg' or 'dust'.
      format_type    'HTML', 'Text', 'ASN.1', or 'XML'.  Def. 'HTML'

    """

    F_RESULT_FOLDER = '/sequences'

    ## standard file name for unclusterd sequences
    F_FASTA_ALL = F_RESULT_FOLDER + '/all.fasta'
    
    ## standard file name for non-redundant seqs
    F_FASTA_NR =  F_RESULT_FOLDER + '/nr.fasta'

    ## clusters
    F_CLUSTER_RAW = F_RESULT_FOLDER + '/cluster_raw.out'
    F_CLUSTER_LOG = F_RESULT_FOLDER + '/cluster_result.out'

    ## blast alignments
    F_BLAST_OUT = F_RESULT_FOLDER + '/blast.out'
    F_CLUSTER_BLAST_OUT = F_RESULT_FOLDER + '/cluster_blast.out'
    
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
        self.ex_gi    = re.compile( '^>ref|gb|ngb|emb|dbj|prf\|{1,2}([A-Z_0-9.]{4,9})\|' )
        self.ex_swiss = re.compile( '^>sp\|([A-Z0-9_]{5,7})\|' )
        self.ex_pdb   = re.compile( '^>pdb\|([A-Z0-9]{4})\|' )

        ##
        ##      RegExp for the remoteBlast searches. Somehow the
        ##      NCBI GI identifier is written before the identifier
        ##      we want (even though this is explicitly turned off in
        ##      SequenceSearcher.remoteBlast).
        ##
        gi = '^gi\|[0-9]+\|'
        self.ex_Rgi = re.compile( gi+'ref|gb|ngb|emb|dbj|prf\|([A-Z_0-9.]+)\|')
        self.ex_Rswiss = re.compile( gi+'sp\|([A-Z0-9_]{5,7})\|' )
        self.ex_Rpdb   = re.compile( gi+'pdb\|([A-Z0-9]{4})\|' )

        self.verbose = verbose
        self.log = log or StdLog()

        self.record_dic = None ## dict ID - Bio.Fasta.Record
        self.clusters = None ## list of lists of ids
        self.bestOfCluster = None ## list of non-redundant ids

        ## the maximal number of clusters to return
        self.clusterLimit =  clusterLimit
        
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
        ncbi_gi- TRUE/FALSE whether to give 'gi' identifier.  Def FALSE.
        **kw   - blast options

        NOTE: Using the remoteBlast is asking for trouble, as
              every change in the output file might kill the
              parser. If you still want to use remoteBlast we
              strongly recomend that you install BioPython
              from CVS. Information on how to do this you
              will find on the BioPython homepage.
        """
        fasta = Fasta.Iterator( open(seqFile) )
        query = fasta.next()
        
        blast_result = NCBIWWW.qblast( program=method, database=db,
                                       sequence=query, expect=e,
                                       ncbi_gi='FALSE', **kw)

        p = NCBIXML.BlastParser()

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


    def localPSIBlast( self, seqFile, db, rounds=3, resultOut=None,
                       e=0.01, **kw ):
        """
        seqFile- str, file name with search sequence as FASTA
        db     - list of str, e.g. ['swissprot', 'pdb']
        e - float, expectation value cutoff, def 0.001
        resultOut - str, save blast output to this new file

        e_init - float, threshold for the initial BLAST search, def. 10
        """
        results = err = p = None
        try:
            results, err = NCBIStandalone.blastpgp( settings.psi_blast_bin,
                                                    db, seqFile,
                                                    npasses=rounds,
                                                    expectation=e, **kw)
            
            p = NCBIStandalone.PSIBlastParser()

            parsed = p.parse( results )

            self.__blast2dict( parsed.rounds[-1], db )

        except Exception, why:
            self.error = { 'results':results,'err':err.read(), 'parser':p,
                           'blast_bin':settings.blast_bin, 'seqFile':seqFile,
                           'db':db,'kw':kw,'exception':str(why) }
            tools.errWriteln( "BlastError " + tools.lastErrorTrace() )
            print tools.lastErrorTrace()
            raise BlastError( str(why) + "\n" +
                              str( self.error ) ) 



    def getSequenceIDs( self, blast_records ):
        """
        getSequenceIDs(  Bio.Blast.Record.Blast ) -> [ str ]
        extract sequence ids from BlastParser result.
        """
        result = []
        ids = []
        for a in blast_records.alignments:
            for pattern in [self.ex_gi,  self.ex_swiss,  self.ex_pdb,
                            self.ex_Rgi, self.ex_Rswiss, self.ex_Rpdb]:
                ids = pattern.findall( a.title )
                if ids and ids[0]!='':
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
                tools.errWriteln("ERROR (ignored): couldn't fetch %s"%str(i) )

        return result


    def copyClusterOut( self, raw=None):
        if self.verbose and raw:
            f = open( self.outFolder + self.F_CLUSTER_RAW, 'w', 1)

            f.write(raw)
            f.close()


    def reportClustering( self, raw=None,  ):
        """
        """
        try:
            if self.verbose:
                f = open( self.outFolder + self.F_CLUSTER_LOG, 'w', 1)

                for cluster in self.clusters:

                    f.write( "%i\t%s\n" % ( len( cluster ), str( cluster )))

                f.close()

                ## write blast records of centers to disc
                centers = [ c[0] for c in self.clusters ]
                
                self.writeClusteredBlastResult( \
                    self.outFolder + self.F_BLAST_OUT,
                    self.outFolder + self.F_CLUSTER_BLAST_OUT, centers )

        
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
        
        self.log.add( "\nClustering sequences:\n%s"%('-'*20) )

        cmd = settings.blastclust_bin + ' -i %s -S %f -L %f -a %i' %\
              (fastaIn, simCut, lenCut, ncpu)

        if self.verbose:
            self.log.add("- Command: %s"%cmd)

        ## bugfix: at all cost prevent blastclust from using shared temp folder
        tmp = os.environ.get( 'TMPDIR', None )
        if tmp:
            del os.environ['TMPDIR']

        err, o = commands.getstatusoutput( cmd )
        if err:
            raise BlastError( err )

        if tmp:
            os.environ['TMPDIR'] = tmp

        ## blastclust might write errors to file, if so the errors
        ## occur before the dateline
        lines = [ l.split() for l in o.split('\n') ]
        dateline = [ l[-1] for l in lines ].index('queries')
        self.clusters = lines[dateline+1:]
        
        self.reportClustering( raw=o )

        self.bestOfCluster = [ self.selectFasta( ids )
                                for ids in self.clusters ]


    def clusterFastaIterative(self, fastaIn=None, simCut=1.75, lenCut=0.9,
                              ncpu=1 ):
        """
        Run cluterFasta iteratively, with tighter clustering settings, until
        the number of clusters are less than self.clusterLimit.

        Older versions of T-Coffee can't handle more that approx.
        50 sequences (out of memory). This problem is taken care of
        in versions > 3.2 of T-Coffee.
        """
        iter = 1
        while self.clustersCurrent > self.clusterLimit \
                  or self.clustersCurrent == None:

            self.clusterFasta( fastaIn, simCut, lenCut, ncpu )
            self.clustersCurrent = len( self.clusters )
        
            self.log.add( "- Clustering iteration %i produced %i clusters." \
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


    def __writeBlastResult( self, parsed_blast, outFile):
        """
        writeBlastResult( parsed_blast, outFile )
        parsed_blast -  Bio.Blast.Record.Blast
        outFile  - str
        """
        f = open( tools.absfile( outFile ), 'w' )

        i=1
        for alignment in parsed_blast.alignments:
            for hsp in alignment.hsps:
                s = string.replace(alignment.title,'\n',' ')
                s = string.replace(s, 'pdb|',  '\npdb|')
                f.write('Sequence %i: %s\n'%(i,s))                
                f.write('Length: %i \tScore: %3.1f \tE-value: %2.1e\n'\
                        %(hsp.identities[1], hsp.score, hsp.expect))
                f.write( 'Identities: %i \tPositives: %i \tGaps: %i\n'\
                         %(hsp.identities[0], hsp.positives[0],
                           hsp.gaps[0] or 0 ))

                f.write( '%s\n'%hsp.query  )
                f.write( '%s\n'%hsp.match )
                f.write( '%s\n\n'%hsp.sbjct )
                i += 1
        f.close()


    def writeClusteredBlastResult( self, allFile, clustFile, selection ):
        """
        Reads the blast.out file and keeps only centers.
        allFile   - file, all blast results
        clustFile - str, output file name
        selection - list, write only sequences in list
        """
        f = open( clustFile, 'w' )

        b  = open( allFile, 'r' )
        line = b.readlines()

        j=0
        sel = copy.copy(selection)
        for l in line:
            for s in sel:
                ## search for center ID, remove from search list if found
                if re.search( s, l[:30]):
                    sel.remove(s)
                    f.write( '\nCluster center %i:\n'%(selection.index(s)+1) )
                    f.writelines( l )

                    ## add following lines containing alignment info
                    while j < len(line)-1:
                        j+=1
                        ## break if empty line
                        if len( line[j] )<4:
                            break
                        ## skipp if additional pdb identifier
                        if line[j][:4] == 'pdb|':
                            pass
                        else:
                            f.writelines( line[j] )  
                
        f.close()        
        b.close()

            
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
    f_target = searcher.outFolder + searcher.F_FASTA_TARGET
    searcher.localBlast( f_target, db, 'blastp', alignments=500, e=0.01 )

    searcher.clusterFasta()  ## expects all.fasta

    searcher.writeFastaClustered()


# fastacmd -d pdbaa -s 2098544
# blast db has to be built with -o potion
#    cat db1.dat db2.dat | formatdb -i stdin -o T -n indexed_db



