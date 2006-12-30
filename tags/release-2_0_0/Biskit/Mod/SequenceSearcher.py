##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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

import Biskit.tools as T
from Biskit.Errors import *
from Biskit import StdLog, EHandler

from Biskit.Mod import settings

import re, os
import string, copy
import commands
import copy
import sys

class BlastError( BiskitError ):
    pass

class InternalError( BiskitError ):
    pass

class SequenceSearcher:
    """
    Take a sequence and return a list of nonredundant homolog
    sequences as fasta file. The search can be performed using three
    different methods:

     1. localBlast
        Uses Bio.Blast.NCBIStandalone.blastall (Biopython) to perform
        the search.

     2. localPSIBlast
        Uses Bio.Blast.NCBIStandalone.blastpgp (Biopython)
        
     3. remoteBlast
        Uses Bio.Blast.NCBIWWW.qblast (Biopython) which performs
        a BLAST search using the QBLAST server at NCBI.

    @note: the blast sequence database has to be built with -o potion
           e.g. cat db1.dat db2.dat | formatdb -i stdin -o T -n indexed_db
    
    @todo: copy blast output
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
        @param outFolder: project folder (results are put into subfolder) ['.']
        @type  outFolder: str
        @param clusterLimit: maximal number of returned sequence clusters
                             (default: 50)
        @type  clusterLimit: int
        @param verbose: keep temporary files (default: 0)
        @type  verbose: 1|0
        @param log: log file instance, if None, STDOUT is used (default: None)
        @type  log: LogFile
        """
        self.outFolder = T.absfile( outFolder )

        ##
        ## NOTE: If you get errors of the type "Couldn't find ID in"
        ##       check these regexps (see getSequenceIDs function)!
        ##

        self.ex_gi    = re.compile( '^>(?P<db>ref|gb|ngb|emb|dbj|prf)\|{1,2}(?P<id>[A-Z_0-9.]{4,9})\|' )
        self.ex_swiss = re.compile( '^>sp\|(?P<id>[A-Z0-9_]{5,7})\|' )
        self.ex_swiss2= re.compile( '^>.*swiss.*\|(?P<id>[A-Z0-9_]{5,7})\|' )
        self.ex_pdb   = re.compile( '^>pdb\|(?P<id>[A-Z0-9]{4})\|' )
        self.ex_all   = re.compile( '^>.*\|(?P<id>[A-Z0-9_.]{4,14})\|' )
        
        ##
        ##      RegExp for the remoteBlast searches. Somehow the
        ##      NCBI GI identifier is written before the identifier
        ##      we want (even though this is explicitly turned off in
        ##      SequenceSearcher.remoteBlast).
        ##
        
        gi = '^gi\|(?P<gi>[0-9]+)\|'
        self.ex_Rgi_1 = re.compile( gi+'(?P<db>ref|sp|pdb|gb|ngb|emb|dbj|prf|pir)\|(?P<id>[A-Z_0-9.]+)\|' )
        self.ex_Rgi_2 = re.compile( gi+'(?P<db>ref|sp|pdb|gb|ngb|emb|dbj|prf|pir)\|{2}(?P<id>[A-Z_0-9.]+)' )
        self.ex_Rswiss = re.compile( gi+'sp\|(?P<id>[A-Z0-9_]{5,7})\|' )
        self.ex_Rpdb   = re.compile( gi+'pdb\|(?P<id>[A-Z0-9]{4})\|' )

  

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
        """
        Create needed output folders if not there.
        """
        if not os.path.exists(self.outFolder + self.F_RESULT_FOLDER):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def getRecords( self ):
        """
        Get all homologues.
        
        @return: list of Bio.Fasta.Record with all found protein homologues
        @rtype: [Bio.Fasta.Record]

        @raise BlastError: if no sequences found
        """
        if not self.record_dic:
            raise BlastError( "No sequences found (yet)." )

        return self.record_dic.values()


    def getClusteredRecords( self ):
        """
        Get representative of each cluster.
        
        @return: list with Bio.Fasta.Record with best record of each cluster
        @rtype: [Bio.Fasta.Record]
        
        @raise BlastError: if called before clustering
        """
        if not self.clusters:
            raise BlastError( "Sequences are not yet clustered." )

        return [ self.record_dic[id] for id in self.bestOfCluster ]


    def __blast2dict( self, parsed_blast, db ):
        """
        Convert parsed blast result into dictionary of FastaRecords indexed
        by sequence ID.   
        Writes all the sequences in fasta format to L{F_FASTA_ALL} and the
        raw blast result to L{F_BLAST_OUT}.

        @param parsed_blast: parsed blast result
        @type  parsed_blast: Bio.Blast.Record.Blast
        @param db: database
        @type  db: str
        """
        ## fastaFromIds( db, ids ) -> { str:Bio.Fasta.Record }
        ids = self.getSequenceIDs( parsed_blast )
        self.record_dic = self.fastaFromIds( db, ids )

        if self.verbose:
            self.writeFasta( self.record_dic.values(),
                             self.outFolder + self.F_FASTA_ALL )

            self.__writeBlastResult( parsed_blast,
                                   self.outFolder + self.F_BLAST_OUT)



    def remoteBlast( self, seqFile, db, method, e=0.01, **kw ):
        """
        Perform a remote BLAST search using the QBLAST server at NCBI.
        Uses Bio.Blast.NCBIWWW.qblast (Biopython) for the search
        
        @param seqFile: file name with search sequence as FASTA
        @type  seqFile: str
        @param db: database(s) to search in, e.g. ['swissprot', 'pdb']
        @type  db: [str]
        @param method: search method, e.g. 'blastp', 'fasta'
        @type  method: str
        @param e: expectation value cutoff
        @type  e: float
        @param kw: optional keywords::
               program        BLASTP, BLASTN, BLASTX, TBLASTN, or TBLASTX.
               database       Which database to search against.
               sequence       The sequence to search.
               ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
                              (default: FALSE)
               descriptions   Number of descriptions to show.  Def 500.
               alignments     Number of alignments to show.  Def 500.
               expect         An expect value cutoff.  Def 10.0.
               matrix         Specify an alt. matrix
                              (PAM30, PAM70, BLOSUM80, BLOSUM45).
               filter         'none' turns off filtering. Default uses 'seg'
                              or 'dust'.
               format_type    'HTML', 'Text', 'ASN.1', or 'XML'.  Def. 'HTML
        @type  kw: any
        

        @note: Using the remoteBlast is asking for trouble, as
               every change in the output file might kill the
               parser. If you still want to use remoteBlast we
               strongly recomend that you install BioPython
               from CVS. Information on how to do this you
               will find on the BioPython homepage.

        @todo: Remote Blasting is running as expected but sequences
               are still retrieved from a local database. Implement remote
               collection of fasta seuqences from NCBI (there should be
               something like in Biopython). Otherwise something like this
               will also work::

                 ## collect the entry with gi 87047648
                 url = 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?CMD=text&db=protein&dopt=FASTA&dispmax=20&uid=87047648'
                 import urllib
                 handle = urllib.urlopen(url)
                 lines = handle.readlines()
                 seq = ''
                 for i in range(len(lines)):
                     if re.match( "[<]pre[>][<]div class='recordbody'[>]{2}gi", l[i] ):
                         i+= 1
                         while not re.match( "^[<]/pre[>][<]/div[>]", l[i]):
                             seq += l[i][:-1]
                             i+= 1
                 print seq       
        """
        fasta = Fasta.Iterator( open(seqFile) )
        query = fasta.next()

        blast_result = NCBIWWW.qblast( program=method, database=db,
                                       sequence=query, expect=e,
                                       ncbi_gi='FALSE', **kw )

        p = NCBIXML.BlastParser()

        parsed = p.parse( blast_result )

##        self.__remote_blast2dict( parsed, db )
        self.__blast2dict( parsed, db )

    def localBlast( self, seqFile, db, method='blastp',
                    resultOut=None, e=0.01, **kw ):
        """
        Performa a local blast search (requires that the blast binaries
        and databases are installed localy).
        Uses Bio.Blast.NCBIStandalone.blastall (Biopython) for the search.
        
        @param seqFile: file name with search sequence in FASTA format
        @type  seqFile: str
        @param db: database(s) to search, e.g. ['swissprot', 'pdb']
        @type  db: [str]
        @param method: search program to use, e.g. 'blastp', 'fasta'
                       (default: blastp)
        @type  method: str
        @param e: expectation value cutoff
        @type  e: float
        @param resultOut: save blast output to this new file
        @type  resultOut: str
        @param kw: optional keywords::
                --- Scoring ---
                matrix         Matrix to use (default BLOSUM62).
                gap_open       Gap open penalty (default 0).
                gap_extend     Gap extension penalty (default 0).

                --- Algorithm ---
                gapped         Whether to do a gapped alignment. T/F 
                                (default T)
                wordsize       Word size (blastp default 11).
                keep_hits      Number of best hits from a region to keep
                                (default off).
                xdrop          Dropoff value (bits) for gapped alignments
                                (blastp default 25).
                hit_extend     Threshold for extending hits (blastp default 11)

                --- Processing ---
                filter         Filter query sequence? (T/F, default F)
                restrict_gi    Restrict search to these GI's.
                believe_query  Believe the query defline? (T/F, default F)
                nprocessors    Number of processors to use (default 1).

                --- Formatting ---
                alignments     Number of alignments. (default 250)
        @type  kw: any
        
        @raise BlastError: if program call failes
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

##             print "Kicking out"
##             raise Exception, 'kick out to recover full result file'
        
            parsed = p.parse( results )

            self.__blast2dict( parsed, db )

        except Exception, why:
            self.error = { 'results':results,
##                           'err':err.read(),
                           'parser':p,
                           'blast_bin':settings.blast_bin, 'seqFile':seqFile,
                           'method':method, 'db':db,'kw':kw,
                           'exception':str(why) }
            print T.lastErrorTrace()
            raise BlastError( str(why) + "\n" +
                              str( self.error ) ) 


    def localPSIBlast( self, seqFile, db, resultOut=None, e=0.01, **kw ):
        """
        Performa a local psi-blast search (requires that the blast binaries
        and databases are installed localy).
        Uses Bio.Blast.NCBIStandalone.blastpgp (Biopython) for the search
        
        @param seqFile: file name with search sequence in FASTA format
        @type  seqFile: str
        @param db: database(s) to search e.g. ['swissprot', 'pdb']
        @type  db: [str]
        @param e: expectation value cutoff (default: 0.001)
        @type  e: float
        @param resultOut: save blast output to this new file
        @type  resultOut: str

        @param kw: optional keywords::
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
            keep_hits        Number of beset hits from a region to keep (def 0)
            xdrop            Dropoff value (bits) for gapped alignments
                             (def 15)
            hit_extend       Threshold for extending hits (default 11).
            nbits_gapping    Number of bits to trigger gapping (default 22).
            pseudocounts     Pseudocounts constants for multiple passes
                             (def 9).
            xdrop_final      X dropoff for final gapped alignment (default 25).
            xdrop_extension  Dropoff for blast extensions (default 7).
            model_threshold  E-value threshold to include in multipass model
                             (default 0.005).
            required_start   Start of required region in query (default 1).
            required_end     End of required region in query (default -1).

            --- Processing --- 
            filter           Filter query sequence with SEG? (T/F, default F)
            believe_query    Believe the query defline? (T/F, default F)
            nprocessors      Number of processors to use (default 1).

            --- Formatting --- 
            alignments       Number of alignments (default 250).
        @type  kw: any

        @raise BlastError: if program call failes
        """
        results = err = p = None
        try:
            results, err = NCBIStandalone.blastpgp( settings.psi_blast_bin,
                                                    db, seqFile,
                                                    expectation=e, **kw)

            p = NCBIStandalone.PSIBlastParser()

            parsed = p.parse( results )

            self.__blast2dict( parsed.rounds[-1], db )

        except Exception, why:
            self.error = { 'results':results,'err':err.read(), 'parser':p,
                           'blast_bin':settings.blast_bin, 'seqFile':seqFile,
                           'db':db,'kw':kw,'exception':str(why) }
            T.errWriteln( "BlastError " + T.lastErrorTrace() )
            print T.lastErrorTrace()
            raise BlastError( str(why) + "\n" +
                              str( self.error ) )
        


    def getSequenceIDs( self, blast_records ):
        """
        Extract sequence ids from BlastParser result.

        Use::
           getSequenceIDs( Bio.Blast.Record.Blast ) -> [ str ]
        
        @param blast_records: blast search result
        @type  blast_records: Bio.Blast.Record.Blast

        @return: list of sequence IDs
        @rtype: [str]
        
        @raise BlastError: if can't find ID
        """
        result = []

        for a in blast_records.alignments:
            for pattern in [self.ex_gi,  self.ex_swiss, self.ex_swiss2,
                            self.ex_pdb, self.ex_all,
                            self.ex_Rgi_1, self.ex_Rgi_2,
                            self.ex_Rswiss, self.ex_Rpdb]:

                match = re.search( pattern, a.title )

                if match:
                    break

            if (not match) or (not match.group('id')):
                raise BlastError( "Couldn't find ID in " + a.title +\
                                  "with pattern " + str(pattern))

            result += [ match.group('id') ]

        return result


    def fastaRecordFromId( self, db, id ):
        """
        Use::
           fastaRecordFromId( db, id ) -> Bio.Fasta.Record

        @param db: database
        @type  db: str
        @param id: sequence database ID
        @type  id: str

        @return: fasta record
        @rtype: Bio.Fasta.Record
        
        @raise BlastError: if can't fetch fasta record from database
        """
        cmd = settings.fastacmd_bin + ' -d %s -s %s' % (db, id)

        err, o = commands.getstatusoutput( cmd )
        if err:
            EHandler.warning('%s returned error: %r' % (cmd, err) )
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
        Use::
           fastaFromIds( id_lst, fastaOut ) -> { str: Bio.Fasta.Record }

        @param db: database
        @type  db: str
        @param id_lst: sequence database IDs
        @type  id_lst: [str]
        
        @return: dictionary mapping IDs to Bio.Fasta.Records
        @rtype: {str: Bio.Fasta.Record}

        @raise BlastError: if couldn't fetch record
        """
        result = {}
        for i in id_lst:
            try:
                r = self.fastaRecordFromId( db, i )
                result[i] = r
            except BlastError, why:
                EHandler.warning("couldn't fetch %s"%str(i),trace=0 )

        return result


    def copyClusterOut( self, raw=None):
        """
        Write clustering results to file.

        @param raw: write raw clustering result to disk (default: None)
        @type  raw: 1|0        
        """
        if self.verbose and raw:
            f = open( self.outFolder + self.F_CLUSTER_RAW, 'w', 1)

            f.write(raw)
            f.close()


    def reportClustering( self, raw=None ):
        """
        Report the clustering result.
        
        Writes:
         - clustering results to L{F_CLUSTER_LOG}
         - blast records to L{F_BLAST_OUT}
         - blast records of centers to L{F_CLUSTER_BLAST_OUT}
         - raw clustering results to L{F_CLUSTER_RAW} if raw not None

        @param raw: write raw clustering result to disk (default: None)
        @type  raw: 1|0         
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
        Cluster sequences. The input fasta titles must be the IDs.
        fastaClust( fastaIn [, simCut, lenCut, ncpu] )

        @param fastaIn: name of input fasta file
        @type  fastaIn: str
        @param simCut: similarity threshold (score < 3 or %identity)
                       (default: 1.75)
        @type  simCut: double
        @param lenCut: length threshold (default: 0.9)
        @type  lenCut: double
        @param ncpu: number of CPUs
        @type  ncpu: int
        
        @raise BlastError: if fastaIn is empty
        """
        fastaIn = fastaIn or self.outFolder + self.F_FASTA_ALL

        if T.fileLength( fastaIn ) < 1:
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

        @note: Older versions of T-Coffee can't handle more that approx.
               50 sequences (out of memory). This problem is taken care of
               in versions > 3.2 of T-Coffee.

        @param fastaIn: name of input fasta file
        @type  fastaIn: str
        @param simCut: similarity threshold (score < 3 or %identity)
                       (default: 1.75)
        @type  simCut: double
        @param lenCut: length threshold (default: 0.9)
        @type  lenCut: double
        @param ncpu: number of CPUs
        @type  ncpu: int  
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
        Select one member of cluster of sequences.
        
        @param ids_from_cluster: list of sequence ids defining the cluster
        @type  ids_from_cluster: [str]
        @return: id of best in cluster
        @rtype: str
        """
        return ids_from_cluster[0]


    def writeFasta( self, frecords, fastaOut ):
        """
        Create fasta file for given set of records.
        
        @param frecords: list of Bio.Blast.Records
        @type  frecords: [Bio.Blast.Record]
        @param fastaOut: file name
        @type  fastaOut: str
        """
        f = open( T.absfile(fastaOut), 'w' )
        for r in frecords:
            f.write( str( r ) + '\n' )
        f.close()


    def writeFastaAll( self, fastaOut=None ):
        """
        Write all found template sequences to fasta file.
        
        @param fastaOut: write all fasta records to file
                         (default: L{F_FASTA_ALL})
        @type  fastaOut: str OR None
        """
        fastaOut = fastaOut or self.outFolder + self.F_FASTA_ALL
        self.writeFasta( self.frecords, fastaOut )


    def writeFastaClustered( self, fastaOut=None ):
        """
        Write non-redundant set of template sequences to fasta file.
        
        @param fastaOut: write non-redundant fasta records to file
                         (default: L{F_FASTA_NR})
        @type  fastaOut: str
        """
        fastaOut = fastaOut or self.outFolder + self.F_FASTA_NR

        self.writeFasta( self.getClusteredRecords(), fastaOut )


    def __writeBlastResult( self, parsed_blast, outFile):
        """
        Write the result from the blast search to file (similar to the
        output produced by a regular blast run).
        
        writeBlastResult( parsed_blast, outFile )
        
        @param parsed_blast: Bio.Blast.Record.Blast
        @type  parsed_blast: Bio.Blast.Record.Blast
        @param outFile: file to write the blast result to
        @type  outFile: str
        """
        f = open( T.absfile( outFile ), 'w' )

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
        
        @param allFile: all blast results
        @type  allFile: file
        @param clustFile: output file name
        @type  clustFile: str
        @param selection: write only sequences in list
        @type  selection: list
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




#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0, flavour='blastp' ):
        """
        run function test

        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        @param flavour: flavour of blast to test, blastp for localBlast
                          OR blastpgp for localPSIBlast (default: blastp)
        @type  flavour: str

        @return: 1
        @rtype:  int
        """
        import tempfile
        import shutil
        from Biskit.LogFile import LogFile
        
        query = T.testRoot() + '/Mod/project/target.fasta'
        outfolder = tempfile.mkdtemp( '_test_SequenceSearcher' )
        shutil.copy( query, outfolder )

        ## log file
        f_out = outfolder+'/SequenceSearcher.log'
        l = LogFile( f_out, mode='w')

        searcher = SequenceSearcher( outFolder=outfolder, verbose=1, log=l )

        ## blast db has to be built with -o potion
        ## e.g. cat db1.dat db2.dat | formatdb -i stdin -o T -n indexed_db
        db = settings.db_swiss

        f_target = searcher.outFolder + searcher.F_FASTA_TARGET

        if local:
            globals().update( locals() )

        ## skipp remote blast for now
        # self.searcher.remoteBlast( f_target, 'nr', 'blastp', alignments=50 )

        if flavour == 'blastpgp':
            searcher.localPSIBlast( f_target, db, 'blastpgp', npasses=2 )
        else:
            searcher.localBlast( f_target, db, 'blastp', alignments=500, e=0.01 )

        searcher.clusterFasta()  ## expects all.fasta

        searcher.writeFastaClustered()

        if local:
            print '\nThe clustered result from the search can be found in %s'%f_out
            print '\nSequenceSearche log file written to: %s'%f_out
            globals().update( locals() )

        ## cleanup
        T.tryRemove( outfolder, tree=1 )
        
        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    
        

if __name__ == '__main__':

    test = Test()
    
    ## test localBlast
    assert test.run( local=1 ) ==  test.expected_result()
    
    ## test localPSIBlast
    assert test.run( local=1, flavour='blastpgp') ==  test.expected_result()


