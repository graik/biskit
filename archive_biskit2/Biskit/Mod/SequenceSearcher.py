##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
## Contributions: Marc-Michael Blum (transition to Bio.SeqIO)
##
"""
Search one or more sequence databases (SwissProt, Tremble) with a sequence
using the method of choice (Blast , FastA)
"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import Biskit.tools as T
import Biskit.Errors as E
from Biskit import StdLog, EHandler

from Biskit.Mod import settings

import re, os
import string
import commands
import copy
import cStringIO

class BlastError( E.BiskitError ):
    pass

class InternalError( E.BiskitError ):
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
    F_BLAST_OUT   = F_RESULT_FOLDER + '/blast.out'
    F_CLUSTER_BLAST_OUT = F_RESULT_FOLDER + '/cluster_blast.out'

    F_BLAST_RAW_OUT =  F_RESULT_FOLDER + '/blast_raw.out'
    F_BLAST_ERROR = F_RESULT_FOLDER + '/blast.error'

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

        self.verbose = verbose
        self.log = log or StdLog()

        self.record_dic = None #: dict ID - Bio.Fasta.Record
        self.clusters = None #: list of lists of ids
        self.bestOfCluster = None #: list of non-redundant ids

        #: the maximal number of clusters to return
        self.clusterLimit =  clusterLimit

        self.clustersCurrent = None #: current number of clusters

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


    def __blast2dict( self, parsed_blast, db, remote=False ):
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
        self.record_dic = self.fastaFromIds( db, ids, remote=remote )

        self.writeFasta( self.record_dic.values(),
                         self.outFolder + self.F_FASTA_ALL )

        self.__writeBlastResult( parsed_blast,
                                 self.outFolder + self.F_BLAST_OUT)


    def __copyFileHandle(self, result_handle, fname ):
        """
        Copy the blast result to a file. The input file handle will be
        empty afterwards! Use the returned string handle instead.
    
        @param result_handle: file handle returned by NCBI* blast runners
        @type  result_handle: file
        @param fname : absolute path to output file
        @type  fname : str
    
        @return: a fresh 'file' (string) handle with the output
        @rtype : cStringIO.StringI
        """
        str_results = result_handle.read()
        try:
            save_file = open( fname, 'w' )
            save_file.write( str_results )
            save_file.close()
        finally:
            return cStringIO.StringIO( str_results )

    def __dictvalues2str( self, kw ):
        """convert all dict values to strings; see bug #2019987"""
        r = {}
        for k,v in kw.items():
            r[k] = str( v )
        return r

    def remoteBlast( self, seqFile, db, method, e='0.01', **kw ):
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
               from CVS. Information on how to do this can be
               found on the BioPython homepage.

        @todo: Remote Blasting is running as expected but sequences
               are still retrieved from a local database. Implement remote
               collection of fasta seuqences from NCBI (there should be
               something alike in Biopython). Otherwise something like this
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
        kw = self.__dictvalues2str( kw )
        e = str(e)

        try:
            fasta = SeqIO.parse( open(seqFile), "fasta" )
            query = fasta.next()

            if self.verbose:
                self.log.add('starting blast query...')

            results = NCBIWWW.qblast( program=method, database=db,
                                      sequence=query, expect=e,
                                      ncbi_gi='FALSE', **kw )

            results = self.__copyFileHandle(results,
                                            self.outFolder+self.F_BLAST_RAW_OUT)
            if self.verbose:
                self.log.writeln('Raw blast output copied to: ' +\
                                 self.outFolder + self.F_BLAST_RAW_OUT  )

            parsed = NCBIXML.parse( results ).next()

            self.__blast2dict( parsed, db, remote=True )

        except Exception, why:
            self.log.add( T.lastErrorTrace() )
            globals().update( locals() )
            self.log.writeln('local namespace is pushed into global ')
            raise BlastError( str(why) ) 


    def localBlast( self, seqFile, db, method='blastp',
                    resultOut=None, e='0.01', **kw ):
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
        resultOut = resultOut or self.outFolder+ self.F_BLAST_RAW_OUT
        kw = self.__dictvalues2str( kw )
        e = str(e)

        try:
            if self.verbose:
                self.log.add('running blast...')

            results, err = NCBIStandalone.blastall( settings.blast_bin,
                                                    method, db, seqFile,
                                                    expectation=e,
                                                    align_view='7', ## XML output
                                                    **kw)

            results = self.__copyFileHandle(results, resultOut)
            err = self.__copyFileHandle(err, self.outFolder+self.F_BLAST_ERROR)

            if self.verbose:
                self.log.writeln('Raw blast output copied to: ' + resultOut  )

            parsed = NCBIXML.parse( results ).next()

            self.__blast2dict( parsed, db )

        except Exception, why:
            self.log.add( T.lastErrorTrace() )
            globals().update( locals() )
            self.log.writeln('local namespace is pushed into global ')
            raise BlastError( str(why) ) 


    def localPSIBlast( self, seqFile, db, method='blastp',
                       resultOut=None, e='0.001', **kw ):
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
            --- New Blast+ routine ---
            (see NcbipsiblastCommandline)

            num_iterations   Number of passes (default 1).
            matrix           Matrix to use (default BLOSUM62).
            
            --- old blastall routine ---
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
        ## the following should work for new Blast+ tools:
        
        #from Bio.Blast.Applications import NcbipsiblastCommandline

        #resultOut = resultOut or self.outFolder+ self.F_BLAST_RAW_OUT
        #blastx_cline = NcbipsiblastCommandline(query=seqFile, 
                                               #db=db, 
                                               #evalue=e,
                                               #outfmt=5, 
                                               #out=resultOut,
                                               #**kw)
        #stdout, stderr = blastx_cline()
        #parsed = NCBIXML.parse( results ).next()
        #self.__blast2dict( parsed, db )
        
        results = err = None
        resultOut = resultOut or self.outFolder+ self.F_BLAST_RAW_OUT
        kw = self.__dictvalues2str( kw )
        e = str(e)

        try:
            results, err = NCBIStandalone.blastpgp( settings.psi_blast_bin,
                                                    db, seqFile,
                                                    program='blastpgp',
                                                    align_view='7', ## XML output
                                                    expectation=e, **kw)

            results = self.__copyFileHandle(results,resultOut )
            err = self.__copyFileHandle(err, self.outFolder+self.F_BLAST_ERROR)

            if self.verbose:
                self.log.writeln('Raw blast output copied to: ' + resultOut )

            parsed = NCBIXML.parse( results ).next()

            self.__blast2dict( parsed, db )

        except Exception, why:
            self.log.add( T.lastErrorTrace() )
            globals().update( locals() )
            self.log.writeln('local namespace is pushed into global ')
            raise BlastError( str(why) ) 


    def getSequenceIDs( self, blast_records ):
        """
        Extract sequence ids from BlastParser result.

        Use::
           getSequenceIDs( Bio.Blast.Record.Blast ) -> [ str ]

        @param blast_records: blast search result
        @type  blast_records: Bio.Blast.Record.Blast

        @return: list of sequence IDs
        @rtype: [str]
        """
        result = [ a.accession for a in blast_records.alignments ]

        return result
    
    def fastaRecordFromId_remote( self, id ):
        """
        experimental: fetch fasta records from remote database
        """
        from Bio import Entrez
        from Bio import SeqIO
        Entrez.email = "A.N.Other@example.com"
        if self.verbose: self.log.add_nobreak( 'r' )
        handle = Entrez.efetch(db="protein", rettype="fasta", id=id)
        frecord = SeqIO.read(handle, "fasta")
        frecord.id = str(id)
        handle.close()
        return frecord

    def fastaRecordFromId( self, db, id ):
        """
        Use::
           fastaRecordFromId( db, id ) -> Bio.Fasta.Record

        @param db: database
        @type  db: str
        @param id: sequence database ID
        @type  id: str

        @return: fasta record
        @rtype: Bio.SeqRecord.SeqRecord

        @raise BlastError: if can't fetch fasta record from database
        """
        cmd = settings.fastacmd_bin + ' -d %s -s %s' % (db, id)

        err, o = commands.getstatusoutput( cmd )
        if err:
            EHandler.warning('%s returned error: %r' % (cmd, err) )
            raise BlastError( 'fastacmd failed. Error code: ' + str(err) )

        try:
            frecord = SeqIO.parse( cStringIO.StringIO(o), 'fasta').next()
            frecord.id = str(id)

        except StopIteration:
            raise InternalError, \
                  "Couldn't fetch fasta record %s from database %s" % (id,db)

        return frecord


    def fastaFromIds( self, db, id_lst, remote=False ):
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
        if self.verbose:
            s = 'from local %s using fastacmd'
            if remote:
                s = 'remotely from Entrez'
            self.log.add('Fetching %i fasta records %s...\n'% (len(id_lst), s))
            
        for i in id_lst:
            try:
                if remote:
                    r = self.fastaRecordFromId_remote( i )
                else:
                    r = self.fastaRecordFromId( db, i )
                result[i] = r
            except BlastError, why:
                EHandler.warning("couldn't fetch %s"%str(i),trace=0 )
                
        if self.verbose:
            self.log.add('done.\n')

        return result


    def copyClusterOut( self, raw=None):
        """
        Write clustering results to file.

        @param raw: write raw clustering result to disk (default: None)
        @type  raw: str       
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

        if self.verbose:
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
            raise BlastError( "blastclust failed. Error code: " + str(err) )

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
        while (self.clustersCurrent > self.clusterLimit \
              or self.clustersCurrent is None) \
              and (simCut > 0 and lenCut > 0):

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
            f.write( r.format('fasta') )  ## note better use direct SeqIO
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
        try:
            f = open( T.absfile( outFile ), 'w' )

            i=1
            for alignment in parsed_blast.alignments:
                for hsp in alignment.hsps:
                    s = string.replace(alignment.title,'\n',' ')
                    s = string.replace(s, 'pdb|',  '\npdb|')
                    f.write('Sequence %i: %s\n'%(i,s))                
                    f.write('Score: %3.1f \tE-value: %2.1e\n'\
                            %(hsp.score, hsp.expect))
                    f.write('Lenght/Identities: %r\tPositives: %r\tGaps: %r\n'\
                            %(hsp.identities, hsp.positives, hsp.gaps))

                    f.write( '%s\n'%hsp.query  )
                    f.write( '%s\n'%hsp.match )
                    f.write( '%s\n\n'%hsp.sbjct )
                    i += 1
            f.close()
        except Exception, why:
            EHandler.warning("Error while writing blast result to %s" %outFile)
            globals().update(locals())
            EHandler.warning("function namespace published to globals")


    def writeClusteredBlastResult( self, allFile, clustFile, selection ):
        """
        Reads the blast.out file and keeps only centers.

        @param allFile: c!all blast results
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
import Biskit.test as BT

class TestBase(BT.BiskitTest):
    """Prepare all Test cases (without providing any itself)"""

    def prepare(self):
        self.query = T.testRoot() + '/Mod/project/target.fasta'

    def prepareOutput(self, ident):
        """Set up an output folder for a blast run and point
        self.outfolder and log files to the right location
        """
        import tempfile
        import shutil
        from Biskit.LogFile import LogFile

        self.outfolder = tempfile.mkdtemp( '_test_SequenceSearcher_%s'%ident )

        shutil.copy( self.query, self.outfolder )

        ## log file
        self.f_out = self.outfolder+'/SequenceSearcher.log'


    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )

    def searchLocalBlast(self):
        self.searcher.localBlast( self.f_target, self.db, 'blastp',
                                  alignments=500, e=0.01 )

    def searchLocalPsiBlast(self):
        self.searcher.localPSIBlast( self.f_target, self.db, npasses=2, 
                                     e=0.0000000001 )

    def searchRemoteBlast(self):
        self.searcher.remoteBlast( self.f_target, self.db, 'blastp',
                                   alignments=50)

    def search(self, searchFunction ):
        from Biskit import LogFile

        self.prepareOutput( searchFunction.__name__ )

        l = self.log
        if self.local:
            l = None   ## use STDOUT

        self.searcher = SequenceSearcher( outFolder=self.outfolder,
                                          verbose=self.local, log=l )

        ## blast db has to be built with -o potion
        ## e.g. cat db1.dat db2.dat | formatdb -i stdin -o T -n indexed_db
        self.db = settings.db_swiss

        self.f_target = self.searcher.outFolder + self.searcher.F_FASTA_TARGET

        ## perform search variant
        searchFunction()

        self.searcher.clusterFasta()  ## expects all.fasta

        self.searcher.writeFastaClustered()

        if self.local and self.DEBUG:
            print '\nThe clustered result from the search can be found in %s'%\
                  self.outfolder
            print '\nSequenceSearche log file written to: %s' % self.f_out

        self.assert_( len(self.searcher.clusters) > 0 )
        self.assertEqual( len(self.searcher.clusters),
                          len(self.searcher.bestOfCluster) )

        ## make sure we have a fasta record for each cluster center
        for id in self.searcher.bestOfCluster:
            self.assert_( id in self.searcher.record_dic )


class TestLocal( TestBase ):
    """Test case for local Blast and PSI-Blast"""

    TAGS = [ BT.EXE ]

    def test_blast(self):
        """Mod.SequenceSearcher.localBlast test"""
        self.search( self.searchLocalBlast )

    def skippedtest_psiblast(self):
        """Mod.SequenceSearcher.localPSIBlast test"""
        self.search( self.searchLocalPsiBlast )


class TestRemote( TestBase ):
    """Test case for remote (WWW-)Blast""" 

    TAGS = [ BT.LONG, BT.EXE ]

    def test_remoteBlast(self):
        """Mod.SequenceSearcher.remoteBlast test"""
        self.search( self.searchRemoteBlast )


if __name__ == '__main__':

    BT.localTest(debug=0)
