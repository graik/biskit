## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## last $Author$
## last $Date$
## $Revision$

"""
Search for templates.
"""

from SequenceSearcher import SequenceSearcher, BlastError
import Biskit.settings as bisSettings
import settings
import Biskit.tools as T
from Biskit import StdLog, EHandler

import re
import os, shutil
import numpy as N

import urllib
import string
import gzip
import subprocess
import  cStringIO
import commands
from Bio import File
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class PickyURLopener(urllib.FancyURLopener):
    """Override default urllib error handling -- always raise an exception."""

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        """
        Default error handling.
        @raise IOError: if something went wrong
        """
        raise IOError, 'Cannot open %r. Error %r (%s)' % (url,errcode,errmsg)

urllib._urlopener = PickyURLopener()


class TemplateSearcher( SequenceSearcher ):
    """
    Take a sequence and return a list of files of nonredundant PDB homologues
    (selecting for best resolution).
    """

    ## resolution assigned to NMR structures
    NMR_RESOLUTION = 3.5

    ## default folder names (created if verbose=1)
    F_RESULT_FOLDER = '/templates'        ## default parent folder

    ## standard file name for unclustered sequences
    F_FASTA_ALL = F_RESULT_FOLDER + '/all.fasta'

    ## standard file name for non-redundant seqs
    F_FASTA_NR =  F_RESULT_FOLDER + '/nr.fasta'

    ## clusters
    F_CLUSTER_RAW = F_RESULT_FOLDER + '/cluster_raw.out'
    F_CLUSTER_LOG = F_RESULT_FOLDER + '/cluster_result.out'

    ## pseudo blast output
    F_BLAST_OUT = F_RESULT_FOLDER + '/blast.out'
    F_CLUSTER_BLAST_OUT = F_RESULT_FOLDER + '/cluster_blast.out'

    ## folders for PDB files
    F_ALL       = F_RESULT_FOLDER + '/all'## all PDB homologues
    F_NR        = F_RESULT_FOLDER + '/nr' ## best PDB homologue of each cluster

    ## default name for file : chain dic 
    F_CHAIN_INDEX = '/chain_index.txt'


    def __init__( self, outFolder='.', clusterLimit=20, verbose=1,
                  log=None, silent=0 ):
        """
        @param outFolder: project folder (results are put into subfolder) ['.']
        @type  outFolder: str
        @param clusterLimit: maximal number of returned sequence clusters
                             (default: 20)
        @type  clusterLimit: int
        @param verbose: keep temporary files (default: 1)
        @type  verbose: 1|0
        @param log: log file instance, if None, STDOUT is used (default: None)
        @type  log: LogFile
        @param silent: don't print messages to STDOUT is used (default: 0)
        @type  silent: 1|0
        """
        SequenceSearcher.__init__( self, outFolder=outFolder, verbose=1 )

        self.ex_pdb   = re.compile( 'pdb\|([A-Z0-9]{4})\|([A-Z]*)' )
        self.ex_gb    = re.compile( 'gi\|([0-9]+)\|')

        self.ex_resolution = re.compile(\
            'REMARK   2 RESOLUTION\. *([0-9\.]+|NOT APPLICABLE)' )

        self.prepareFolders()

        self.verbose = verbose
        self.silent = silent

        self.log = log or StdLog()

        ## the maximal number of clusters to return
        self.clusterLimit = clusterLimit


    def prepareFolders( self ):
        """
        Create folders needed by this class.
        """
        SequenceSearcher.prepareFolders( self )

        if not os.path.exists( self.outFolder + self.F_ALL ):
            os.mkdir( self.outFolder + self.F_ALL )
        if not os.path.exists( self.outFolder + self.F_NR ):
            os.mkdir( self.outFolder + self.F_NR )


    def getSequenceIDs( self, blast_records ):
        """
        Extract sequence ids (pdb codes and chain ID) from BlastParser result.

        @param blast_records: result from BlastParser
        @type  blast_records: Bio.Blast.Record.Blast

        @return: list of dictionaries mapping pdb codes and chain IDs
        @rtype: [ {'pdb':str, 'chain':str } ]

        @raise BlastError: if couldn't find ID
        """
        result = []
        for a in blast_records.alignments:

            ids_pdb = self.ex_pdb.findall( a.title )
            ids_gb  = self.ex_gb.findall( a.title )

            if not ids_pdb:
                raise BlastError( "Couldn't find PDB ID in " + a.title)
            if not ids_gb:
                raise BlastError( "Couldn't find genebank ID in " + a.title)
            assert len( ids_pdb ) == len(ids_gb)

            for pdb, gb in zip( ids_pdb, ids_gb):
                result += [ { 'pdb':pdb[0], 'chain':pdb[1], 'gb':gb } ]

        return result

    def fastaRecordFromId( self, db, id, chain='' ):
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
        if chain == '':
            return SequenceSearcher.fastaRecordFromId( self, db, id )

        cmd = settings.fastacmd_bin + " -d %s -s 'pdb|%s|%s'" \
            % (db, id, chain)

        err, o = commands.getstatusoutput( cmd )
        if err:
            EHandler.warning('%s returned error: %r' % (cmd, err) )
            raise BlastError( err )

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

        @param db: database name
        @type  db: str
        @param id_lst: list of dictionaries with pdb codes and chain IDs
        @type  id_lst: [{'pdb':str, 'chain':str}]

        @return: Dictionary mapping pdb codes to Bio.Fasta.Records. The
                 returned records have an additional field: chain.
        @rtype: { str: Bio.Fasta.Record }        
        """
        result = {}
        if self.verbose:
            s = 'from local %s using fastacmd' % db
            if remote:
                s = 'remotely from Entrez'
            self.log.add('Fetching %i fasta records %s...\n'% (len(id_lst), s))

        for i in id_lst:
            try:
                if remote:
                    r = self.fastaRecordFromId_remote( i['gb'] )
                    r.id = i['pdb']  ## clustering expects PDB, not gb ID
                else:
                    r = self.fastaRecordFromId( db, i['pdb'], i['chain'] )
                r.chain = i['chain']
                result[ i['pdb'] ] = r
            except BlastError, why:
                EHandler.warning("ERROR (ignored): couldn't fetch "+ str(i) )

        return result


    def getLocalPDBHandle( self, id, db_path=bisSettings.pdb_path ):
        """
        Get the coordinate file from a local pdb database.

        @param id: pdb code, 4 characters
        @type  id: str
        @param db_path: path to local pdb database
                        (default: L{settings.pdb_path})
        @type  db_path: str

        @return: the requested pdb file as a file handle
        @rtype: open file handle

        @raise BlastError: if couldn't find PDB file
        """
        id = string.lower( id )
        filenames = [os.path.join( db_path, '%s.pdb' % id),
                     db_path + '/pdb%s.ent' % id,
                     db_path + '/%s/pdb%s.ent.Z' %( id[1:3], id ) ]

        for f in filenames:
            if os.path.exists( f ):
                ## gzipped pdb file
                if f[-3:]=='.gz':
                    return gzip.open(f)
                ## the gzip module doesn't handle .Z files
                ## doesn't return open file handle 
                elif f[-2:]=='.Z':
                    p = subprocess.Popen( [ 'gunzip', '-c', f ],
                                          stdout=subprocess.PIPE )
                    return p.communicate()[0]
                ## uncompressed
                else:
                    return open(f)

        raise BlastError( "Couldn't find PDB file.")


    def getRemotePDBHandle( self, id, rcsb_url=bisSettings.rcsb_url ):
        """
        Get the coordinate file remotely from the RCSB.

        @param id: pdb code, 4 characters
        @type  id: str
        @param rcsb_url: template url for pdb download
                         (default: L{settings.rcsb_url})
        @type  rcsb_url: str

        @return: the requested pdb file as a file handle
        @rtype: open file handle

        @raise BlastError: if couldn't retrieve PDB file
        """
        handle = urllib.urlopen( rcsb_url% (id,id) )

        uhandle = File.UndoHandle(handle)

        if not uhandle.peekline():
            raise BlastError( "Couldn't retrieve ", rcsb_url )

        return uhandle


    def parsePdbFromHandle(self, handle, first_model_only=True ):
        """
        Parse PDB from file/socket or string handle into memory.

        @param handle: fresh open file/socket handle to PDB ressource or string
        @type  handle: open file-like object or str
        @param first_model_only: only take first of many NMR models [True]
        @type  first_model_only: bool

        @return: pdb file as list of strings, dictionary with resolution
        @rtype: [str], {'resolution':float }
        @raise BlastError: if passed in string is too short
        """
        lines = []
        res_match = None
        infos = {}

        if type( handle ) is str:
            if len(handle) < 5000:
                raise BlastError( "Couldn't extract PDB Info." )
            handle =  cStringIO.StringIO( handle )

## if handle.peekline()[:6] != 'TITLE':
##     raise BlastError, 'Ressource does not seem to be a PDB:\n%r' %\
##   handle.peekline()

        for l in handle:
            lines += [ l ]

            res_match = res_match or self.ex_resolution.search( l )

            if first_model_only and l[:6] == 'ENDMDL':
                break

        if res_match:
            if res_match.groups()[0] == 'NOT APPLICABLE':
                infos['resolution'] = self.NMR_RESOLUTION
            else:
                infos['resolution'] = float( res_match.groups()[0] )
        else:
            raise BlastError, 'No resolution record found in PDB.'

        return lines, infos


    def retrievePDBs( self, outFolder=None, pdbCodes=None ):
        """
        Get PDB from local database if it exists, if not try to
        download the coordinartes drom the RSCB.
        Write PDBs for given fasta records. Add PDB infos to internal
        dictionary of fasta records. NMR structures get resolution 3.5.

        @param outFolder: folder to put PDB files into (default: L{F_ALL})
        @type  outFolder: str OR None
        @param pdbCodes: list of PDB codes [all previously found templates]
        @type  pdbCodes: [str]

        @return: list of PDB file names
        @rtype: [str]

        @raise BlastError: if can't write file
        """
        outFolder = outFolder or self.outFolder + self.F_ALL
        pdbCodes = pdbCodes or self.record_dic.keys()
        result = []
        i = 0
        if not self.silent:
            T.flushPrint("fetching %i PDBs (l=local, r=remotely)..." % \
                         len( pdbCodes ) )

        for c in pdbCodes:

            i += 1

            fname = '%s/%s.pdb' % (outFolder, c)

            try:
                if os.path.exists( fname ):
                    h = open( fname, 'r' )
                else:
                    h = self.getLocalPDBHandle( c )
                if not self.silent:
                    T.flushPrint('l')
            except:
                h = self.getRemotePDBHandle( c )
                if not self.silent:
                    T.flushPrint('r')

            try:
                lines, infos = self.parsePdbFromHandle( h, first_model_only=1 )
                infos['file'] = fname

                if c in self.record_dic:
                    self.record_dic[ c ].__dict__.update( infos )

                ## close if it is a handle
                try: h.close()
                except:
                    pass

                if not os.path.exists( fname ):
                    f = open( fname, 'w', 1 )
                    f.writelines( lines )
                    f.close()

                result += [ fname ]

            except IOError, why:
                raise BlastError( "Can't write file "+fname )

        if self.verbose:
            self.log.add('\n%i PDB files written to %s\n' %(i,outFolder) )

        return result


    def selectFasta( self, ids_in_cluster ):
        """
        select one member of cluster of sequences.

        @param ids_in_cluster: list of sequence ids defining the cluster
        @type  ids_in_cluster: [str]

        @return: Bio.Fasta.Record
        @rtype: Bio.Fasta.Record
        """
        resolutions = []
        for id in ids_in_cluster:
            resolutions += [ getattr( self.record_dic[id],'resolution', 99. ) ]

        return ids_in_cluster[ N.argmin( resolutions ) ]


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
                f = open( self.outFolder +self.F_CLUSTER_LOG, 'w', 1)

                for cluster in self.clusters:
                    f.write( "%i\t" % ( len( cluster )))
                    for id in cluster:
                        f.write("%s:%.2f "%(id,self.record_dic[id].resolution))
                    f.write( "\n")

                f.close()

                ## write blast records of centers to disc
                centers = [ c[0] for c in self.clusters ]

                self.writeClusteredBlastResult( \
                    self.outFolder + self.F_BLAST_OUT,
                    self.outFolder + self.F_CLUSTER_BLAST_OUT, centers )

                self.copyClusterOut( raw=raw )

        except IOError, why:
            T.errWriteln( "Can't write cluster report." + str(why) )


    def saveClustered( self, outFolder=None ):
        """
        Copy best PDB of each cluster into another folder.
        Create index file in same folder.
        The returned dictionary or index file (L{F_CHAIN_INDEX}) is
        used as input to TemplateCleaner.

        @param outFolder: folder to write files to (default: L{F_NR})
        @type  outFolder: str OR None

        @return: { str_filename : str_chain_id }, file names and chain ids
        @rtype: {str, str}

        @raise BlastError: if sequences are not clustered
        """
        result = {}
        outFolder = outFolder or self.outFolder + self.F_NR

        f_index = open( outFolder +self.F_CHAIN_INDEX, 'w')
        f_index.write('## template files mapped to matching chain ID\n')

        if not self.bestOfCluster:
            raise BlastError( 'Sequences are not yet clustered.' )

        for id in self.bestOfCluster:
            fn = self.record_dic[ id ].file
            fnNew = outFolder+'/'+id+'.pdb'
            shutil.copy( fn, fnNew )

            self.record_dic[ id ].file = fnNew
            result[ fnNew ] = self.record_dic[ id ].chain

            f_index.write( '%s\t%s\n' % (fnNew, self.record_dic[ id ].chain) )

        f_index.close()

        return result



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """

    TAGS = [ BT.EXE, BT.LONG ]

    def prepare(self):
        import tempfile
        import shutil
        from Biskit.LogFile import LogFile

        self.query = T.testRoot() + '/Mod/project/target.fasta'
        self.outfolder = tempfile.mkdtemp( '_test_TemplateSearcher' )
        shutil.copy( self.query, self.outfolder )

        ## log file
        self.f_out = self.outfolder + '/TemplateSearcher.log'

        self.l = LogFile( self.f_out, mode='w')
        if self.local:  self.l = StdLog()

        self.f_target = self.outfolder + '/target.fasta'


    def test_TemplateSearcher(self):
        """Mod.TemplateSearcher test"""
        silent = 0
        if self.local: silent=0

        self.searcher = TemplateSearcher( outFolder=self.outfolder,
                                          verbose=1, log=self.l,
                                          silent=silent )

        db = settings.db_pdbaa
        self.searcher.localBlast( self.f_target, db, 'blastp',
                                  alignments=200, e=0.0001)

        ## first tries to collect the pdb files from a local db and if
        ## that fails it tries to collect them remotely
        self.searcher.retrievePDBs()

        self.searcher.clusterFasta()  ## expects all.fasta

        self.searcher.writeFastaClustered()

        fn = self.searcher.saveClustered()

        if self.local and self.DEBUG:
            self.log.add(
                '\nThe set of clustered template files from the search' +\
                '    can be found in %s/templates' % self.outfolder +\
                '\nTemplateSearcher log file written to: %s' % self.f_out)

    def test_TemplateSearcherRemote(self):
        """Mod.TemplateSearcher remote blast test"""
        silent = 0
        if self.local: silent=0

        self.searcher = TemplateSearcher( outFolder=self.outfolder,
                                          verbose=1, log=self.l,
                                          silent=silent )

        db = settings.db_pdbaa
        self.searcher.remoteBlast( self.f_target, db, 'blastp',
                                  alignments=200, e=0.0001)

        ## first tries to collect the pdb files from a local db and if
        ## that fails it tries to collect them remotely
        self.searcher.retrievePDBs()

        self.searcher.clusterFasta()  ## expects all.fasta

        self.searcher.writeFastaClustered()

        fn = self.searcher.saveClustered()

        if self.local and self.DEBUG:
            self.log.add(
                '\nThe set of clustered template files from the search' +\
                '    can be found in %s/templates' % self.outfolder +\
                '\nTemplateSearcher log file written to: %s' % self.f_out)

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )


if __name__ == '__main__':

    BT.localTest(debug=False)
