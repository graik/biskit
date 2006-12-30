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
Search for templates.
"""

from SequenceSearcher import SequenceSearcher, BlastError
import settings
import Biskit.tools as T
from Biskit import StdLog, EHandler

from Bio import Fasta
import re
import os, shutil
import Numeric

import urllib
import string
import gzip
import subprocess
import types
from Bio import File

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

            ids = self.ex_pdb.findall( a.title )

            if not ids:
                raise BlastError( "Couldn't find ID in " + a.title)

            for id in ids:
                result += [ { 'pdb':id[0], 'chain':id[1] } ]

        return result


    def fastaFromIds( self, db, id_lst, fastaOut=None ):
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
        for i in id_lst:
            try:
                r = self.fastaRecordFromId( db, i['pdb'] )
                r.chain = i['chain']
                result[ i['pdb'] ] = r
            except BlastError, why:
                T.errWriteln("ERROR (ignored): couldn't fetch "+ str(i) )

        return result


    def getLocalPDB( self, id, db_path=settings.pdb_path ):
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
        filenames = ['%s.pdb' % id,
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


    def getRemotePDB( self, id, rcsb_url=settings.rcsb_url ):
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


    def __extractPDBInfos( self, handle ):
        """
        Extract extra infos from PDB file.
        NMR files get resolution 3.5.
        
        @param handle: open file handle OR string of file to examine
        @type  handle: open file handle OR strings
        
        @return: pdb file as list of strings, dictionary with resolution
        @rtype: [str], {'resolution':float }

        @raise BlastError: if couldn't extract PDB Info
        """
        infos = {}
        if type( handle ) == types.FileType:
            lines = handle.readlines()
        elif type( handle ) == types.StringType and len(handle) > 5000:
            lines = handle.splitlines( True )
        elif type( handle ) == types.InstanceType:
            lines = handle.readlines()
        else:
            raise BlastError( "Couldn't extract PDB Info." )

        for l in lines:
            found = self.ex_resolution.findall( l )

            if found:
                if found[0] == 'NOT APPLICABLE':
                    infos['resolution'] = self.NMR_RESOLUTION
                else:
                    infos['resolution'] = float( found[0] )
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
            T.flushPrint("retrieving %i PDBs..." % len( pdbCodes ) )
        for c in pdbCodes:

            i += 1
            if i%10 == 0 and not self.silent:
                T.flushPrint('#')

            fname = '%s/%s.pdb' % (outFolder, c)

            try:
                if os.path.exists( fname ):
                    h = open( fname, 'r' )
                else:
                    h = self.getLocalPDB( c )
            except:
                h = self.getRemotePDB( c )

            try:
                lines, infos = self.__extractPDBInfos( h )
                infos['file'] = fname

                if c in self.record_dic:
                    self.record_dic[ c ].__dict__.update( infos )

                ## close if it is a handle
                try:
                    h.close()
                except:
                    pass

                if not os.path.exists( fname ):
                    f = open( fname, 'w', 1 )
                    f.writelines( lines )
                    f.close()

                result += [ fname ]

            except IOError, why:
                raise BlastError( "Can't write file "+fname )

        if not self.silent:
            T.flushPrint('\n%i files written to %s\n' %(i,outFolder) )

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

        return ids_in_cluster[ Numeric.argmin( resolutions ) ]


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
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype:  int
        """
        import tempfile
        import shutil
        from Biskit.LogFile import LogFile
        
        query = T.testRoot() + '/Mod/project/target.fasta'
        outfolder = tempfile.mkdtemp( '_test_TemplateSearcher' )
        shutil.copy( query, outfolder )

        ## log file
        f_out = outfolder + '/TemplateSearcher.log'
        l = LogFile( f_out, mode='w')

        f_target = outfolder + '/target.fasta'

        silent = 1
        if local: silent=0

        searcher = TemplateSearcher( outFolder=outfolder,
                                     verbose=1, log=l,
                                     silent=silent )

        db = settings.db_pdbaa
        searcher.localBlast( f_target, db, 'blastp', alignments=200, e=0.0001)

        ## first tries to collect the pdb files from a local db and if
        ## that fails it tries to collect them remotely
        searcher.retrievePDBs()

        searcher.clusterFasta()  ## expects all.fasta

        searcher.writeFastaClustered()

        fn = searcher.saveClustered()
        
        if local:
            print '\nThe set of clustered template files from the search'
            print '    can be found in %s/templates'%outfolder
            print '\nTemplateSearcher log file written to: %s'%f_out
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
    
    assert test.run( local=1 ) ==  test.expected_result()



