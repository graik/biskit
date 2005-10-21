##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$

from SequenceSearcher import SequenceSearcher, BlastError
import settings
import Biskit.tools as tools

from Bio import Fasta
import re
import os, shutil
import Numeric

import urllib
import string
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
    F_CLUSTER_LOG = F_RESULT_FOLDER + '/cluster_result.out'

    ## pseudo blast output
    F_BLAST_OUT = F_RESULT_FOLDER + '/blast.out'

    ## folders for PDB files
    F_ALL       = F_RESULT_FOLDER + '/all'## all PDB homologues
    F_NR        = F_RESULT_FOLDER + '/nr' ## best PDB homologue of each cluster

    ## default name for file : chain dic 
    F_CHAIN_INDEX = '/chain_index.txt'
    

    def __init__( self, outFolder='.', verbose=1 ):
        """
        outFolder - str, project folder (results are put into subfolder) ['.']
        """
        SequenceSearcher.__init__( self, outFolder, verbose=verbose )

        self.ex_pdb   = re.compile( 'pdb\|([A-Z0-9]{4})\|([A-Z]*)' )

        self.ex_resolution = re.compile(\
             'REMARK   2 RESOLUTION\. *([0-9\.]+|NOT APPLICABLE)' )

        self.prepareFolders()


    def prepareFolders( self ):
        """
        """
        SequenceSearcher.prepareFolders( self )
        
        if not os.path.exists( self.outFolder + self.F_ALL ):
            os.mkdir( self.outFolder + self.F_ALL )
        if not os.path.exists( self.outFolder + self.F_NR ):
            os.mkdir( self.outFolder + self.F_NR )


    def getSequenceIDs( self, blast_records ):
        """
        getSequenceIDs( Bio.Blast.Record.Blast )
        -> [ {'pdb':str, 'chain':str } ]

        Extract sequence ids from BlastParser result.
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
        fastaFromIds( id_lst, fastaOut ) -> { str: Bio.Fasta.Record }
        db      - str
        id_lst  - [ {'pdb':str, 'chain':str} ]
        The returned records have additional fields: annotation and chain.
        """
        result = {}
        for i in id_lst:
            try:
                r = self.fastaRecordFromId( db, i['pdb'] )
                r.chain = i['chain']
                result[ i['pdb'] ] = r
            except BlastError, why:
                tools.errWriteln("ERROR (ignored): couldn't fetch "+ str(i) )

        return result


    def getLocalPDB( self, id, db_path=settings.pdb_path ):
        """
        id - str, pdb code, 4 characters
        db_path - str, path to pdb database
        -> File handle
        """
        id = string.lower( id )
        filenames = ['%s.pdb' % id,
                    db_path + '/pdb%s.ent' % id ]

        for f in filenames:
            if os.path.exists( f ):
                return open(f)

        raise BlastError( "Couldn't find PDB file.")


    def getRemotePDB( self, id, rcsb_url=settings.rcsb_url ):
        """
        id - str, pdb code, 4 characters
        rcsb_url - str, template url for pdb download
        -> File handle
        """
        handle = urllib.urlopen( rcsb_url% (id,id) )
        
        uhandle = File.UndoHandle(handle)
        
        if not uhandle.peekline():
            raise BlastError( "Couldn't retrieve ", rcsb_url )
    
        return uhandle


    def __extractPDBInfos( self, handle ):
        """
        Extract extra infos from PDB file line. NMR files get resolution 3.5.
        handle - open file handle
        infos  - dict getting additional infos, if any
        -> [ str ], { 'resolution':float }
        """
        infos = {}
        
        lines = handle.readlines()
        for l in lines:
            found = self.ex_resolution.findall( l )

            if found:
                if found[0] == 'NOT APPLICABLE':
                    infos['resolution'] = self.NMR_RESOLUTION
                else:
                    infos['resolution'] = float( found[0] )

        return lines, infos


    def retrievePDBs( self,  outFolder=None, pdbCodes=None,
                      db_path=settings.pdb_path,
                      rcsb_url=settings.rcsb_url ):
        """
        Write PDBs for given fasta records. Add PDB infos to internal
        dictionary of fasta records. NMR structures get resolution 3.5.
        outFolder - str, folder to put PDB files into [templates/all]
        pdbCodes  - [ str ], list of PDB codes [all previously found templates]
        -> [ str ], list of PDB file names
        """
        outFolder = outFolder or self.outFolder + self.F_ALL
        pdbCodes = pdbCodes or self.record_dic.keys()
        result = []
        i = 0
        tools.flushPrint("retrieving %i PDBs..." % len( pdbCodes ) )
        for c in pdbCodes:

            i += 1
            if i%10 == 0:
                tools.flushPrint('#')

            fname = '%s/%s.pdb' % (outFolder, c)
            try:
                if os.path.exists( fname ):
                    h = open( fname, 'r' )
                else:
                    h = self.getLocalPDB( c, db_path )
            except:
                h = self.getRemotePDB( c, rcsb_url )

            try:
                lines, infos = self.__extractPDBInfos( h )
                infos['file'] = fname

                if c in self.record_dic:
                    self.record_dic[ c ].__dict__.update( infos )

                h.close()

                if not os.path.exists( fname ):
                    f = open( fname, 'w', 1 )
                    f.writelines( lines )
                    f.close()
                
                result += [ fname ]
                
            except IOError, why:
                raise BlastError( "Can't write file "+fname )

        tools.flushPrint('\n%i files written to %s\n' %(i,outFolder) )

        return result


    def selectFasta( self, ids_in_cluster ):
        """
        select one member of cluster of sequences.
        ids - [ str ], list of sequence ids defining the cluster
        -> Bio.Fasta.Record
        """
        resolutions = []
        for id in ids_in_cluster:
            resolutions += [ getattr( self.record_dic[id],'resolution', 99. ) ]

        return ids_in_cluster[ Numeric.argmin( resolutions ) ]


    def reportClustering( self, raw=None ):
        """
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

                self.copyClusterOut( raw=raw )
                
        except IOError, why:
            tools.errWriteln( "Can't write cluster report." + str(why) )


    def saveClustered( self, outFolder=None ):
        """
        Copy best PDB of each cluster into another folder.
        Create index file in same folder.
        The returned dictionary or  index file can be used as input to
        TemplateCleaner.
        -> { str_filename : str_chain_id }, file names and chain ids
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

    db = 'pdbaa'
    outfolder = tools.projectRoot()+ '/test/Mod/project'
    f_target = outfolder + '/target.fasta'

    searcher = TemplateSearcher(outFolder=outfolder,
                                verbose=1 )

    searcher.localBlast( f_target, db, 'blastp', alignments=200, e=0.0001)

    searcher.retrievePDBs()

    searcher.clusterFasta()  ## expects all.fasta

    searcher.writeFastaClustered()

    fn = searcher.saveClustered()


