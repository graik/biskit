##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$

import Biskit.tools as tools
from Biskit import PDBModel, PDBCleaner, CleanerError

from TemplateSearcher import TemplateSearcher
import modUtils

import os, string

class TemplateCleaner:
    """
    Takes a list of PDB-files and chains identifiers.
    -> cleaned PDB files with all chains
    -> cleaned PDB files with needed chain for modeller.
    -> CA trace PDB for T_coffee
    -> sequences in fasta format for T_coffee

    mkdir cleaned  .. complete PDBs but cleaned
          modeller .. only chains needed
          t_coffee .. CA traces
    templates.fasta
    """ 

    F_RESULT_FOLDER = TemplateSearcher.F_RESULT_FOLDER

    F_CLEANED = F_RESULT_FOLDER + '/nr_cleaned/'
    F_MODELLER= F_RESULT_FOLDER + '/modeller/'
    F_COFFEE  = F_RESULT_FOLDER + '/t_coffee/'
    F_FASTA   = F_RESULT_FOLDER + '/templates.fasta'

    def __init__( self, outFolder, log=None ):
        """
        pdbFiles - [ str ]
        chain_ids- [ str ], chain ids for each file, '' means all
        log      - LogFile instance or None, None reports to STDOUT
        """
        self.outFolder = tools.absfile( outFolder )
        self.log = log

        self.prepareFolders()


    def prepareFolders( self ):
        """
        """
        if not os.path.exists( self.outFolder + self.F_CLEANED ):
            os.mkdir( self.outFolder + self.F_CLEANED )
        if not os.path.exists( self.outFolder + self.F_MODELLER ):
            os.mkdir( self.outFolder + self.F_MODELLER )
        if not os.path.exists( self.outFolder + self.F_COFFEE ):
            os.mkdir( self.outFolder + self.F_COFFEE )


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg


    def __chain_by_id( self, model, chainId ):
        """
        chainId - str, chain identifier
        -> PDBModel
        """
        return model.compress( model.mask( \
            lambda a, id=chainId: a['chain_id'] == id ))


    def fasta_sequence( self, header, s ):

        n_chunks = len( s ) / 80

        result = ">%s\n" % header 

        for i in range(0, n_chunks+1):

            if i * 80 + 80 < len( s ):
                chunk = s[i * 80 : i * 80 + 80]
            else:
                chunk = s[i * 80 :]

            result += chunk + '\n'

        return result


    def write_tcoffee_pdb( self, model, fname ):
        """
        Write CA trace PDB with SEQRES records for t_coffee.
        model  - PDBModel
        fname  - str
        """
        ## write CA trace
        m_ca = model.compress( model.maskCA() )

        ## prepare sequence block

        n_lines = len( m_ca ) / 13
        if not len( m_ca ) % 13 == 0:
            n_lines += 1

        chain_id = m_ca.getAtoms()[0]['chain_id']

        ## structure alignment program needs chain identifier
        if not chain_id:
            for a in m_ca.atoms:
                a['chain_id'] = 'A'

        ## create SEQRES record
        n_res = m_ca.lenResidues()

        head = []

        res_from = 0
        res_to = 13
        for i in range( n_lines ):
            s = "%4i %s %4i " % (i+1, chain_id, n_res)

            if res_to > n_res:
                res_to = n_res

            for ca in m_ca.atoms[ res_from: res_to ]:
                s+= " %3s" % ca['residue_name']

            head += [ ( 'SEQRES', s ) ]

            res_from += 13
            res_to += 13

        m_ca.writePdb( fname, headlines=head )
        


    def write_modeller_pdb( self, model, fname ):
        """
        """
        model = model.clone( deepcopy=1 )
        for a in model.atoms:
            a['chain_id'] = ''
            
        model.writePdb( fname )
        

    def process_all( self, file_dic, keep_hetatoms=0 ):
        """
        Process PDB files in file_dic. The PDB is read from /templates/nr and
        cleaned files are written to /templates/nr_cleaned, /templates/t_coffee
        and /templates/modeller. If the /templates/nr_cleaned already exists,
        this file is used to write modeller and t-coffee pdbs.
        
        file_dic      - { fname : chain_id, }
        keep_hetatoms - 0||1, [0]
        """
        fasta = ""

        for f, id in file_dic.items():

            self.logWrite( '\nCLEANING ' + f + '...')

            try:
                code = string.split(f, '/')[-1][:4]

                ## clean to standard PDB
                c = PDBCleaner( f, self.log )

                ## if /templates/nr_cleaned/pfb-file exists
                if os.path.exists(
                    self.outFolder+self.F_CLEANED + code + '.pdb' ):

                    model = PDBModel( self.outFolder+self.F_CLEANED \
                                      + code + '.pdb' )

                ## read pdb from /templates/nr
                else:
                    model = c.process( keep_hetatoms=keep_hetatoms )

                    ## write complete cleaned PDB
                    model.writePdb( self.outFolder+self.F_CLEANED \
                                    + code + '.pdb')

                code = model.pdbCode
                title = code + '_' + id

                ## extract needed chain
                if len( id ) > 0:
                    model = self.__chain_by_id( model, id )

                fname =  "%s%s.pdb" % (self.outFolder + self.F_MODELLER, title)
                self.write_modeller_pdb( model, fname )

                ## write T_coffee CA trace with res sequence
                fname = "%s%s.alpha"%(self.outFolder + self.F_COFFEE, title)
                self.write_tcoffee_pdb( model, fname )

                ## append fasta
                fasta += self.fasta_sequence( title, model.sequence() )

            except:
                self.logWrite( 'Error cleaning ' + f)
                self.logWrite( tools.lastError() )
                self.err_cleaner = c
                self.err_file = f

        fasta_out = open( self.outFolder + self.F_FASTA, 'w' )
        fasta_out.write( fasta )
        fasta_out.close()


if __name__ == '__main__':

    outfolder = tools.projectRoot()+ '/test/Mod/project'

    c = TemplateCleaner( outfolder )

    inp_dic = modUtils.parse_tabbed_file(
        tools.absfile( outfolder + '/templates/nr/chain_index.txt' ) )

    c.process_all( inp_dic )
