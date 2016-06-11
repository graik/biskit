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
##
"""
Prepare template coordinates for modelling.
"""

import Biskit.tools as T
from Biskit import PDBModel, PDBCleaner, CleanerError

from TemplateSearcher import TemplateSearcher
import modUtils

import os, string
import numpy as N

class TemplateCleaner:
    """
    Takes a list of PDB-files and chains identifiers.

    Returns:
      - cleaned PDB files with all chains
      - cleaned PDB files with needed chain for modeller.
      - CA trace PDB for T_coffee
      - sequences in fasta format for T_coffee

    Creates (folders and files):
      - mkdir cleaned   .. complete PDBs but cleaned
      - mkdir  modeller .. only the chains needed
      - mkdir  t_coffee .. CA traces
      - templates.fasta
    """ 

    F_RESULT_FOLDER = TemplateSearcher.F_RESULT_FOLDER

    F_CLEANED = F_RESULT_FOLDER + '/nr_cleaned/'
    F_MODELLER= F_RESULT_FOLDER + '/modeller/'
    F_COFFEE  = F_RESULT_FOLDER + '/t_coffee/'
    F_FASTA   = F_RESULT_FOLDER + '/templates.fasta'

    def __init__( self, outFolder, log=None ):
        """
        @param outFolder: output folder
        @type  outFolder: str
        @param log: None reports to STDOUT (drfault: None)
        @type  log: LogFile instance or None
        """
        self.outFolder = T.absfile( outFolder )
        self.log = log

        self.prepareFolders()


    def prepareFolders( self ):
        """
        Create folders needed by this class.
        """
        if not os.path.exists( self.outFolder + self.F_CLEANED ):
            os.mkdir( self.outFolder + self.F_CLEANED )
        if not os.path.exists( self.outFolder + self.F_MODELLER ):
            os.mkdir( self.outFolder + self.F_MODELLER )
        if not os.path.exists( self.outFolder + self.F_COFFEE ):
            os.mkdir( self.outFolder + self.F_COFFEE )


    def logWrite( self, msg, force=1 ):
        """
        Write message to log.

        @param msg: message
        @type  msg: str
        @param force: if no log, print message (default: 1)
        @type  force: 1|0 
        """
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg


    def __chain_by_id( self, model, chainId ):
        """
        Get a PDBModel with only requested chains.

        @param model: original PDBModel
        @type  model: PDBModel
        @param chainId: chain identifier
        @type  chainId: str
        
        @return: PDBModel with only the specified chain
        @rtype: PDBModel
        """
        chain_mask = model.mask( lambda a, id=chainId: a['chain_id'] == id )

        if N.sum( chain_mask ) == 0:
            raise CleanerError( "Cannot find chain %s in PDB '%s'." % \
                                (repr(chainId), model.pdbCode) )
        
        return model.compress( chain_mask )


    def fasta_sequence( self, header, s ):
        """
        Convert sequence to fasta format.
        
        @param header: fasta header
        @type  header: str
        @param s: sequence
        @type  s: str

        @return: fasta formated sequence
        @rtype: str   
        """
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
        
        @param model: PDBModel
        @type  model: PDBModel
        @param fname: filename of new PDB file
        @type  fname: str
        """
        ## write CA trace
        m_ca = model.compress( model.maskCA() )

        ## prepare sequence block

        n_lines = len( m_ca ) / 13
        if not len( m_ca ) % 13 == 0:
            n_lines += 1

        chain_id = m_ca['chain_id'][0]

        ## structure alignment program needs chain identifier
        if not chain_id:
            m_ca['chain_id'] = ['A'] * m_ca.lenAtoms()

            chain_id = 'A'  ## for sequence record

        ## create SEQRES record
        n_res = m_ca.lenResidues()

        head = []

        res_from = 0
        res_to = 13
        for i in range( n_lines ):
            s = "%4i %s %4i " % (i+1, chain_id, n_res)

            if res_to > n_res:
                res_to = n_res

            for resname in m_ca.atoms['residue_name'][ res_from: res_to ]:
                s+= " %3s" % resname

            head += [ ( 'SEQRES', s ) ]

            res_from += 13
            res_to += 13

        m_ca.writePdb( fname, headlines=head )



    def write_modeller_pdb( self, model, fname ):
        """
        Write a PDB file for modeller.

        @param model: PDBModel
        @type  model: PDBModel
        @param fname: filename of new PDB file
        @type  fname: str        
        """
        model = model.clone()
        model['chain_id'] = [''] * len(model)

        model.writePdb( fname )


    def process_all( self, file_dic, keep_hetatoms=0 ):
        """
        Process PDB files in file_dic.
        The PDB is read from:
          - L{TemplateSearcher.F_NR}

        and cleaned files are written to:
          - L{F_CLEANED}
          - L{F_COFFEE}
          - L{F_MODELLER}

        If the file L{F_CLEANED} already exists, this file is
        used to write modeller and t-coffee pdbs.

        @param file_dic: dictionary mapping filenames of pdb files to
                         the chains of interest, e.g. { fname : chain_id, }
        @type  file_dic: {str:str}
        @param keep_hetatoms: keep hetatoms (default: 0)
        @type  keep_hetatoms: 0|1
        """
        fasta = ""
        c = None

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
                self.logWrite( T.lastError() )
                self.err_cleaner = c
                self.err_file = f

        fasta_out = open( self.outFolder + self.F_FASTA, 'w' )
        fasta_out.write( fasta )
        fasta_out.close()



#############
##  TESTING        
#############
import Biskit.test as BT    

class Test(BT.BiskitTest):
    """
    Test class
    """

    def prepare(self):
        import tempfile
        import shutil
        from Biskit.LogFile import LogFile

        ## temp output directory
        self.outfolder = tempfile.mkdtemp( '_test_TemplateCleaner' )
        os.mkdir( self.outfolder +'/templates' )

        ## log file
        self.f_out = self.outfolder + '/TemplateCleaner.log'
        self.l = None
        if not self.local:
            self.l = LogFile( self.f_out, mode='w')
    
        shutil.copytree( T.testRoot() + '/Mod/project/templates/nr',
                         self.outfolder + '/templates/nr' )

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )

    def test_TemplateCleaner(self):
        """Mod.TemplateCleaner test"""
        import glob

        self.c = TemplateCleaner( self.outfolder, log=self.l)

        inp_dic = modUtils.parse_tabbed_file(
            T.absfile( self.outfolder + '/templates/nr/chain_index.txt' ) )

        self.c.process_all( inp_dic )
        
        if self.local:
            print 'TemplateCleaner log file written to: %s'% self.f_out
            globals().update( locals() )

        ## check that cleaned pdbs have been written
        f_cleaned = glob.glob(self.outfolder+ '/templates/nr_cleaned/*pdb')
        self.assert_( len(f_cleaned) > 0 )

if __name__ == '__main__':

    BT.localTest(debug=0)



