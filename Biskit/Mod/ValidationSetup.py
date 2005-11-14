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
## Contributions: olivier PERIN
## last $Author$
## last $Date$
## $Revision$

import Biskit.tools as tools
from Biskit.PDBModel import PDBModel
import re
import glob
import copy

import modUtils as MU

from TemplateSearcher import TemplateSearcher
from SequenceSearcher import SequenceSearcher
from TemplateCleaner import TemplateCleaner


import os, string

class ValidationSetup:
    """
    Takes a TemplateSearcher result folder and creates sub-projects
    with each template cluster center as target sequence to be modeled.
    In each sub-project folder a folder structures analogue to the main
    project is set up.
    The real structure is linked into the sub-project folder as reference.pdb
    """ 

    F_RESULT_FOLDER = '/validation'
    F_NR_FOLDER = SequenceSearcher.F_RESULT_FOLDER

    F_ALPHA_FOLDER = TemplateCleaner.F_COFFEE
    F_PDB_FOLDER = TemplateCleaner.F_MODELLER
    F_PDB_LINK = F_PDB_FOLDER
    F_ALPHA_LINK = F_ALPHA_FOLDER
    
    F_TEMPLATE_SEQUENCE = '/target.fasta'
    F_TCOFFEE = '/t_coffee_template_files'
    F_TEMPLATES_FASTA = '/templates.fasta'
    F_KNOWN_STRUCTURE = '/reference.pdb'


    def __init__( self, outFolder, log=None ):
        """
        log      - LogFile instance or None, None reports to STDOUT
        """
        self.outFolder = tools.absfile( outFolder )
        self.log = log

        self.prepareFolders()


    def prepareFolders( self ):
        """
        Check that all needed folders exist, if not create them
        """
        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg

        
    def cluster_result(self, chain_index = None):
        """
        Take clustering result from the file 'chain_index.txt'
        -> [str], pdb codes of templates
        """
        chain_index = chain_index or self.outFolder + \
                      TemplateSearcher.F_NR + TemplateSearcher.F_CHAIN_INDEX

        r1 = re.compile( r'([A-Z0-9]{4}).pdb' )
        index = open( "%s"%chain_index, 'r' )

        cluster_list = []
        
        string_lines = index.readlines()
        for i in string_lines:
            if( r1.search(i) ):
                code = r1.findall(i)[0]
                cluster_list.append(code)
                
        index.close()
        
        return cluster_list
    

    def createTemplatesFolder(self, validation_folder, cluster):
        os.mkdir( '%s/%s'%(validation_folder,cluster) )


    def prepare_alpha(self, cluster_list, alpha_folder = None,
                      output_folder = None):
        """
        Create a dictionary where the keys are template pdb codes and
        the value are the corresponding file names of the carbon alpha
        pdb files for ALIGNER (.alpha).
        
        cluster_list -[str], pdb codes of templates
        -> {str:str}
        """
        alpha_folder = alpha_folder or self.outFolder + self.F_ALPHA_FOLDER
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        alpha_path = glob.glob('%s/*.alpha'%alpha_folder)
        
        alpha_files = []
        for i in alpha_path:
            alpha_files.append(os.path.split(i)[1])

        alpha_dictionary = {}
        for cluster in cluster_list:
            alpha_tmp = copy.copy(alpha_path)

            for i in range( len(alpha_files) ):
                ## remove current cluster center
                if alpha_files[i][0:4] == cluster:
                    alpha_tmp.remove(alpha_tmp[i])

            alpha_dictionary.update({'%s'%cluster : alpha_tmp})
            output = open("%s/%s"%(output_folder,cluster +self.F_TCOFFEE),'w')

            for line in alpha_tmp:
                output.write(line + "\n")

            output.close()

        return alpha_dictionary
    

    def prepare_pdb(self, cluster_list, pdb_folder = None,
                    output_folder = None):
        """
        Create a dictionary which keys are templates pdb code and the value
        the different file names of pdb files for MODELLER
        cluster_list -[str], pdb codes of templates
        -> {str:str}
        """
        pdb_folder = pdb_folder or self.outFolder + self.F_PDB_FOLDER
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER
        self.pdb_path = glob.glob('%s/*.pdb'%pdb_folder)
        
        pdb_files = []
        for i in self.pdb_path:
            pdb_files.append(os.path.split(i)[1])

        pdb_dictionary = {}
        for cluster in cluster_list:
            pdb_tmp = copy.copy(self.pdb_path)

            for i in range(len(pdb_files)):
                if(pdb_files[i][0:4] == cluster):
                    pdb_tmp.remove(pdb_tmp[i])

            pdb_dictionary.update({'%s'%cluster : pdb_tmp})

        return pdb_dictionary
            

    def prepare_templatesfasta(self, cluster_list, pdb_dictionary,
                               output_folder = None):
        """
        Create 'templates.fasta' file for each template to validate
        """
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER
        
        for cluster in cluster_list:
            os.system('mkdir %s/%s'%(output_folder, cluster + \
                                     TemplateSearcher.F_RESULT_FOLDER))
            pdb_path = pdb_dictionary["%s"%cluster]
            PDBModels_list = []
            pdb_name = []
                        
            for pdb in pdb_path:
                PDBModels_list.append(PDBModel('%s'%pdb))
                pdb_name.append(os.path.split(pdb)[1][:-4])

            input_file = self.outFolder + self.F_RESULT_FOLDER + \
                         '/%s'%cluster + TemplateSearcher.F_RESULT_FOLDER \
                         + self.F_TEMPLATES_FASTA
            
            templatesfasta = open("%s"%input_file,'w')
            
            for i in range(len(PDBModels_list)):
                templatesfasta.write(">%s\n"%pdb_name[i])
                sequence = PDBModels_list[i].sequence()
                sequence = MU.format_fasta(seq = sequence)
                templatesfasta.write("%s\n"%sequence)

            templatesfasta.close()
            

    def link_pdb(self, cluster_list, pdb_dictionary, alpha_dictionary,
                 output_folder = None):
        """
        Create link in each template folder to the pdb files for MODELLER
        and for the alpha files for T-Coffee.
        """
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER
        
        for cluster in cluster_list:
            ## Modellar pdb links
            os.mkdir('%s/%s'%(output_folder, cluster + self.F_PDB_LINK))
            pdb_path = pdb_dictionary[cluster]
            
            for pdb in pdb_path:
                os.link('%s'%pdb, '%s/%s/%s'%(output_folder,
                                              cluster + self.F_PDB_LINK,
                                              os.path.split(pdb)[1]))

            ## T-Coffee alpha links
            os.mkdir('%s/%s'%(output_folder, cluster + self.F_ALPHA_FOLDER))
            alpha_path =  alpha_dictionary[cluster]
   
            for alpha in alpha_path:
                os.link('%s'%alpha, '%s/%s/%s'%(output_folder,
                                              cluster + self.F_ALPHA_LINK,
                                              os.path.split(alpha)[1]))


    def prepare_target(self,cluster, output_folder = None):
        """
        Create the 'target.fasta' file for each template to validate
        """
        output_folder = output_folder or self.outFolder + \
                        self.F_RESULT_FOLDER + '/%s/'%cluster
        target = open("%s"%(output_folder + self.F_TEMPLATE_SEQUENCE),'w')
        target.write(">target\n")

        for pdb in self.pdb_path:
            if(cluster == os.path.split(pdb)[1][0:4]):

                model = PDBModel('%s'%pdb)
                sequence = model.sequence()
                sequence = MU.format_fasta(seq = sequence)
                target.write("%s"%sequence)
                    
        target.close()


    def prepare_sequences(self, cluster, sequences_folder = None,
                          output_folder = None):
        """
        Link the 'sequences' directory from the project directory
        in each template folder
        """
        sequences_folder = sequences_folder or self.outFolder + \
                           SequenceSearcher.F_RESULT_FOLDER

        output_folder = output_folder or self.outFolder + \
                        self.F_RESULT_FOLDER + '/%s'%cluster
                  
        os.system('ln -s %s %s'%(sequences_folder , output_folder))
        

    def link_reference_pdb(self, cluster, input_folder = None,
                           output_file = None):
        """
        Create a link in each template folder with their own known
        structure 'reference.pdb'
        """
        input_folder = input_folder or self.outFolder + self.F_PDB_FOLDER
        output_file = output_file or self.outFolder + self.F_RESULT_FOLDER +\
                      '/%s/'%cluster + self.F_KNOWN_STRUCTURE
        
        files = os.listdir('%s'%input_folder)
        for pdb in files:
            if(cluster == pdb[0:4]):
                os.system('ln %s %s'%(input_folder + pdb, output_file))


    def go(self, validation_folder = None):
        """
        """
        validation_folder = validation_folder + self.F_RESULT_FOLDER or \
                            self.outFolder + self.F_RESULT_FOLDER

        cluster_list = self.cluster_result()
        for cluster in cluster_list:
            self.createTemplatesFolder(validation_folder, cluster)
            self.prepare_sequences(cluster)

        alpha_dictionary = self.prepare_alpha(cluster_list)
        pdb_dictionary = self.prepare_pdb(cluster_list)
        self.prepare_templatesfasta(cluster_list, pdb_dictionary)
        self.link_pdb(cluster_list, pdb_dictionary, alpha_dictionary)
        
        for cluster in cluster_list:
            self.prepare_target(cluster)
            self.link_reference_pdb(cluster)
      

#######
## TEST
#######
if __name__ == '__main__':

    base_folder = tools.testRoot() + '/Mod/project'

    v = ValidationSetup(outFolder = base_folder)    

    v.go(validation_folder = base_folder)  
    
