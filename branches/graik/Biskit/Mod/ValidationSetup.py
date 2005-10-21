##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: olivier PERIN
## last $Author$
## last $Date$


import Biskit.tools as tools
from Biskit.PDBModel import PDBModel
import re
import glob
import copy

from TemplateSearcher import TemplateSearcher
from SequenceSearcher import SequenceSearcher
from TemplateCleaner import TemplateCleaner


import os, string

class ValidationSetup:
    """
    Take a TemplateSearcher result folder. Create sub-projects to model each
    template cluster as target. That is, the cluster has to be removed from
    the template files and the representative structure of the cluster should
    become target.
    """ 

    F_RESULT_FOLDER = '/validation'
    F_NR_FOLDER = SequenceSearcher.F_RESULT_FOLDER

    F_ALPHA_FOLDER = TemplateCleaner.F_COFFEE
    F_PDB_FOLDER = TemplateCleaner.F_MODELLER
    F_PDB_LINK = F_PDB_FOLDER
    
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
        Take clustering result from the file
        'chain_index.txt'
        -> [str], pdb codes of templates
        """

        chain_index = chain_index or self.outFolder + TemplateSearcher.F_NR + TemplateSearcher.F_CHAIN_INDEX

        index = open("%s"%chain_index,'a+')

        string_lines = index.readlines()
        
        cluster_list = []

        r1 = re.compile(r'([A-Z0-9]{4}).pdb')

        for i in string_lines:

            if(r1.search(i)):
                code = r1.findall(i)[0]
                cluster_list.append(code)

        return cluster_list



    def createTemplatesFolder(self, validation_folder, cluster):

        os.system('mkdir %s/%s'%(validation_folder,cluster))
        


    def prepare_alpha(self, cluster_list, alpha_folder = None, output_folder = None):
        """
        Create a dictionnary which keys are templates pdb code and the value
        the different file names of pdb alpha files for ALIGNER
        cluster_list -[str], pdb codes of templates
        -> {str:str}
        """
        
        alpha_folder = alpha_folder or self.outFolder + self.F_ALPHA_FOLDER

        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        alpha_path = glob.glob('%s/*.alpha'%alpha_folder)
        
        alpha_files = []

        for i in alpha_path:
            alpha_files.append(os.path.split(i)[1])

        alpha_dictionnary = {}

        for cluster in cluster_list:

            alpha_tmp = copy.copy(alpha_path)

            for i in range(len(alpha_files)):

                if(alpha_files[i][0:4] == cluster):
                    alpha_tmp.remove(alpha_tmp[i])

            alpha_dictionnary.update({'%s'%cluster : alpha_tmp})

            output = open("%s/%s"%(output_folder,cluster + self.F_TCOFFEE),'a+')

            for line in alpha_tmp:
                output.write(line + "\n")

            output.close()

        return alpha_dictionnary
    

    def prepare_pdb(self, cluster_list, pdb_folder = None, output_folder = None):
        """
        Create a dictionnary which keys are templates pdb code and the value
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


        pdb_dictionnary = {}
        
        for cluster in cluster_list:

            pdb_tmp = copy.copy(self.pdb_path)

            for i in range(len(pdb_files)):

                if(pdb_files[i][0:4] == cluster):
                    pdb_tmp.remove(pdb_tmp[i])

            pdb_dictionnary.update({'%s'%cluster : pdb_tmp})

           ##  output = open("%s/%s"%(output_folder, cluster + self.F_MODLLER),'a+')
            

##             for line in pdb_tmp:
##                 output.write(line + "\n")

##             output.close()

        return pdb_dictionnary
            


    def prepare_templatesfasta(self, cluster_list, pdb_dictionnary, output_folder = None):
        """
        Create 'templates.fasta' file for each template to validate
        """

        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER
        
        for cluster in cluster_list:

            os.system('mkdir %s/%s'%(output_folder, cluster + TemplateSearcher.F_RESULT_FOLDER))
            
            pdb_path = pdb_dictionnary["%s"%cluster]
            PDBModels_list = []
            pdb_name = []

                        
            for pdb in pdb_path:

                PDBModels_list.append(PDBModel('%s'%pdb))

                pdb_name.append(os.path.split(pdb)[1][:-4])
        

            input_file = self.outFolder + self.F_RESULT_FOLDER + '/%s'%cluster + TemplateSearcher.F_RESULT_FOLDER + self.F_TEMPLATES_FASTA
            
            templatesfasta = open("%s"%input_file,'a+')
    
            
            for i in range(len(PDBModels_list)):

                templatesfasta.write(">%s\n"%pdb_name[i])

                sequence = PDBModels_list[i].sequence()
                sequence = self.write_fasta(seq = sequence)

                templatesfasta.write("%s\n"%sequence)

            templatesfasta.close()

            

    def link_pdb(self, cluster_list, pdb_dictionnary, output_folder = None):
        """
        Create link in each template folder to the pdb files for MODELLER
        """

        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        for cluster in cluster_list:

            os.system('mkdir %s/%s'%(output_folder, cluster + self.F_PDB_LINK))
            
            pdb_path = pdb_dictionnary[cluster]
            
            for pdb in pdb_path:

                #os.system('ln  %s %s/%s/%s'%(pdb, output_folder, cluster +self.F_PDB_LINK,os.path.split(pdb)[1]))
                os.link('%s'%pdb, '%s/%s/%s'%(output_folder, cluster +self.F_PDB_LINK,os.path.split(pdb)[1]))


    def prepare_target(self,cluster, output_folder = None):
        """
        Create the 'target.fasta' file for each template to validate
        """

        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER + '/%s/'%cluster

        target = open("%s"%(output_folder + self.F_TEMPLATE_SEQUENCE),'a+')

        target.write(">target\n")

        for pdb in self.pdb_path:
                
            if(cluster == os.path.split(pdb)[1][0:4]):

                model = PDBModel('%s'%pdb)
                sequence = model.sequence()
                sequence = self.write_fasta(seq = sequence)
                target.write("%s"%sequence)
                    
        target.close()



    def prepare_sequences(self, cluster, sequences_folder = None, output_folder = None):
        """
        Link the 'sequences' directory from the project directory in each template folder
        """

        sequences_folder = sequences_folder or self.outFolder + SequenceSearcher.F_RESULT_FOLDER

        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER + '/%s'%cluster
                  
        os.system('ln -s %s %s'%(sequences_folder , output_folder))
        



    def write_fasta(self, seq, width=60):
        """
        Transform a given sequence in fasta format
        seq -str, sequence
        -> str, string sequence in fasta format
        """

        fasta_sequence = ""
        
        for i in xrange(0,len(seq),width):
            
            fasta_sequence += seq[i:i+width]

            if(i+width>=len(seq)):
                pass

            else:
                fasta_sequence += "\n"

        return fasta_sequence



    def link_reference_pdb(self, cluster, input_folder = None, output_file = None):
        """
        Create a link in each template folder with their own known structure 'reference.pdb'
        """

        input_folder = input_folder or self.outFolder + self.F_PDB_FOLDER

        output_file = output_file or self.outFolder + self.F_RESULT_FOLDER + '/%s/'%cluster + self.F_KNOWN_STRUCTURE
        
        files = os.listdir('%s'%input_folder)

        for pdb in files:

            if(cluster == pdb[0:4]):

                os.system('ln %s %s'%(input_folder + pdb, output_file))



    def go(self, validation_folder = None):
    
        validation_folder = validation_folder + self.F_RESULT_FOLDER or self.outFolder + self.F_RESULT_FOLDER

        cluster_list = self.cluster_result()
        
        for cluster in cluster_list:

            self.createTemplatesFolder(validation_folder, cluster)

            self.prepare_sequences(cluster)


        alpha_dictionnary = self.prepare_alpha(cluster_list)

        pdb_dictionnary = self.prepare_pdb(cluster_list)

        self.prepare_templatesfasta(cluster_list, pdb_dictionnary)

        self.link_pdb(cluster_list, pdb_dictionnary)

        
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
    
