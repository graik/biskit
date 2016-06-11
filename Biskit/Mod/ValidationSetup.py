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
## Contributions: olivier PERIN

"""
Setup directory structure for the validation.  
"""

import Biskit.tools as T
from Biskit.PDBModel import PDBModel
import re
import glob
import copy

import modUtils as MU

from TemplateSearcher import TemplateSearcher
from SequenceSearcher import SequenceSearcher
from TemplateCleaner import TemplateCleaner

import os.path

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
        @param outFolder: base folder 
        @type  outFolder: str
        @param log: None reports to STDOUT
        @type  log: LogFile instance or None
        """
        self.outFolder = T.absfile( outFolder )
        self.log = log

        self.prepareFolders()


    def prepareFolders( self ):
        """
        Check that all needed folders exist, if not create them
        """
        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def logWrite( self, msg, force=1 ):
        """
        Write message to log.

        @param msg: message to print
        @type  msg: str
        """
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg


    def cluster_result(self, chain_index = None):
        """
        Take clustering result from the file 'chain_index.txt'

        @param chain_index: file with clustering results
                            (default: None-> L{TemplateSearcher.F_CHAIN_INDEX})
        @type  chain_index: 

        @return: pdb codes of templates
        @rtype: [str]
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
        """
        Create folders for the templates to be used for the validation.

        @param validation_folder: top folder for the validation
        @type  validation_folder: str
        @param cluster: name for validation subfolder
                        (e.g. pdb code of cluster center)
        @type  cluster: str
        """
        try:
            os.mkdir( '%s/%s'%(validation_folder,cluster) )
        except:
            print 'Folder %s/%s already exists.'\
                  %(self.F_RESULT_FOLDER, cluster) 


    def prepare_alpha(self, cluster_list, alpha_folder = None,
                      output_folder = None):
        """
        Create a dictionary where the keys are template pdb codes and
        the value are the corresponding file names of the carbon alpha
        pdb files for ALIGNER (.alpha).

        @param cluster_list: pdb codes of templates
        @type  cluster_list: [str]
        @param alpha_folder: folder with template CA-trace files
                             (default: None -> L{F_ALPHA_FOLDER})
        @type  alpha_folder: str
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str

        @return: dictionary mapping pdb code to CA-trace files
        @rtype: {str:str}
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

        @param cluster_list: pdb codes of templates
        @type  cluster_list: [str]
        @param pdb_folder: folder with Modeller pdb files
                             (default: None -> L{F_PDB_FOLDER})
        @type  pdb_folder: str
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str

        @return: dictionary mapping pdb code to pdb files used by  Modeller
        @rtype: {str:str}
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

        @param cluster_list: pdb codes of templates
        @type  cluster_list: [str]
        @param pdb_dictionary: dictionary mapping pdb code to pdb files
                               used by Modeller
        @type  pdb_dictionary: {str:str}
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str        
        """
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        for cluster in cluster_list:
            folder = '%s/%s'%(output_folder, cluster + \
                              TemplateSearcher.F_RESULT_FOLDER)
            if not os.path.exists( folder ):
                os.mkdir( folder)
            else:
                print 'Directory %s exists, skipping'%( cluster + \
                                                        TemplateSearcher.F_RESULT_FOLDER)

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

        @param cluster_list: pdb codes of templates
        @type  cluster_list: [str]
        @param pdb_dictionary: dictionary mapping pdb code to pdb files
                               used by Modeller
        @type  pdb_dictionary: {str:str}
        @param alpha_dictionary: dictionary mapping pdb code to CA-trace files
        @type  alpha_dictionary: {str:str}
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str            
        """
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        for cluster in cluster_list:
            ## Modellar pdb links
            folder = '%s/%s'%(output_folder, cluster + self.F_PDB_LINK)
            if not os.path.exists( folder ):
                os.mkdir( folder )
            else:
                print 'Directory %s exists, skipping'%\
                      (cluster + self.F_PDB_LINK)

            pdb_path = pdb_dictionary[cluster]

            for pdb in pdb_path:
                target = '%s/%s/%s'%(output_folder, cluster + self.F_PDB_LINK,
                                     os.path.split(pdb)[1])
                if not os.path.exists( target ):
                    os.link('%s'%pdb, target )
                else:
                    print 'File exists %s/%s no link made.'%\
                          (self.F_PDB_LINK,
                           os.path.split(pdb)[1])

            ## T-Coffee alpha links
            folder = '%s/%s'%(output_folder, cluster + self.F_ALPHA_FOLDER)
            if not os.path.exists( folder ):
                os.mkdir( folder )
            else:
                print '##'

            alpha_path =  alpha_dictionary[cluster]

            for alpha in alpha_path:
                target = '%s/%s/%s'%(output_folder,
                                     cluster + self.F_ALPHA_LINK,
                                     os.path.split(alpha)[1])
                if not os.path.exists( target ):
                    os.link('%s'%alpha, target )
                else:
                    print 'File exists %s/%s no link made.'%\
                          (self.F_ALPHA_LINK,
                           os.path.split(pdb)[1])                    


    def prepare_target(self, cluster, output_folder = None):
        """
        Create the 'target.fasta' file for each template to validate

        @param cluster: name of the cluster which is used for the
                        foldder name in which the validation is run.
        @type  cluster: str
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str   
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

        @param cluster: name of the cluster which is used for the
                        folder name in which the validation is run.
        @type  cluster: str
        @param sequences_folder: folder with sequences (default: None ->
                                 L{SequenceSearcher.F_RESULT_FOLDER})
        @type  sequences_folder: str
        @param output_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  output_folder: str         
        """
        sequences_folder = sequences_folder or self.outFolder + \
                         SequenceSearcher.F_RESULT_FOLDER

        output_folder = output_folder or self.outFolder + \
                      self.F_RESULT_FOLDER + '/%s'%cluster + \
                      SequenceSearcher.F_RESULT_FOLDER

        if not os.path.exists( output_folder ):
            ## os.link doesn't seem to work with folders
            os.system('ln -s %s %s'%(sequences_folder , output_folder))

        else:
            print 'Folder %s already exists, linking skipped.\
            '%(self.F_RESULT_FOLDER + '/%s'%cluster +
               SequenceSearcher.F_RESULT_FOLDER )


    def link_reference_pdb(self, cluster, input_folder = None,
                           output_file = None):
        """
        Create a link in each template folder with their own known
        structure 'reference.pdb'

        @param cluster: name of the cluster which is used for the
                        foldder name in which the validation is run.
        @type  cluster: str
        @param input_folder: folder with pdb files
                             (default: None -> L{F_PDB_FOLDER})
        @type  input_folder: str           
        @param output_file: target file
        @type  output_file: str           
        """
        input_folder = input_folder or self.outFolder + self.F_PDB_FOLDER
        output_file = output_file or self.outFolder + self.F_RESULT_FOLDER +\
                    '/%s/'%cluster + self.F_KNOWN_STRUCTURE

        files = os.listdir('%s'%input_folder)
        for pdb in files:
            if(cluster == pdb[0:4]):
                if not os.path.exists( output_file ):
                    os.link( input_folder + pdb, output_file)
#                os.system('ln %s %s'%(input_folder + pdb, output_file))


    def go(self, validation_folder = None):
        """
        Greate validation directory setup.

        @param validation_folder: top output folder
                             (default: None -> L{F_RESULT_FOLDER})
        @type  validation_folder: str         
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

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_ValidationSetup' )
        os.mkdir( self.outfolder +'/templates' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/nr',
                         self.outfolder + '/templates/nr' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/modeller',
                         self.outfolder + '/templates/modeller' )    



    def test_ValidationSetup(self):
        """Mod.ValidationSetup test"""
        v = ValidationSetup( outFolder = self.outfolder )    

        v.go( validation_folder =self.outfolder )  

        if self.local and self.DEBUG:
            print 'The validation project can be found in %s/validation'%\
                  self.outfolder


    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1)

if __name__ == '__main__':

    BT.localTest()
