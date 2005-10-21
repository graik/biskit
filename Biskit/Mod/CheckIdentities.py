##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004, Raik Gruenberg & Johan Leckner; All rights reserved
##
## contributing authors: Olivier Perin
##
## last $Author$
## $Date$
## $Revision$

import os
import re
import linecache
import Numeric as N


import Biskit.tools as T

from Biskit import StdLog, EHandler



class Check_Identities:
    """
    This class regroups the different methods to prevent
    the use of templates with a too high percentage of identities
    """
    # standard directory for input
    F_INPUT_FOLDER = '/t_coffee'
    
    # Standard Input Files
    # t_coffee final alignment
    F_INPUT_ALNS= F_INPUT_FOLDER +'/final.pir_aln'
    

    # Standard Ouput File
    F_OUTPUT_IDENTITIES = '/identities.out'
    F_OUTPUT_IDENTITIES_INF = '/identities_info.out' 
    F_OUTPUT_IDENTITIES_COV = '/identities_cov.out'
    
    def __init__(self, outFolder='.', verbose=1):
        """
        verbose - 1||0, write intermediary files [0]
        """
        self.outFolder = T.absfile( outFolder )
        
        self.verbose = verbose

        self.sequences_name=["target"]

        
    def get_lines(self, aln_file = None):
        """
        Get all the lines from the aln file
        """
        
        aln_file = aln_file or self.outFolder + self.F_INPUT_ALNS
        
        file = open('%s'%aln_file, 'a+')

        string_lines = file.readlines() 

        return string_lines

    
    def search_length(self, string_lines, aln_file = None):
        """
        This function is able to find the length an alignment takes
        in the file 'final.pir_aln'
        """
       
        end_of_sequence=1
        
        for x in string_lines:

            if(x == '*\n'):
                endof_1st_sequence = end_of_sequence
                break

            end_of_sequence+=1

        return end_of_sequence - 3    
        

    def get_aln_sequences(self, string_lines, aln_length, aln_file = None):
        """
        Build a dictionnary with templates name and their aln
        sequences
        """

        aln_dictionnary = {}
        target_sequence = ""

        y=0
        
        for x in string_lines:

            if(x == '>P1;target\n'):
                break
            y+=1

        for s in string_lines[y+2:y+2+aln_length]:
            target_sequence += s[:-1]

        aln_dictionnary["target"] = {'nom':'target','seq':target_sequence}

        return aln_dictionnary

    
        

    def get_aln_templates(self, string_lines, aln_dict, aln_length):
        """
        This function is able to identify the different templates
        which will be used to the next step 'model.py' and to
        prepare a matrix_identities between all these sequences

        grace to the updating of the dictionnary 'self.sequences'
        is updated with those sequences
        """
       
        r1 = re.compile(r'>P1;\d')

        z = 0

        for string in string_lines[:]:

            if(r1.match(string)):

                pdb_name =""
                
                for w in string[4:]:

                    for letters in w:

                        if(letters!='\n'):
                            pdb_name += w

                        else: break

                template_sequence = ""

                for lines in string_lines[z+2:z+2+aln_length]:
                    template_sequence += lines[:-1]

                
                ## self.sequences["%s"%pdb_name] = template_sequence

                self.sequences_name.append("%s"%pdb_name)
                aln_dict["%s"%pdb_name] = {'nom':'%s'%pdb_name,'seq':template_sequence}
        

            z+=1
            
        return aln_dict


    def identities(self, aln_dictionnary):
        """
        Create a list that contains the global and template information
        identities for all the sequences in the aln_dictionnary
        """  
        for i in self.sequences_name:

            templates_name = []
            
            for name in self.sequences_name:
                if(name is not i):
                    templates_name.append(name)

            info_ID = {}
            ID = {}
            cov_ID = {}
            

            for y in self.sequences_name:

                identity = 0
                info_identity = 0
                cov_identity = 0
                nb_of_identities = 0
                nb_of_template = 0 
                template_info = []
                nb_of_residues = 0

                 
                for w in range(len(aln_dictionnary["target"]["seq"])):

                    nb_of_info_res=0

                    if(aln_dictionnary["%s"%i]["seq"][w] is not '-'):

                        nb_of_residues += 1

                        if(aln_dictionnary["%s"%i]["seq"][w] == aln_dictionnary["%s"%y]["seq"][w]):
                            nb_of_identities += 1

                        if(aln_dictionnary["%s"%y]["seq"][w] is not '-'):
                            nb_of_template += 1

    
                        for z in templates_name:

                            if(aln_dictionnary["%s"%z]["seq"][w] is not '-'):
                                nb_of_info_res += 1

                        template_info.append(nb_of_info_res)

                nb_cov_res = 0
                
                for z in range(len(template_info)):
                    if(template_info[z] != 0):
                        nb_cov_res += 1

                identity = float(float(nb_of_identities) /float( nb_of_residues)) * 100.
                info_identity = float(float(nb_of_identities) / float(nb_of_template)) * 100.
                cov_identity = float(float(nb_of_identities) / float(nb_cov_res)) * 100.
                
                
                info_ID["%s"%y] = info_identity
                ID["%s"%y] = identity
                cov_ID["%s"%y] = cov_identity
                

            aln_dictionnary["%s"%i]["info_ID"] = info_ID 
            aln_dictionnary["%s"%i]["ID"] = ID
            aln_dictionnary["%s"%i]["cov_ID"] = cov_ID

            aln_dictionnary["%s"%i]["template_info"] = template_info

        
        return aln_dictionnary        
                        
                        
    def output_identities(self, aln_dictionnary, identities_file = None, identities_info_file = None, identities_cov_file = None):

        identities_file = identities_file or self.outFolder + self.F_OUTPUT_IDENTITIES
        identities_info_file =  identities_info_file or self.outFolder + self.F_OUTPUT_IDENTITIES_INF
        identities_cov_file = identities_cov_file or self.outFolder + self.F_OUTPUT_IDENTITIES_COV

        file1 = open(identities_file, 'a+')
        file2 = open(identities_info_file, 'a+')
        file3 = open(identities_cov_file, 'a+')
    
        for i in self.sequences_name:

            file1.write(i+'\t')
            file2.write(i+'\t')
            file3.write(i+'\t')
            
            for y in (self.sequences_name):

                file1.write("%s\t"%aln_dictionnary["%s"%i]["ID"]["%s"%y])
                file2.write("%s\t"%aln_dictionnary["%s"%i]["info_ID"]["%s"%y])
                file3.write("%s\t"%aln_dictionnary["%s"%i]["cov_ID"]["%s"%y])

            file1.write('\n')
            file2.write('\n')
            file3.write('\n')
            
        file1.close()
        file2.close()
        file3.close()
        
    def go(self, output_folder = None):

        output_folder = output_folder or self.outFolder

        string_lines = self.get_lines()
        
        aln_length = self.search_length(string_lines)

        aln_dictionnary = self.get_aln_sequences(string_lines, aln_length)

        aln_dictionnary = self.get_aln_templates(string_lines, aln_dictionnary, aln_length)

        aln_dictionnary = self.identities(aln_dictionnary)

        self.output_identities(aln_dictionnary)

        return aln_dictionnary
        
        
	    
		

############
### TEST ###
############

if __name__ == '__main__':
    
    folders = os.listdir(T.absfile('~/Homstrad_final/'))


    for folder in folders:

        folders_to_validate =  os.listdir(T.absfile('~/Homstrad_final/%s/validation/'%folder))
        
        for fold in folders_to_validate:

            base_folder = T.absfile('~/Homstrad_final/%s/validation/%s'%(folder, fold))
    
            m = Check_Identities( outFolder = base_folder )

            m.go()
