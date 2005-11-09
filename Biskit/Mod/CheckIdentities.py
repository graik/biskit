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

## Contributions: Olivier PERIN
## last $Author$
## last $Date$
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
    
    # Standard Ouput Files
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
        Retrieve the lines from an aln file
        """
        aln_file = aln_file or self.outFolder + self.F_INPUT_ALNS
        file = open('%s'%aln_file, 'r')
        string_lines = file.readlines()
        file.close()

        return string_lines

    
    def search_length(self, string_lines, aln_file = None):
        """
        This function returns the length an alignment in the
        file 'final.pir_aln'
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
        Create a dictionary with the name of the target (i.e 'target')
        and sequence from the final output from T-Coffee (final.pir_aln).
        -> {'name':'target, 'seq': 'sequence of the target'}
        """
        aln_dictionary = {}
        target_sequence = ""

        y=0
        for x in string_lines:           
            if(x == '>P1;target\n'):
                break
            y+=1

        for s in string_lines[y+2:y+2+aln_length]:
            target_sequence += s[:-1]

        aln_dictionary["target"] = {'name':'target','seq':target_sequence}

        return aln_dictionary


    def get_aln_templates(self, string_lines, aln_dict, aln_length):
        """
        Add information about the name of the template sequences and
        the sequence that is aligned to the template. Data taken from
        the T-Coffee alignment (final.pir_aln).
        -> { str :{'name':str, 'seq': str} }
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
                    
                self.sequences_name.append("%s"%pdb_name)
                aln_dict["%s"%pdb_name] = {'name':'%s'%pdb_name,
                                           'seq':template_sequence}
        
            z+=1
            
        return aln_dict


    def identities(self, aln_dictionary):
        """
        Create a dictionary that contains information about all the
        alignments in the aln_dictionar using pairwise comparisons.

        -> a dictionary of dictionaries with the sequence name as the
           top key. Each sub dictionary then has the keys: 
             'name' - str, sequence name
             'seq' - str, sequence of
             'template_info' - list of the same length as the 'key'
                  sequence excluding deletions. The number of sequences
                  in tha multiple alignmentthat contain information at
                  this position.
             'ID' - dict, sequence identity in precent comparing the
                  'key'  sequence all other sequences (excluding deletions)
             'info_ID' - dict, same as 'ID' but compared to the template
                  sequence length (i.e excluding deletions and insertions
                  in the 'key' sequence )
             'cov_ID' - dict, same as 'info_ID' but insertions are defined
                  comparing to all template sequences (i.e where
                  'template_info' is zero )
        """
        ## loop over all sequences in alignment
        for i in self.sequences_name:
            template_names = []

            ## don't compare to self, remove current sequence
            for name in self.sequences_name:
                if(name is not i):
                    template_names.append(name)

            ## loop over all sequences in alignment
            info_ID, ID, cov_ID  = {}, {}, {}
            for y in self.sequences_name:
                identity = 0
                info_identity = 0
                cov_identity = 0
                nb_of_identities = 0
                nb_of_template = 0 
                template_info = []
                nb_of_residues = 0

                ## loop over the full length of the alignment
                for w in range(len(aln_dictionary["target"]["seq"])):

                    ## skipp deletions
                    nb_of_info_res=0
                    if(aln_dictionary[i]["seq"][w] is not '-'):
                        nb_of_residues += 1
                        
                        ## count identities
                        if(aln_dictionary[i]["seq"][w] == \
                           aln_dictionary[y]["seq"][w]):
                            nb_of_identities += 1

                        ## length excluding insertions
                        if(aln_dictionary[y]["seq"][w] is not '-'):
                            nb_of_template += 1

                        ## loop over all sequences but self
                        for z in template_names:
                            ## count how many sequences contain alignment
                            ## information at this position
                            if(aln_dictionary[z]["seq"][w] is not '-'):
                                nb_of_info_res += 1

                        template_info.append(nb_of_info_res)

                ## number of positions in which any other sequence
                ## contains alignment information
                nb_cov_res = N.sum( N.greater(template_info, 0) )

                ## calculate identities
                info_ID[y] = 100. * nb_of_identities / nb_of_template
                ID[y]      = 100. * nb_of_identities / nb_of_residues
                cov_ID[y]  = 100. * nb_of_identities / nb_cov_res
                
            aln_dictionary[i]["info_ID"] = info_ID 
            aln_dictionary[i]["ID"] = ID
            aln_dictionary[i]["cov_ID"] = cov_ID
            aln_dictionary[i]["template_info"] = template_info

        return aln_dictionary        


    def __writeId( self, name, dic, key, description ):
        """
        Write an identity matrix to file.
        """
        f = open( name, 'w' )
        f.write( description +'\n\n')
        
        ## write header
        f.write('%s'%(' '*5))
        for s in self.sequences_name:
            f.write('%8s'%s)
        f.write('\n')

        ## write matrix
        for s in self.sequences_name:
            f.write('%5s'%s)
            for t in self.sequences_name:
                f.write("%8.2f"%dic[s][key][t])
            f.write('\n')
                    
            
    def output_identities(self, aln_dictionary, identities_file = None,
                          identities_info_file = None,
                          identities_cov_file = None):
        """
        Writes three files to disk with identity info about the current
        multiple alignment.
        """
        ## filenames to create
        identities_file = identities_file or \
                          self.outFolder + self.F_OUTPUT_IDENTITIES
        identities_info_file =  identities_info_file or \
                               self.outFolder + self.F_OUTPUT_IDENTITIES_INF
        identities_cov_file = identities_cov_file or \
                              self.outFolder + self.F_OUTPUT_IDENTITIES_COV
        
        head_ID = "Sequence identity in precent comparing a sequence to another \n(excluding deletions in the first sequence)"

        head_info_ID ="Sequence identity in precent comparing a sequence to another \n(excluding deletions and insertions in the first sequence)"

        head_conv_ID ="Sequence identity in precent comparing a sequence to another \n(excluding deletions and insertions in the first sequence but only \nwhen the first sequence doesn't match any other sequence in the \nmultiple alignment )"
                  
        self.__writeId( identities_file, aln_dictionary,
                        'ID' , head_ID )
        self.__writeId( identities_info_file, aln_dictionary,
                        'info_ID' , head_info_ID )
        self.__writeId( identities_cov_file, aln_dictionary,
                        'cov_ID', head_conv_ID )


    def go(self, output_folder = None):

        output_folder = output_folder or self.outFolder

        string_lines = self.get_lines()
        
        aln_length = self.search_length( string_lines )

        ## get information about the target sequence from the alignment
        aln_dictionary = self.get_aln_sequences( string_lines,
                                                 aln_length )

        ## add information about the aligned templates from the alignment
        aln_dictionary = self.get_aln_templates( string_lines,
                                                 aln_dictionary,
                                                 aln_length )

        aln_dictionary = self.identities( aln_dictionary )

        self.output_identities( aln_dictionary )

        return aln_dictionary
        
        
	    
		

############
### TEST ###
############

if __name__ == '__main__':
    
    m = Check_Identities( T.testRoot() + '/Mod/project')
    m.go()

    val_root =  T.testRoot() + '/Mod/project/validation'
    folders = os.listdir( val_root )

    for f in folders:
        m = Check_Identities( outFolder = val_root + '/' + f )
        m.go()
