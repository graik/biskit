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

## Contributions: Olivier PERIN (first version)
## last $Author$
## last $Date$
## $Revision$

"""
Check sequence identity between templates.
"""

import os
import re
import numpy as N

import Biskit.tools as T

from Biskit import StdLog, EHandler


class CheckIdentities:
    """
    This class provides different methods to prevent
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


    def __init__(self, outFolder='.', alignment=None, verbose=1):
        """
        @param outFolder: base folder
        @type  outFolder: str
        @param alignment: path to alignment file in non-default location
        @type  alignment: str
        @param verbose: write intermediary files (default: 0)
        @type  verbose: 1|0
        """
        self.outFolder = T.absfile( outFolder )
        self.verbose = verbose
        self.f_aln = alignment or self.outFolder + self.F_INPUT_ALNS
        self.sequences_name=["target"]

        #: will hold result dict after go()
        self.result = None


    def get_lines(self, aln_file = None):
        """
        Retrieve the lines from an aln file

        @param aln_file: aln file (default: None -> L{F_INPUT_ALNS})
        @type  aln_file: str

        @return: file as list of strings
        @rtype: [str]
        """
        return open( self.f_aln, 'r').readlines()


    def search_length(self, string_lines):
        """
        Return the number of sequences in the alignment 'final.pir_aln'.

        @param string_lines: aln file as list of string
        @type  string_lines: [str]

        @return: length of alignment
        @rtype: int
        """
        end_of_sequence=1

        for x in string_lines:
            if(x == '*\n'):
                break

            end_of_sequence+=1

        return end_of_sequence - 3    


    def get_aln_sequences(self, string_lines, aln_length):
        """
        Create a dictionary with the name of the target (i.e 'target')
        and sequence from the final output from T-Coffee (final.pir_aln).

        @param string_lines: aln file as list of string
        @type  string_lines: [str]
        @param aln_length: length of alignment
        @type  aln_length: int

        @return: alignment dictionary
                 e.g. {'name':'target, 'seq': 'sequence of the target'}
        @rtype: dict
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

        @param string_lines: aln file as list of string
        @type  string_lines: [str]
        @param aln_dict: alignment dictionary
        @type  aln_dict: dict

        @return: template alignment dictionary
                 e.g. { str :{'name':str, 'seq': str} }
        @rtype: dict
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
        alignments in the aln_dictionary using pairwise comparisons.

        @param aln_dictionary: alignment dictionary
        @type  aln_dictionary: dict

        @return: a dictionary of dictionaries with the sequence name as the
        top key. Each sub dictionary then has the keys: 
         - 'name' - str, sequence name
         - 'seq' - str, sequence of
         - 'template_info' - list of the same length as the 'key'
             sequence excluding deletions. The number of sequences
             in the multiple alignment that contain information at
             this position.
         - 'ID' - dict, sequence identity in percent comparing the
            'key'  sequence to all other sequences (excluding deletions)
         - 'info_ID' - dict, same as 'ID' but compared to the template
             sequence length (i.e excluding deletions and insertions
             in the 'key' sequence )
         - 'cov_ID' - dict, same as 'info_ID' but insertions are defined
             comparing to all template sequences (i.e where
             'template_info' is zero )
        @rtype: dict
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
##                identity = 0
##                info_identity = 0
##                cov_identity = 0
                nb_of_identities = 0
                nb_of_template = 0 
                template_info = []
                nb_of_residues = 0

                ## loop over the full length of the alignment
                for w in range(len(aln_dictionary["target"]["seq"])):

                    ## skip deletions
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
                info_ID[y] = ID[y] = cov_ID[y] = 0
                ## RAIK: Hack, nb_of_... can turn 0 for fragmented alignments
                if nb_of_template:
                    info_ID[y] = 100. * nb_of_identities / nb_of_template
                if nb_of_residues:
                    ID[y]      = 100. * nb_of_identities / nb_of_residues
                if nb_cov_res:
                    cov_ID[y]  = 100. * nb_of_identities / nb_cov_res

            aln_dictionary[i]["info_ID"] = info_ID 
            aln_dictionary[i]["ID"] = ID
            aln_dictionary[i]["cov_ID"] = cov_ID
            aln_dictionary[i]["template_info"] = template_info

        return aln_dictionary        


    def __writeId( self, name, dic, key, description ):
        """
        Write an sequence identity matrix to file.

        @param name: file name 
        @type  name: str
        @param dic: alignment dictionary
        @type  dic: dict
        @param key: key in dictionary to write
        @type  key: key
        @param description: description to go into file (first line)
        @type  description: str
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


    def write_identities(self, identities_file = None,
                         identities_info_file = None,
                         identities_cov_file = None):
        """
        Writes three files to disk with identity info about the current
        multiple alignment.

        @param aln_dictionary: alignment dictionary
        @type  aln_dictionary: dict
        @param identities_file: name for file with sequence identity
                                in percent comparing a sequence to another
                                (excluding deletions in the first sequence)
                                (default: None -> L{F_OUTPUT_IDENTITIES})
        @type  identities_file: str
        @param identities_info_file: name for file with sequence identity in
                                     percent comparing a sequence to another
                                     (excluding deletions and insertions in
                                     the first sequence)
                                     (default: None -> L{F_OUTPUT_IDENTITIES_INF})
        @type  identities_info_file: str
        @param identities_cov_file: name for file with sequence identity in
                                    percent comparing a sequence to another
                                    (excluding deletions and insertions in
                                    the first sequence but only when the
                                    first sequence doesn't match any other
                                    sequence in the multiple alignment)
                                    (default: None -> L{F_OUTPUT_IDENTITIES_COV})
        @type  identities_cov_file: str
        """
        ## filenames to create
        identities_file = identities_file or \
                        self.outFolder + self.F_OUTPUT_IDENTITIES
        identities_info_file =  identities_info_file or \
                             self.outFolder + self.F_OUTPUT_IDENTITIES_INF
        identities_cov_file = identities_cov_file or \
                            self.outFolder + self.F_OUTPUT_IDENTITIES_COV

        head_ID = """Pairwise sequence identity in percent (excluding deletions
in the reference) -- for that reason, the matrix is *not* symetric and
should be read row-wise (the first line is the most interesting one).
"""

        head_info_ID ="Sequence identity in percent comparing a sequence to another \n(excluding deletions and insertions in the first sequence)"

        head_conv_ID ="Sequence identity in percent comparing a sequence to another \n(excluding deletions and insertions in the first sequence but only \nwhen the first sequence doesn't match any other sequence in the \nmultiple alignment )"

        assert self.result is not None, 'no alignment result to write, call go() first'

        self.__writeId( identities_file, self.result,
                        'ID' , head_ID )
        self.__writeId( identities_info_file, self.result,
                        'info_ID' , head_info_ID )
        self.__writeId( identities_cov_file, self.result,
                        'cov_ID', head_conv_ID )


    def go(self, output_folder = None):
        """
        Perform sequence comparison.

        @param output_folder: output folder
        @type  output_folder: str
        """
        output_folder = output_folder or self.outFolder

        string_lines = self.get_lines()

        aln_length = self.search_length( string_lines )

        ## get information about the target sequence from the alignment
        self.result = self.get_aln_sequences( string_lines,
                                              aln_length )

        ## add information about the aligned templates from the alignment
        self.result = self.get_aln_templates( string_lines,
                                              self.result,
                                              aln_length )

        self.result = self.identities( self.result )

##         self.write_identities()

        return self.result



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
        self.outfolder = tempfile.mkdtemp( '_test_CheckIdentities' )
        os.mkdir( self.outfolder +'/t_coffee' )

        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     self.outfolder + '/t_coffee' )    

    def test_CheckIdentities(self):
        """Mod.CheckIdentities test"""
        self.m = CheckIdentities( self.outfolder )
        self.result = self.m.go()
        self.m.write_identities()

        if self.local and self.DEBUG:
            self.log.add("""The result from the template comparison can be found in the three files %s, %s and %s that reside in the folder %s"""\
                         %(self.m.F_OUTPUT_IDENTITIES[1:],
                           self.m.F_OUTPUT_IDENTITIES_INF[1:],
                           self.m.F_OUTPUT_IDENTITIES_COV[1:],
                           self.outfolder ) )

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )


if __name__ == '__main__':

    BT.localTest()
