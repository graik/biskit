## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

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
##

## Contributions: olivier PERIN
##
## last $Author$
## last $Date$
## $Revision$

"""
Analyze model quality
"""

import Biskit.tools as T
from Biskit.PDBModel import PDBModel
from Biskit.Mod.Benchmark import Benchmark
from Biskit.Mod.ValidationSetup import ValidationSetup as VS
from Biskit.Mod.CheckIdentities import CheckIdentities as CI
from Biskit.Mod.Modeller import Modeller

import os
from string import *
import numpy as N


class Analyse:
    """
    Create a folder named analyse in the project root folder
    that contains the following results:

    GLOBAL: global rmsd: all atoms, c-alpha only, percentage of
    identities, Modeller score, and the number of templates

    LOCAL: results of the cross validation, rmsd per residue c-alpha
    only for each templates and the mean rmsd

    3D structure: pickle down the final.pdb that is the best model of
    the project and the mean rmsd is set inside (temperature_factor)
    """

    F_RESULT_FOLDER = '/analyse'
    F_TEMPLATE_FOLDER = VS.F_RESULT_FOLDER

    F_PDBModels = Benchmark.F_PDBModels_OUT
    F_MODELS = Modeller.F_RESULT_FOLDER + Modeller.F_PDBModels

    F_INPUT_ALNS= '/t_coffee/final.pir_aln'

    F_INPUT_RMSD = Benchmark.F_RESULT_FOLDER
    F_RMSD_AA = Benchmark.F_RMSD_AA
    F_RMSD_CA = Benchmark.F_RMSD_CA


    F_OUTPUT_VALUES = F_RESULT_FOLDER + '/global_results.out'
    F_CROSS_VAL = F_RESULT_FOLDER + '/local_results.out'
    F_FINAL_PDB = F_RESULT_FOLDER + '/final.pdb'


    def __init__( self, outFolder, log=None ):
        """
        @param outFolder: base folder for output
        @type  outFolder: str
        @param log: None reports to STDOUT
        @type  log: LogFile instance or None
        """
        self.outFolder = T.absfile( outFolder )
        self.log = log

        self.prepareFolders()


    def prepareFolders( self ):
        """
        Create folders needed by this class.
        """
        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def parseFile( self, name ):
        """
        Parse a identity matrix file

        @param name: file to parse
        @type  name: str
        
        @return: contents of parsed file
        @rtype: [[str]]
        """
        f = open( name, 'r')
        result = []
        lines = f.readlines()

        for l in lines:
            if not l[0] == '#':
                r =[]
                for s in l.split():
                    try:
                        r += [float(s)]
                    except:
                        pass
                if len(r)>=1:
                    result += [ r ]

        f.close()
        return result


    def __listDir( self, path ):
        """
        List all the files and folders in a directory
        with the exceprion of ...
        
        @param path: dir to list
        @type  path: str

        @return: list of files
        @rtype: [str]
        """
        files = os.listdir( path )
        if 'CVS' in files:
            files.remove('CVS')
            
        return files
    

######################################################################
####  GLOBAL RESTULTS: RMSD_AA, RMSD_CA,
####                   %ID(mean of the templates),
####                   Nb of Templates 

    def global_rmsd_aa(self, validation_folder= None):
        """
        Global RMSD values.
        
        @param validation_folder: folder vith validation data
                      (defult: None S{->} outFolder/L{F_TEMPLATE_FOLDER})
        @type  validation_folder: str

        @return: two dictionaries:
                  - rmsd_aa_wo_if: global all atom rmsd for each
                    template without iterative fitting
                  - rmsd_aa_if: global all atom rmsd for each
                    templates with iterative fitting
        @rtype: dict, dict      
        """
        validation_folder = validation_folder or self.outFolder + \
                            self.F_TEMPLATE_FOLDER

        folders = self.__listDir(validation_folder)

        rmsd_aa_wo_if = {}
        rmsd_aa_if = {}

        for folder in folders:

            file = "%s/%s"%(validation_folder, folder + self.F_RMSD_AA)
            lst = self.parseFile( file )

            rmsd_aa_wo_if[folder] = [ lst[0][0] ]
            rmsd_aa_if[folder]    = [ lst[0][1], lst[0][2]*100.]

        return rmsd_aa_wo_if, rmsd_aa_if


    def global_rmsd_ca(self, validation_folder= None):
        """
        Global RMSD CA values.
        
        @param validation_folder: folder vith validation data
                       (defult: None S{->} outFolder/L{F_TEMPLATE_FOLDER})
        @type  validation_folder: str

        @return: two dictionaries:
                  - rmsd_ca_wo_if: global CA rmsd for each template
                    without iterative fitting
                  - rmsd_ca_if: global CA rmsd for each template
                    with iterative fitting
        @rtype: dict, dict  
        """
        validation_folder = validation_folder or self.outFolder + \
                            self.F_TEMPLATE_FOLDER

        folders = self.__listDir(validation_folder)

        rmsd_ca_wo_if = {}
        rmsd_ca_if = {}

        for folder in folders:

            file = "%s/%s"%(validation_folder, folder + self.F_RMSD_CA)

            lst = self.parseFile( file )

            rmsd_ca_wo_if[folder] = [ lst[0][0] ]
            rmsd_ca_if[folder]    = [ lst[0][1], lst[0][2]*100.]

        return rmsd_ca_wo_if, rmsd_ca_if


    def get_identities(self, nb_templates, validation_folder = None):
        """
        Calculate the mean of the percentage of identities for each
        template with the others.
        
        @param nb_templates: number of templates used in the cross-validation
        @type  nb_templates: int
        @param validation_folder: folder vith validation data
                           (defult: None S{->} outFolder/L{F_TEMPLATE_FOLDER})
        @type  validation_folder: str
        
        @return: dictionary with mean percent identities for each template
        @rtype: {str:float}
        """

        validation_folder = validation_folder or self.outFolder + \
                            self.F_TEMPLATE_FOLDER

        folders = self.__listDir(validation_folder)
        identities = {}

        for folder in folders:
            file = "%s/%s"%(validation_folder, folder + \
                            CI.F_OUTPUT_IDENTITIES_COV)

            lst = self.parseFile( file )

            ## identity to mean template
            identities[folder] = N.sum(lst[0][1:])/nb_templates

        return identities


    def get_score(self, validation_folder = None):
        """
        Get the best global modeller score for each template re-modeled

        @param validation_folder: folder vith validation data
                             (defult: None S{->} outFolder/L{F_TEMPLATE_FOLDER})
        @type  validation_folder: str
        
        @return: dictionary with modeller score for each template
        @rtype: {str:float}
        """
        validation_folder = validation_folder or self.outFolder + \
                            self.F_TEMPLATE_FOLDER

        folders = self.__listDir(validation_folder)
        score = {}

        for folder in folders:
            file = "%s/%s"%(validation_folder, folder + Modeller.F_SCORE_OUT)

            file = open(file, 'r')
            string_lines = file.readlines()[3]
            score[folder] = float( split(string_lines)[1] )

        return score


    def output_values(self, rmsd_aa_wo_if, rmsd_aa_if,  rmsd_ca_wo_if,
                      rmsd_ca_if, identities, score, nb_templates,
                      output_file = None):
        """
        Write result to file.
        
        @param rmsd_aa_wo_if: Rmsd for heavy atoms and normal fit.
                              Data should be a dictionary mapping pdb
                              codes to a list containing the rmsd value
                              and the percent of discharded atoms in
                              the rmsd calculation.
        @type  rmsd_aa_wo_if: {str:[float,float]}
        @param rmsd_aa_if: Rmsd for heavy atoms, iterative fit.
        @type  rmsd_aa_if: {str:[float,float]}
        @param rmsd_ca_wo_if: rmsd for only CA, normal fit.
        @type  rmsd_ca_wo_if: {str:[float,float]}
        @param rmsd_ca_if: Rmsd for only CA, iterative fit.
        @type  rmsd_ca_if: {str:[float,float]}
        @param identities: mean identity to template, dictionary
                           mapping pdb codes to identity values
        @type  identities: {str:float}
        @param score: score calculated by  Modeller
        @type  score: float
        @param nb_templates: number of templates used for re-modeling
        @type  nb_templates: int
        @param output_file: file to write
                            (default: None S{->} outFolder/L{F_OUTPUT_VALUES})
        @type  output_file: str
        """
        output_file = output_file or self.outFolder + self.F_OUTPUT_VALUES

        file = open( output_file, 'w' )
        file.write("PROPERTIES OF RE-MODELED TEMPLATES:\n\n")
        file.write("     | NORMAL FIT  |          ITERATIVE FIT          | IDENTITY  | SCORE | NR\n"  )
        file.write("PDB  | heavy  CA   | heavy percent    CA   percent   | mean to   | mod8  | of\n")
        file.write("code | rmsd   rmsd | rmsd  discarded  rmsd discarded | templates |       | templates\n")

        for key, value in rmsd_aa_wo_if.items():

            file.write("%4s %6.2f %6.2f %6.2f %8.1f %7.2f %7.1f %10.1f %9i %6i\n"%\
                       (key, value[0], rmsd_ca_wo_if[key][0], rmsd_aa_if[key][0],
                        rmsd_aa_if[key][1], rmsd_ca_if[key][0], rmsd_ca_if[key][1],
                        identities[key], score[key], nb_templates))

        file.close()


########################################################################
#########  LOCAL RESULTS: Cross Validation ---- RMSD / Res  

    def get_aln_info(self, output_folder = None):
        """
        Collect alignment information.
        
        @param output_folder: output folder (default: None S{->} outFolder)
        @type  output_folder: str
        
        @return: aln_dictionary, contains information from the alignment
                 between the target and its templates
                 e.g. {'name':'target, 'seq': 'sequence of the target'}
        @rtype: dict
        """
        output_folder = output_folder or self.outFolder

        ci = CI(outFolder=output_folder)

        string_lines = ci.get_lines()
        aln_length = ci.search_length(string_lines)
        aln_dictionnary = ci.get_aln_sequences(string_lines, aln_length)
        aln_dictionnary = ci.get_aln_templates(string_lines, aln_dictionnary,
                                               aln_length)
        aln_dictionnary = ci.identities(aln_dictionnary)

        return aln_dictionnary


    def get_templates_rmsd(self, templates):
        """
        Collect RMSD values between all the templates.

        @param templates: name of the different templates
        @type  templates: [str]
        
        @return: template_rmsd_dic, contains all the rmsd per residues
                 of all the templates
        @rtype: dict
        """
        template_rmsd_dic = {}
        for template in templates:

            pdb_list = self.outFolder + self.F_TEMPLATE_FOLDER \
                       + "/%s"%template + self.F_PDBModels

            pdb_list = T.load(pdb_list)

            template_rmsd_dic[template] = \
              pdb_list[0].compress(pdb_list[0].maskCA()).atoms["rmsd2ref_if"]

        return template_rmsd_dic


    def templates_profiles(self, templates, aln_dic, template_rmsd_dic):
        """
        Collect RMSD profiles of each template with the target and their %ID.
        
        @param templates: name of the different templates
        @type  templates: [str]
        @param aln_dic: contains all the informations between the
                        target and its templates from the alignment
        @type  aln_dic: dict
        @param template_rmsd_dic: contains all the rmsd per residues of all
                                  the templates
        @type  template_rmsd_dic: dict
        
        @return: template_profiles, contains all the profile rmsd of
                 each template with the target and their %ID
        @rtype: dict
        """
        templates_profiles = {}

        target_aln = aln_dic["target"]["seq"]
        for template in templates:
            template_aln = []
            template_profile = []
            template_info = {}

            template_rmsd = template_rmsd_dic[template]
            for key in aln_dic:
                if(key[:4] == template):
                    template_aln = aln_dic[key]["seq"]

            no_res = -1
            for i in range(len(target_aln)):
                if(template_aln[i] is not '-'):
                    no_res += 1

                if(target_aln[i] != '-' and template_aln[i] != '-'):
                    template_profile.append(template_rmsd[no_res])

                if(target_aln[i] != '-' and template_aln[i] == '-'):
                    template_profile.append(-1)

            template_info["rProfile"] = template_profile

            for key in aln_dic["target"]["cov_ID"]:
                if(key[:4] == template):
                    template_info["cov_ID"] = \
                                      aln_dic["target"]["cov_ID"][key]

            templates_profiles[template] = template_info

        return templates_profiles



    def output_cross_val(self, aln_dic, templates_profiles,
                         templates, model, output_file=None):
        """
        Calculates the mean rmsd of the model to the templates and
        write the result to a file.
        
        @param aln_dic: contains all the informations between the
                        target and its templates from the alignment
        @type  aln_dic: dict         
        @param templates_profiles: contains all the profile rmsd of
                                   each template with the target and their %ID
        @type  templates_profiles: dict
        @param templates: name of the different templates
        @type  templates: [str]        
        @param model: model
        @type  model: PDBModel
        @param output_file: output file
                            (default: None S{->} outFolder/L{F_CROSS_VAL})
        @type  output_file: str
        
        @return: mean_rmsd, dictionary with the mean rmsd of the model
                 to the templates.
        @rtype: dict        
        """
        output_file = output_file or self.outFolder + self.F_CROSS_VAL

        mean, sum, values = 0, 0, 0
        mean_rmsd = []

        for k,v in aln_dic["target"]["cov_ID"].items():

            if (k != "target"):
                sum += aln_dic["target"]["cov_ID"][k]
                values +=1

##        cov_id_target = float(sum/values)

        for i in range(len(templates_profiles[templates[0]]["rProfile"])):

            mean = 0
            sum = 0
            n_values = 0

            for k in templates_profiles:

                if(templates_profiles[k]["rProfile"][i] != -1):
                    sum +=  templates_profiles[k]["rProfile"][i]
                    n_values += 1

            if(n_values != 0):        
                mean = float(sum) / float(n_values)

            else: mean = -1

            mean_rmsd.append(mean)

        ## write header
        file = open (output_file, 'w')
        file.write("Mean rmsd of model to templates and the residue rmsd.\n")

        ## write pdb code
        file.write(" "*7)
        for k in templates_profiles.keys():
            file.write("%6s"%k)
        file.write("  mean\n")

        ## write mean rmsd
        file.write(" "*7)
        for k in templates_profiles.keys():
            file.write("%6.2f"%templates_profiles[k]["cov_ID"])
        file.write("\n%s\n"%('='*70))

        ## write rmsd residue profiles
        res_nr = model.compress( model.maskCA()).atoms['residue_number']
        res_na = model.compress( model.maskCA()).atoms['residue_name']        
        for i in range(len(templates_profiles[templates[0]]["rProfile"])):
            file.write("%3i %3s"%(res_nr[i],
                                  res_na[i]))
            for k in templates_profiles:
                file.write("%6.2f"%(templates_profiles[k]["rProfile"][i]))

            file.write("%6.2f\n"%(mean_rmsd[i]))

        file.close()

        return mean_rmsd


#################################################
######   3D Structure: mean RMSD  


    def updatePDBs_charge(self, mean_rmsd_atoms, model):
        """
        pickle down the final.pdb which is judged to be the best model
        of the project. The mean rmsd to the templates is written to the
        temperature_factor column.

        @param mean_rmsd_atoms: mean rmsd for each atom of the
                                target's model
        @type  mean_rmsd_atoms: [int]
        @param model: target's model with the highest modeller score
        @type  model: PDBModel
        """
        model['temperature_factor'] = mean_rmsd_atoms

        model.writePdb(self.outFolder + self.F_FINAL_PDB)       



######################
### LAUNCH FUNCTION ##
######################

    def go(self, output_folder = None, template_folder = None):
        """
        Run analysis of models.

        @param output_folder: folder for result files
                         (default: None S{->} outFolder/L{F_RESULT_FOLDER})
        @type  output_folder: str
        @param template_folder: folder with template structures
                         (default: None S{->} outFolder/L{VS.F_RESULT_FOLDER})
        @type  template_folder: str
        """
        ##
        pdb_list = T.load(self.outFolder + self.F_MODELS)
        model = PDBModel(pdb_list[0])

        ## 
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER
        template_folder = template_folder or self.outFolder +VS.F_RESULT_FOLDER

        templates = self.__listDir(template_folder)

        ##
        global_rmsd_aa_wo_if, global_rmsd_aa_if = self.global_rmsd_aa()
        global_rmsd_ca_wo_if, global_rmsd_ca_if = self.global_rmsd_ca()
        nb_templates = len(templates)-1

        identities = self.get_identities(nb_templates)
        score = self.get_score()

        self.output_values(global_rmsd_aa_wo_if, global_rmsd_aa_if,
                           global_rmsd_ca_wo_if, global_rmsd_ca_if,
                           identities, score, nb_templates)

        ##
        aln_dic = self.get_aln_info(output_folder=self.outFolder)

        template_rmsd_dic = self.get_templates_rmsd(templates)
        templates_profiles = self.templates_profiles(templates,
                                                     aln_dic,
                                                     template_rmsd_dic)
        mean_rmsd = self.output_cross_val(aln_dic, templates_profiles,
                                          templates, model)

        ##
        mean_rmsd_atoms = model.res2atomProfile(mean_rmsd) 
        self.updatePDBs_charge(mean_rmsd_atoms, model)



#############
##  TESTING        
#############
import Biskit.test as BT        

class Test( BT.BiskitTest ):
    """
    Test class
    """

    def prepare(self):
        import tempfile
        import shutil

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_Analyse' )

        ## data from validation sub-projects
        dir = '/Mod/project/validation/'
        for v in ['1DT7', '1J55']:

            os.makedirs( self.outfolder +'/validation/%s/modeller'%v )
            shutil.copy( T.testRoot() +dir+'/%s/modeller/Modeller_Score.out'%v,
                         self.outfolder + '/validation/%s/modeller'%v)

            shutil.copy( T.testRoot() + dir + '/%s/identities_cov.out'%v,
                         self.outfolder + '/validation/%s'%v )

            os.mkdir( self.outfolder +'/validation/%s/benchmark'%v )
            for f in [ 'rmsd_aa.out', 'rmsd_ca.out', 'PDBModels.list' ]:
                shutil.copy( T.testRoot() + dir + '/%s/benchmark/%s'%(v,f),
                             self.outfolder + '/validation/%s/benchmark'%v )


        ## data from main project
        os.mkdir( self.outfolder +'/modeller' )
        shutil.copy( T.testRoot() + '/Mod/project/modeller/PDBModels.list',
                     self.outfolder + '/modeller' )  

        os.mkdir( self.outfolder +'/t_coffee' )
        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     self.outfolder + '/t_coffee' )    


    def test_Analyzer(self):
        """Mod.Analyzer test"""

        self.a = Analyse( outFolder = self.outfolder )
        self.a.go()

        if self.local and self.DEBUG:
            self.log.add(
                'The result from the analysis is in %s/analyse'%self.outfolder)
            
    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )
    

class ProjectAnalyzeTest( BT.BiskitTest ):
    """
    Test case for analyzing a complete modeling project in test/Mod/project
    """

    ## rename to test_Analyze to include this test
    def t_Analyze( self ):
        """
        Mod.Analyze full test/Mod/project test
        """
        self.outfolder = T.testRoot() + '/Mod/project'

        self.a = Analyse( outFolder = self.outfolder )
        self.a.go()

        if self.local:
            print 'The result from the analysis can be found in %s/analyse'%outfolder

if __name__ == '__main__':

    BT.localTest()
