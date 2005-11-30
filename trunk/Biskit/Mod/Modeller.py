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
## original authors: Raik Gruenberg, Johan Leckner
## contributing    : Olivier Perin (post-processing)
##
## last $Author$
## last $Date$
"""
Interface to Modeller
"""

import re, os, os.path, subprocess
import commands
import linecache
from string import *

import settings
import Biskit.tools as tools

from Biskit.Mod.TemplateCleaner import TemplateCleaner as TC
from Biskit.Mod.SequenceSearcher import SequenceSearcher as SS
from Biskit.Mod.CheckIdentities import Check_Identities as CI
from Aligner import Aligner
from Biskit.ModelList import ModelList
from Biskit.DictList import DictList


from Biskit import PDBModel
import glob

from Biskit import StdLog, EHandler


class ModellerError( Exception ):
    pass

class Modeller:
    """
    Take Alignment from t_coffee and template PDBs.
    -> create modeller-script and run modeller
    """

    MODELLER_TEMPLATE = """
# Homology modelling by the MODELLER TOP routine 'model'.

INCLUDE                               # Include the predefined TOP routines

SET OUTPUT_CONTROL = 1 1 1 1 1        # uncomment to produce a large log file
SET ALNFILE  = '%(f_pir)s'            # alignment filename
SET KNOWNS   = %(template_ids)s  
SET SEQUENCE = '%(target_id)s'        # code of the target
SET ATOM_FILES_DIRECTORY = '%(template_folder)s'# directories for input atom files
SET STARTING_MODEL= %(starting_model)i 	# index of the first model 
SET ENDING_MODEL  = %(ending_model)i    # index of the last model
                                   #(determines how many models to calculate)

CALL ROUTINE = 'model'             # do homology modelling
"""


    F_RESULT_FOLDER = '/modeller'
    F_MOD_SCRIPT = 'modeller.top'

    F_INPUT_FOLDER = F_RESULT_FOLDER

    ## default modeller inp file
    F_INP = F_RESULT_FOLDER + '/' + F_MOD_SCRIPT

    ## standard file name output for PDBModels
    F_PDBModels = '/PDBModels.list'

    ## standard file name output for Objective Function
    F_SCORE_OUT = F_RESULT_FOLDER + '/Modeller_Score.out'

    
    def __init__( self, outFolder='.', log=None ):

        self.outFolder = tools.absfile( outFolder )
        self.log = log

        self.prepareFolders()

        self.f_inp = None
        
        ## sequence ids from the last prepared alignment
        self.pir_ids = [] 
        

    def prepareFolders( self ):

        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            self.logWrite( 'Creating '+self.outFolder + self.F_RESULT_FOLDER)
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg


    def create_inp( self, f_pir, target_id, template_folder, template_ids,
                    fout=None, alignment=0,
                    starting_model=1, ending_model=10 ):
        """
        f_pir           - str, PIR alignment file
        target_id       - str, ID of target in alignment
        template_ids    - [ str ], ids of all templates in alignment
        starting_model  - int,
        ending_model    - int,
        -> str, inp file name
        """
        template_ids = [ "'%s'"%id  for id in template_ids ]
        template_ids = ' '.join( template_ids )

    
        params = { 'template_ids' : template_ids,
                   'target_id' : target_id,
                   'f_pir' : f_pir,
                   'template_folder' : template_folder,
                   'starting_model' : starting_model,
                   'ending_model' : ending_model }

        fout = fout or self.outFolder + self.F_INP

        result = self.MODELLER_TEMPLATE % params

        try:
            open( fout, 'w').write( result )
        except IOError, why:
            raise ModellerError( "Can't create modeller inp file "+fout )


    def get_template_ids( self, template_folder ):
        """
        Extract names of all PDB files in a folder (without .PDB)
        template_folder - str, folder
        -> [ str ]
        """
        fs = os.listdir( template_folder )
        return [ f[:-4] for f in fs if f[-4:].upper()=='.PDB' ]


    def get_target_id( self, f_fasta ):
        """
        Extract (first) sequence id from fasta file. Make intelligent guess
        if there is no exact match between target fasta and alignment file.
        f_fasta - str, fasta file with target sequence
        """
        ex_fasta_id = re.compile('^>([A-Za-z0-9_]+)')
        
        title = open( f_fasta ).readline()

        try:
            id_fasta = ex_fasta_id.findall( title )[0]
        except IndexError, why:
            raise ModellerError( "Can't extract sequence id from "+f_fasta )

        if self.pir_ids:
            ## the id matches one previously found in the alignment
            if id_fasta in self.pir_ids:
                return id_fasta

            guesses = [ pid for pid in self.pir_ids
                        if pid[:len(id_fasta)] == id_fasta ]

            if len( guesses ) == 1:
                return guesses[0]

            raise ModellerError("Target ID %s is not found in alignment." %\
                                id_fasta )

        self.logWrite('Warning: target sequence ID extracted without check'+
                      ' against alignment file. Use prepare_alignment first.')
        
        title.close()
        
        return id_fasta
    

    def prepare_alignment( self, f_pir, template_ids ):
        """
        Convert t_coffee - generated PIR into PIR for Modeller.
        On the way all sequence ids are extracted into self.pir_ids.
        f_pir_in        - str, PIR alignment file
        f_pir_out       - str, file to create
        template_ids    - [ str ], ids of all templates in alignment
        -> str, pir alignment file name (for modeller)
        """
        ex_title = re.compile( '^>P1;([A-Za-z0-9_]+)' )
        ex_description = re.compile( '^sequence|^structure' )

        self.pir_ids = []

        try:
            bak_fname = f_pir + '_original'
            os.rename( f_pir,  bak_fname )
        except OSError, why:
            raise ModellerError( "Can't rename " + f_pir + " to "+ bak_fname )

        lines = open( bak_fname ).readlines()

        for i in range( len( lines ) ):
            current = lines[ i ]

            if ex_title.match( current ):
                id = ex_title.findall( current )[0]

                self.pir_ids += [ id ]

                prefix = 'sequence'
                if id in template_ids:
                    prefix = 'structure'
##                     if lines[i + 2][0] != '*':
##                         lines[i + 2] = '*' + lines[i + 2]

                lines[i + 1] = "%s:%s: : : : : : : : \n" \
                               % (prefix, id)

            else:
                if not ex_description.match( current ):
                    lines[ i ] = current.upper()

        open( f_pir, 'w' ).writelines( lines )

            


    def prepare_modeller( self, fasta_target=None, f_pir=None,
                          template_folder=None, fout=None,
                          starting_model=1, ending_model=10 ):
        """
        
        """
        fasta_target = fasta_target or self.outFolder + SS.F_FASTA_TARGET
        template_folder = template_folder or self.outFolder + TC.F_MODELLER
        f_pir = f_pir or self.outFolder + Aligner.F_FINAL_ALN

        template_ids = self.get_template_ids( template_folder )
        
        ## add Modeller-style lines and extract sequence ids
        self.prepare_alignment( f_pir, template_ids )

        ## guess target seq id within alignment
        target_id = self.get_target_id( fasta_target )

        self.create_inp( f_pir, target_id, template_folder, template_ids,
                    fout, starting_model=1, ending_model=10)

        
    def go( self, host=None ):
        """
        """
        try:
            cmd = "cd %s%s ; %s %s"% (self.outFolder, self.F_RESULT_FOLDER, \
                                   settings.modeller_bin, self.F_MOD_SCRIPT)
            
            if host:
                cmd = "ssh %s '%s'" % (host, cmd)
            
            self.logWrite( 'Running Modeller .. ')
            self.logWrite( cmd )

            sp = subprocess.Popen(cmd, shell=True, env=os.environ)
            sp.wait()

            ##self.logWrite( '..done: ' + str( output ) + '\n' )

        except EnvironmentError, why:
            self.logWrite("ERROR: Can't run Modeller: "+ str( why ) )
            raise AlignerError( "Can't run Modeller: " + str( why ) )


    ####################
    ## Post processing

    def pdb_list(self, model_folder):
        """
        model_folder - str, folder containing model PDBs (output from Modeller)
        -> [str], pdb files from modeller
        """
        return glob.glob('%s/target.B*'%(model_folder))
          

    def extract_modeller_score(self, f_model, model):
        """
        Extract the modeller score from the models
        f_model - str, file name of model pdb
        """
        r1 = re.compile(r'\w\d.*\d')
        
        line_of_interest = linecache.getline(f_model, 2)
        OF = r1.findall(line_of_interest)

       
        file = open(model.validSource(),'r')
        string = file.readlines()[:2]
        model.info["headlines"] = string
        file.close()

        return  float(OF[0])


    def output_score(self, pdb_list, model_folder):
        """
        pdb_list - ModelList
        model_folder - str, ouput folder
        """
        file_output = open(self.outFolder+self.F_SCORE_OUT,'w')
        file_output.write("The models from modeller with their Score:\n")

        for model in pdb_list:
        
            file_output.write('%s\t%6.2f\n'%(tools.stripFilename(
                model.validSource()), model.info["mod_score"]))
    
        file_output.close()


    def update_PDB(self, pdb_list, cwd, model_folder):
        """
        Extract number of templates per residue position from alignment.
        Write new PDBs with number of templates in occupancy column.
        """

        cwd = cwd or self.outFolder

        ci = CI(cwd)
        aln_dictionnary = ci.go()

        template_info = aln_dictionnary["target"]["template_info"]

        for model in pdb_list:

            model.setResProfile("n_templates", template_info,
                                comment="number of templates for each residue")

        rmask = pdb_list[0].profile2mask("n_templates", 1,1000)
        amask = pdb_list[0].res2atomMask(rmask)



        for i in range(len(pdb_list)):

            for a,v in zip(pdb_list[i].atoms, amask):
                a['occupancy'] = v

            t = []
            
            for string in pdb_list[i].info["headlines"]:

                tuple = ()
                
                tuple = string.split()[0],join(string.split()[1:])
                t.append(tuple)

    
            pdb_list[i].writePdb('%s/model_%02i.pdb'%(model_folder,i),
                                 headlines = t)


    def update_rProfiles(self, pdb_list):
        """
        pdb_list - ModelList
        """
        for model in pdb_list:

            m = model.compress( model.maskCA())

            atoms = DictList( m.atoms )
            prof = atoms.valuesOf( "temperature_factor")

            model.setResProfile('mod_score', prof)
           

            
    def write_PDBModels(self, pdb_list, model_folder = None):
        """
        Dump the list of PDBModels.
        pdb_list - ModelList
        model_folder - str, ouput folder
        """    
        tools.Dump(pdb_list, '%s'%(model_folder + self.F_PDBModels))


    def postProcess(self, model_folder=None):
        """
        Rename model PDBs to 'model_01.pdb' - 'model_xx.pdb' according to their
        score, create and dump a ModelList instance containing these models
        and info-records and profiles with the Modeller score.
        model_folder - str, folder containing model PDBs ['./modeller']
        """
        
        model_folder = model_folder or self.outFolder + self.F_INPUT_FOLDER

        model_files = self.pdb_list(model_folder)
        
        pdb_list = ModelList( model_files )
        
        for model in pdb_list:

            model.info["mod_score"] = self.extract_modeller_score(
                model.validSource(), model)

        pdb_list = pdb_list.sortBy("mod_score")

        self.update_PDB(pdb_list, self.outFolder, model_folder)

        self.output_score(pdb_list,model_folder)

        self.update_rProfiles(pdb_list)

        self.write_PDBModels(pdb_list, model_folder)
    


##########
## TEST ##
##########
if __name__ == '__main__':

    m = Modeller( tools.testRoot() + '/Mod/project')

    r = m.prepare_modeller( )

    m.go()

    m.postProcess()
