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
## $Revision$
"""
Interface to Modeller
"""

import re, os, os.path, subprocess
import commands
import linecache
from string import *

import settings
import Biskit.tools as T

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
    Creates a modeller-script and runs modeller.
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
        """
        @param outFolder: base folder for Modeller output 
                          (default: L{F_RESULT_FOLDER})
        @type  outFolder: str
        @param log: log file instance, if None, STDOUT is used (default: None)
        @type  log: LogFile        
        """
        self.outFolder = T.absfile( outFolder )
        self.log = log

        self.prepareFolders()

        self.f_inp = None

        ## sequence ids from the last prepared alignment
        self.pir_ids = [] 


    def prepareFolders( self ):
        """
        Create needed output folders if they don't exist.
        """
        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            self.logWrite( 'Creating '+self.outFolder + self.F_RESULT_FOLDER)
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


    def create_inp( self, f_pir, target_id, template_folder, template_ids,
                    fout=None, alignment=0,
                    starting_model=1, ending_model=10 ):
        """
        Create a input (.top) file for Modeller.
        See L{MODELLER_TEMPLATE}.
        
        @param f_pir: PIR alignment file
        @type  f_pir: str
        @param target_id: ID of target in the alignment file
        @type  target_id: str
        @param template_folder: folder with template coordinate files
        @type  template_folder: str
        @param template_ids: ids of all templates in the alignment file
        @type  template_ids: [str]
        @param fout: name of modeller script (default: None -> L{F_INP})
        @type  fout: None or str
        @param starting_model: first model
        @type  starting_model: int
        @param ending_model: last model
        @type  ending_model: int
        
        @raise ModellerError: if can't create modeller inp file
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
        
        @param template_folder: folder with template coordinate files
        @type  template_folder: str

        @return: list of PDB names
        @rtype: [str]      
        """
        fs = os.listdir( template_folder )
        return [ f[:-4] for f in fs if f[-4:].upper()=='.PDB' ]


    def get_target_id( self, f_fasta ):
        """
        Extract (first) sequence id from fasta file. Make intelligent guess
        if there is no exact match between target fasta and alignment file.
        
        @param f_fasta: fasta file with target sequence
        @type  f_fasta: str

        @raise ModellerError: if target ID is not found in alignment
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

        @param f_pir: PIR alignment file
        @type  f_pir: str
        @param template_ids: file to create
        @type  template_ids: str

        @return: pir alignment file name (for modeller)
        @rtype: str

        @raise ModellerError: if can't rename file
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
        @param fasta_target: target fasta file
                             (default: None -> L{SS.F_FASTA_TARGET})
        @type  fasta_target: str
        @param f_pir: PIR alignment file
                      (default: None -> L{Aligner.F_FINAL_ALN})
        @type  f_pir: str
        @param template_folder: folder with template coordinate files
                                (default: None -> L{TC.F_MODELLER})
        @type  template_folder: str
        @param fout: name of modeller script (default: None -> L{F_INP})
        @type  fout: None or str
        @param starting_model: first model (default: 1)
        @type  starting_model: int
        @param ending_model: last model (default: 10)
        @type  ending_model: int
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
        Run Modeller job.

        @param host:  name of host to use (default: None -> localhost)
        @type  host: str

        @raise ModellerError: if Modeller job fails
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
            raise ModellerError( "Can't run Modeller: " + str( why ) )


    ####################
    ## Post processing

    def pdb_list(self, model_folder):
        """
        Compile a list with the names of all models created by Modeller.
        
        @param model_folder: folder containing model PDBs
                             (output from Modeller)
        @type  model_folder: str

        @return: pdb files from modeller
        @rtype: [str]
        """
        return glob.glob('%s/target.B*'%(model_folder))


    def extract_modeller_score(self, f_model, model):
        """
        Extract the modeller score from the models

        @param f_model: file name of model pdb
        @type  f_model: str
        @param model: model
        @type  model: PDBModel
        
        @return: modeller score for the models
        @rtype: float
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
        Write a list of model scores to file L{F_SCORE_OUT}
        
        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param model_folder: ouput folder
        @type  model_folder: str
        """
        file_output = open(self.outFolder+self.F_SCORE_OUT,'w')
        file_output.write("The models from modeller with their Score:\n")

        for model in pdb_list:

            file_output.write('%s\t%6.2f\n'%(T.stripFilename(
                model.validSource()), model.info["mod_score"]))

        file_output.close()


    def update_PDB(self, pdb_list, cwd, model_folder):
        """
        Extract number of templates per residue position from alignment.
        Write new PDBs with number of templates in occupancy column.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param cwd: current project directory
        @type  cwd: str
        @param model_folder: folder with models
        @type  model_folder: str
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
        Read Modeller score from the 'temperature_factor' column and
        create a profile from it.
        
        @param pdb_list: list of models
        @type  pdb_list: ModelList
        """
        for model in pdb_list:

            m = model.compress( model.maskCA())

            atoms = DictList( m.atoms )
            prof = atoms.valuesOf( "temperature_factor")

            model.setResProfile('mod_score', prof)



    def write_PDBModels(self, pdb_list, model_folder = None):
        """
        Dump the list of PDBModels.
        
        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param model_folder: ouput folder (default: None -> L{F_PDBModels})
        @type  model_folder: str
        """    
        T.Dump(pdb_list, '%s'%(model_folder + self.F_PDBModels))


    def postProcess(self, model_folder=None):
        """
        Rename model PDBs to 'model_01.pdb' - 'model_xx.pdb' according to their
        score, create and dump a ModelList instance containing these models
        and info-records and profiles with the Modeller score.

        @param model_folder: folder containing model PDBs
                             (default: None -> L{F_INPUT_FOLDER})
        @type  model_folder: str
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



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self ):
        """
        run function test

        @return: 1
        @rtype:  int
        """
        import tempfile
        import shutil

        ## collect the input files needed
        outfolder = tempfile.mkdtemp( '_test_Modeller' )
        os.mkdir( outfolder +'/templates' )
        os.mkdir( outfolder +'/t_coffee' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/modeller',
                         outfolder + '/templates/modeller' )

        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     outfolder + '/t_coffee' )    

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     outfolder  )


        m = Modeller( outfolder )

        r = m.prepare_modeller( )

        m.go()

        m.postProcess()

        print 'The modelling result can be found in %s/modeller'%outfolder


        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    

if __name__ == '__main__':

    test = Test()
    
    assert test.run() ==  test.expected_result()
