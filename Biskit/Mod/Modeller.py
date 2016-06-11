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
## original authors: Raik Gruenberg, Johan Leckner
## contributing    : Olivier Perin (post-processing)
##
"""
Interface to Modeller
"""

import re, os.path
import linecache
import glob

import settings
import Biskit.tools as T

from Biskit.Mod.TemplateCleaner import TemplateCleaner as TC
from Biskit.Mod.SequenceSearcher import SequenceSearcher as SS
from Biskit.Mod.CheckIdentities import CheckIdentities
from Aligner import Aligner
from Biskit.Mod.TemplateFilter import TemplateFilter
from Biskit.ModelList import ModelList
from Biskit.DictList import DictList
from Biskit.Executor import Executor
from Biskit import PDBModel


class ModellerError( Exception ):
    pass

class Modeller( Executor ):
    """
    Interface to the Modeller program for building homology models:

    1) Read: Alignment from t_coffee and template PDBs.

       - t_coffee/final.pir_aln
       - templates/modeller/*pdb

    2) Create a modeller-script -> modeller/modeller.top
    3) Run modeller.
    4) Postprocess the modeller models.

    Result:

    - models sorted by their modeller score into model_00.pdb - model_xx.pdb:
       - the local modeller score is in the temperature_factor column
       - the number of templates for each residue is in the occupancy column

    - a pickled list of corresponding PDBModels:
       - modeller score in residue profile 'mod_score'
       - template coverage is in residue profile 'n_templates'
       - info['mod_score'] contains the overall modeller score
       - info['mod_pdb'] contains the file name of the raw modeller PDB 
    """

    ## non-.py ending to prevent byte-compilation of template
    MODEL_SCRIPT      = 'model_mult.py_'  
    MODEL_SCRIPT_PATH = T.dataRoot() + '/modeller/' + MODEL_SCRIPT

    F_RESULT_FOLDER = '/modeller'
##     F_MOD_SCRIPT = 'modeller.top'

    F_INPUT_FOLDER = F_RESULT_FOLDER

    #: standard file name output for PDBModels
    F_PDBModels = '/PDBModels.list'

    #: standard file name output for Objective Function
    F_SCORE_OUT = F_RESULT_FOLDER + '/Modeller_Score.out'


    def __init__( self, outFolder='.', mod_template=None,
                  fasta_target=None, template_folder=None, f_pir=None,
                  starting_model=1, ending_model=10,
                  zfilter=None, idfilter=None,
                  disc_report=1,
                  **kw ):
        """
        Create Executor instance for one modeller run.

        @param outFolder: base folder for Modeller output 
                          (default: L{F_RESULT_FOLDER})
        @type  outFolder: str
        @param mod_template   : template file/str for modeller input script
        @type  mod_template   : str
        @param fasta_target   : file with target sequence
        @type  fasta_target   : str
        @param template_folder: folder containing template PDBs for modeller
        @type  template_folder: str
        @param f_pir         : file containing alignment
        @type  f_pir         : str
        @param starting_model: first homology model to report [1]
        @type  starting_model: int
        @param ending_model  : last model to report [10]
        @type  ending_model  : int
        @param zFilter: override z-score cutoff of Mod.TemplateFilter
        0..skip this filter
        @type  zfilter: float
        @param idfilter: override sequence identity cutoff of TemplateFilter
        0..skip this filter
        @type  idFilter: float
        @param disc_report : keep logs, write ordered pdbs to model_xx.pdb [1]
        @type  disc_report : bool
        """
        self.outFolder = T.absfile( outFolder )
        cwd = self.outFolder + self.F_RESULT_FOLDER

        mod_template = mod_template or self.MODEL_SCRIPT_PATH

        # don't override existing modeller script
        f_in = os.path.join( cwd, self.MODEL_SCRIPT )
        if mod_template and os.path.isfile( mod_template ) or \
           os.path.islink( mod_template ):
            f_in = os.path.join( cwd, os.path.basename( mod_template ) )

        Executor.__init__(self, 'modeller',
                          template=mod_template,
                          args= f_in, # command line argument
                          f_in= f_in, # target of template completion
                          cwd= cwd,
                          **kw)

        self.fasta_target= fasta_target or self.outFolder + SS.F_FASTA_TARGET
        self.template_folder= template_folder or self.outFolder + TC.F_MODELLER
        self.f_pir = f_pir or self.outFolder + Aligner.F_FINAL_ALN

        self.target_id = None   # will be assigned in prepare()
        self.starting_model = starting_model
        self.ending_model   = ending_model

        #: prepare() fills in the template ids formatted for old top or py
        self.knowns_top = None
        self.knowns_py  = None 

        self.aln_info = CheckIdentities( self.outFolder )
        self.z_filter = zfilter
        self.id_filter = idfilter

        self.disc_report = disc_report

        ## sequence ids from the last prepared alignment
        self.pir_ids = []


    def prepareFolders( self ):
        """
        Create needed output folders if they don't exist.
        """
        if not os.path.exists( self.outFolder + self.F_RESULT_FOLDER ):
            os.mkdir( self.outFolder + self.F_RESULT_FOLDER )
            if self.verbose:
                self.log.add('Creating '+self.outFolder + self.F_RESULT_FOLDER)


    def generateInp( self ):
        """
        Modify Executor.generateInp so that it does NOT override an existing
        input file.
        """
        if os.path.exists( self.f_in ):
            self.log.add('Warning: Using EXISTING modeller input file %s' %\
                         self.f_in )
            self.keep_inp = 1
        else:
            Executor.generateInp( self )


    def get_template_ids( self, template_folder ):
        """
        Extract names of all PDB files in a folder (without .PDB)

        @param template_folder: folder with template coordinate files
        @type  template_folder: str

        @return: list of PDB names
        @rtype: [str]      
        """
        fs = os.listdir( template_folder )
        r = [ f[:-4] for f in fs if f[-4:].upper()=='.PDB' ]

        return r


    def filter_templates( self, target_id='target' ):
        """
        Kick out some templates before the modeller run.
        """

        tf = TemplateFilter( self.aln_info, target_id=target_id,
                             verbose=self.verbose, log=self.log )

        if self.z_filter != 0:
            tf.filter_z( self.z_filter )

        if self.id_filter != 0:
            tf.filter_id( self.id_filter )

        r = tf.get_filtered()

        return r


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

        self.log.add('Warning: target sequence ID extracted without check'+
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




    def prepare( self ):
        """
        Executed before modeller run.
        """
        self.prepareFolders()

        template_ids = self.get_template_ids( self.template_folder )

        ## add Modeller-style lines and extract sequence ids
        self.prepare_alignment( self.f_pir, template_ids )

        ## analyze alignment
        self.aln_info.go()
        ## write identity matrix files into project folder
        self.aln_info.write_identities()

        ## guess target seq id within alignment
        self.target_id = self.get_target_id( self.fasta_target )

        ## remove templates with low or below average similarity to target 
        template_ids = self.filter_templates( target_id=self.target_id )

        ## convert template ids into modeller-formatted input line
        ## for python or old-style top input
        self.knowns_py  = [ "'%s'"%id  for id in template_ids ]
        self.knowns_py  = ', '.join( self.knowns_py )
        self.knowns_top = ' '.join( template_ids )


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
        score = r1.findall(line_of_interest)

        f = open(model.validSource(),'r')
        s = f.readlines()[:2]
        model.info["headlines"] = s
        f.close()

        return  float(score[0])


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


    def update_models(self, pdb_list ):
        """
        Extract number of templates per residue position from alignment
        and put it into a residue profile 'n_templates' of each model.
        Then put the mask of 0 (no template) or 1 (any template) into
        the occupancy column.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param model_folder: folder with models
        @type  model_folder: str
        """
        template_info =  self.aln_info.result["target"]["template_info"]

        for model in pdb_list:

            model.residues.set("n_templates", template_info,
                               comment="number of templates for each residue")

        rmask = pdb_list[0].profile2mask("n_templates", 1,1000)
        amask = pdb_list[0].res2atomMask(rmask)

        for m in pdb_list:
            m['occupancy'] = amask


    def write_Pdbs( self, pdb_list, model_folder):
        """
        Save re-ordered homology models as model_00.pdb - model_xx.pdb.
        The B-factor column contains the modeller score, occupancy
        contains a mask with 0 meaning there was no template for this position.

        @param pdb_list: list of PDBModels
        @type  pdb_list: ModelList
        @param model_folder: folder to save in
        @type  model_folder: str
        """

        for i, m in enumerate( pdb_list ):

            t = []

            for string in m.info["headlines"]:

                pair = string.split()[0], ''.join(string.split()[1:])
                t.append(pair)

            m.writePdb('%s/model_%02i.pdb'%(model_folder,i), headlines = t)


    def update_rProfiles(self, pdb_list):
        """
        Read Modeller score from the 'temperature_factor' column and
        create a profile from it.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        """
        for model in pdb_list:

            m = model.compress( model.maskCA())

            model.residues.set('mod_score', m['temperature_factor'],
                               comment='local modeller score')


    def write_PDBModelList(self, pdb_list, model_folder = None):
        """
        Dump the list of PDBModels.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param model_folder: ouput folder (default: None -> L{F_PDBModels})
        @type  model_folder: str
        """
        ## cut connection to target.Bxxxxx.pdb
        for m in pdb_list:
            m.disconnect()

        T.dump(pdb_list, '%s'%(model_folder + self.F_PDBModels))


    def postProcess( self ):
        """
        @todo: parse modeller log into this instance?
        """
        pass


    def isFailed( self ):
        """
        Detect wether the modeller process failed to finish. (overrides
        Executor method)
        """
        r = self.pdb_list( self.outFolder + self.F_INPUT_FOLDER )
        return len( r ) < 1

    def cleanup( self ):
        Executor.cleanup( self )

        if not self.disc_report:
            T.tryRemove( self.outfolder + self.F_RESULT_FOLDER, tree=1 )


    def finish(self):
        """
        If self.disc_report == 1:
        Rename model PDBs to 'model_01.pdb' - 'model_xx.pdb' according to their
        score, create and dump a ModelList instance containing these models
        and info-records and profiles with the Modeller score.

        Otherwise the resulting models are only put into self.result and
        nothing is written to disc.

        @param model_folder: folder containing model PDBs
                             (default: None -> L{F_INPUT_FOLDER})
        @type  model_folder: str
        """
        model_folder = self.outFolder + self.F_RESULT_FOLDER

        model_files = self.pdb_list(model_folder)

        self.result = ModelList( model_files )

        for model in self.result:

            model.info["mod_score"] = self.extract_modeller_score(
                model.validSource(), model)
            model.info['mod_pdb'] = model.source

        self.result = self.result.sortBy("mod_score")

        self.update_models(self.result )

        if self.disc_report:

            self.output_score(self.result,model_folder)

            self.update_rProfiles(self.result)

            self.write_Pdbs( self.result, model_folder )

            self.write_PDBModelList(self.result, model_folder)



#############
##  TESTING        
#############
import Biskit.test as BT

class TestBase(BT.BiskitTest):
    """
    Test base class
    """

    def prepare(self):
        import tempfile
        import shutil

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_Modeller' )
        os.mkdir( self.outfolder +'/templates' )
        os.mkdir( self.outfolder +'/t_coffee' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/modeller',
                         self.outfolder + '/templates/modeller' )

        shutil.copy( T.testRoot() + '/Mod/project/t_coffee/final.pir_aln',
                     self.outfolder + '/t_coffee' )    

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     self.outfolder  )


    def runmodeller( self, run=0 ):
        """
        @param run: run the full test (call external application) or not
        @type  run: 1|0
        """

        self.m = Modeller( self.outfolder, verbose=self.local,
                           ending_model=2 )

##         self.m.prepare_modeller( )

        if run:

            self.m.run()
##             self.m.postProcess()

            if self.DEBUG and self.VERBOSITY > 2:
                self.log.add('The modelling result id in %s/modeller'\
                             %self.outfolder)

            result = T.load( self.outfolder + '/modeller/PDBModels.list' )
            self.assertEqual( len(result), 2 )

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )

class TestDry( TestBase ):
    """Incomplete test case without running Modeller"""

    TAGS = [BT.EXE]

    def test_modeller(self):
        """Mod.Modeller dry run"""
        self.runmodeller(run=0)

class Test( TestBase ):
    """Complete test case running modeller"""

    TAGS = [BT.EXE, BT.LONG]

    def test_modeller(self):
        """Mod.Modeller complete test (ignore the 'import site' error!)"""
        self.runmodeller(run=1)


if __name__ == '__main__':

    BT.localTest(verbosity=3, debug=1)
