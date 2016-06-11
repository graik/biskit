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
"""
Sequence Alignment
"""

import re, os

from TemplateCleaner import TemplateCleaner as TC
from SequenceSearcher import SequenceSearcher as SS

import Biskit.tools as T
from Biskit.Errors import *

from Biskit import Executor, TemplateError


class AlignerError( BiskitError ):
    pass


class TCoffee( Executor ):
    """
    Execute a t_coffee command
    """

    def __init__(self, cmd, host=None, **kw ):
        """
        @param cmd: command sequence
        @type  cmd: str
        @param host: host name for remote execution
                     (default: None=local execution)
        @type  host: str
        """
        Executor.__init__( self, 't_coffee', cmd, catch_out=1,
                           node=host, **kw )


class Aligner:
    """
    Take template CA traces and template sequences plus non-template
    sequences.

    Returns a T-Coffee script or T-Coffee alignment
    """

    F_RESULT_FOLDER = '/t_coffee'

    F_FAST_LIB   = F_RESULT_FOLDER + '/fast_pair.lib'
    F_LALIGN_LIB = F_RESULT_FOLDER + '/lalign_id_pair.lib'
    F_SAP_LIB    = F_RESULT_FOLDER + '/sap_pair.lib'
    F_STRUCT_ALN = F_RESULT_FOLDER + '/struct.aln'
    F_FINAL      = F_RESULT_FOLDER + '/final'
    F_COFFEE_LOG = F_RESULT_FOLDER + '/t_coffee.log'
    F_COFFEE_INP = F_RESULT_FOLDER + '/t_coffee.inp'
    F_FAST_TREE  = F_RESULT_FOLDER + '/fast_tree.dnd'
    F_SEQ_TREE   = F_RESULT_FOLDER + '/sequence.dnd'

    ## default result alignment file. Not actually needed in Aligner itself.
    F_FINAL_ALN  = F_FINAL + '.pir_aln'

    def __init__( self, outFolder='.', log=None, verbose=1, sap=1,
                  debug=False):
        """
        @param outFolder: base folder for t_coffee output
                          (default: L{F_RESULT_FOLDER})
        @type  outFolder: str
        @param log: log file instance, if None, STDOUT is used (default: None)
        @type  log: LogFile
        @param verbose: be verbose
        @type  verbose: 1 | 0
        @param sap: perform structural alignment
        @type  sap: 1 | 0
        """
        self.log = log
        self.outFolder = T.absfile( outFolder )
        self.verbose = verbose
        self.sap = sap
        self.debug = debug

        ## recognize file types for adding t_coffee type code
        self.ex_inp_types = { 'P' : re.compile('.+\.[Pp][Dd][Bb]$|'+
                                               '.+\.[Aa][Ll][Pp][Hh][Aa]$'),
                              'L' : re.compile('.+\.[Ll][Ii][Bb]$'),
                              'S' : re.compile('.*\.[Ff][Aa][Ss][Tt][Aa]$')
                              }

        self.commands = []  ## list of t_coffee commands to run

        fn = self.outFolder + self.F_RESULT_FOLDER
        if not os.path.exists( fn ):
            os.mkdir( fn )
            if verbose: self.logWrite( 'Creating new output folder '+ fn )


    def logWrite( self, msg, force=1 ):
        """
        Write message to log.

        @param msg: message to print
        @type  msg: str
        """
        if self.log:
            self.log.writeln( msg )
        else:
            if force:
                print msg


    def __add_type_code( self, inp ):
        """
        Try guessing a t_coffee type code for file name and add type code
        (P for structurem S for sequence and L for library).

        @param inp: file name
        @type  inp: str

        @return: type_code + file name OR file name unchanged
        @rtype: str
        """
        for code, pattern in self.ex_inp_types.items():
            if pattern.match( inp ):
                return code + inp
        return inp


    def coffee_align_inp( self, input, method=None, out_lib=None,
                          output=None, outfile=None, **kw ):
        """
        Create a single t_coffee command.

        @param input: list of input files
        @type  input: [str]
        @param method: t_coffee method ('M' is added automatically)
                       (default:None)
        @type  method: str
        @param output: output format(s) (default:None)
        @type  output: [str]
        @param out_lib: output library file if needed (default:None)
        @type  out_lib: str
        @param outfile: output file (base) if needed (default:None)
        @type  outfile: str
        @param kw: additional param=value pairs
        @type  kw: param=value

        @return: t_coffee command
        @rtype: str
        """
        r = ''

        input = [ self.__add_type_code( i ) for i in T.toList(input) ]

        if input or method:
            r += ' -in'
            for i in input:
                r += ' ' + i
            if method:
                r += ' M' + method

        if out_lib:
            r += ' -out_lib ' + out_lib

        if output:
            r += ' -output'
            for o in T.toList( output ):
                r += ' ' + o

            if not outfile:
                raise AlignerError, 'T-Coffee -output also requires -outfile'

        if outfile:
            r += ' -outfile ' + outfile

        if not output:
            r += ' -lib_only'

        r += ' -template_file No'  ## MSAP runs generate unwanted *.template_file

        for param, value in kw.items():
            r += ' -'+param + ' ' + str( value )

        return r


    def __default_input_files( self, pdbFiles=None, fasta_templates=None,
                               fasta_sequences=None, fasta_target=None ):
        """
        Fetch missing fasta and pdb files from standard locations.

        @param pdbFiles: fetch pdb file from L{TC.F_COFFEE} (default:None)
        @type  pdbFiles: [str]
        @param fasta_templates: fetch fasta template file from L{TC.F_FASTA}
                                (default:None)
        @type  fasta_templates: [str]
        @param fasta_sequences: fetch fasta sequence file from
                                L{SS.F_FASTA_NR} (default:None)
        @type  fasta_sequences: [str]
        @param fasta_target: fetch fasta target file from
                             L{SS.F_FASTA_TARGET} (default:None)
        @type  fasta_target: [str]

        @return: dictionary mapping input file class to path(s)
        @rtype:  {str:[str]}
        """
        r = {}

        cwd = TCoffee( '' ).cwd or os.getcwd()

        ## collect CA trace PDBs from default location
        if not pdbFiles:
            folder = self.outFolder + TC.F_COFFEE
            fs = os.listdir( folder )
            pdbFiles= [ folder + f for f in fs if f[-6:].upper()=='.ALPHA' ]

            ## relativize to current working folder to reduce length
            ## of T-Coffee / SAP arguments
            pdbFiles= [ T.relpath( cwd, f ) for f in pdbFiles ]

        r['pdbs'] = T.toList( pdbFiles )
        r['templates'] = fasta_templates or self.outFolder + TC.F_FASTA
        r['sequences'] = fasta_sequences or self.outFolder + SS.F_FASTA_NR
        r['target']    = fasta_target or self.outFolder + SS.F_FASTA_TARGET

        return r


    def align_for_modeller_inp( self, pdbFiles=None, fasta_templates=None,
                                fasta_sequences=None, fasta_target=None,
                                f_fast_tree=None, f_sequence_tree=None ):
        """
        Prepare alignment commands for homology modeling.

        @note: If verbose==1, the commands are mirrored to t_coffee.inp.

        @param pdbFiles: template PDBs from L{TC.F_COFFEE} (default:None)
        @type  pdbFiles: [str]
        @param fasta_templates: template sequences from
                                templates/templates.fasta (default:None)
        @type  fasta_templates: [str]
        @param fasta_sequences: more sequences from L{SS.F_FASTA_NR}
                                (default:None)
        @type  fasta_sequences: [str]
        @param fasta_target: fasta target files
                             L{SS.F_FASTA_TARGET} (default:None)
        @type  fasta_target: str
        @param f_fast_tree: filename of clustalW style guide tree (.dnd)
        @type  f_fast_tree: str
        @param f_sequence_tree: filename of clustalW style guide tree (.dnd)
        @type  f_sequence_tree: str    
        """
        ## Fetch default files where needed
        d = self.__default_input_files( pdbFiles, fasta_templates,
                                        fasta_sequences, fasta_target )

        pdbFiles = d['pdbs']
        f_templates = d['templates']
        f_sequences = d['sequences']
        f_target    = d['target']

##         pdbFiles = [ T.absfile( f ) for f in pdbFiles ]

        f_templates = T.absfile( f_templates )
        f_sequences = T.absfile( f_sequences )

        ## create output file names
        f_fast_lib = self.outFolder + self.F_FAST_LIB
        f_lalign_lib = self.outFolder + self.F_LALIGN_LIB
        f_sap_lib = self.outFolder + self.F_SAP_LIB
        f_struct_aln = self.outFolder + self.F_STRUCT_ALN
        f_final_aln = self.outFolder + self.F_FINAL
        f_coffee_log= self.outFolder + self.F_COFFEE_LOG
        f_fast_tree= f_fast_tree or self.outFolder + self.F_FAST_TREE
        f_sequence_tree= f_sequence_tree or self.outFolder + self.F_SEQ_TREE

        ## fix target file
        self.repair_target_fasta( f_target )

        if len( pdbFiles ) == 0:
            raise AlignerError( "\n No templates avaliable. ABORTING" )

        ## SAP will not run if there is only one template
        if (not self.sap) or ( len( pdbFiles ) == 1 ):
            if self.sap:
                self.logWrite('WARNING! Only one template avaliable:' +\
                              str(pdbFiles) )
            self.logWrite('Structural alignment (SAP) will not be performed.')

            r = [
                ## fast global pair-wise alignment
                ## why not use slow pair? better for distant sequences
                self.coffee_align_inp( [ f_sequences, f_templates, f_target ],
                                       method='fast_pair',
                                       out_lib=f_fast_lib,
                                       quiet=f_coffee_log+'_1', convert='' ),

                ## uses LALIGN to find the best pair-wise local alignments
                self.coffee_align_inp( [ f_sequences, f_templates, f_target ],
                                       method='lalign_id_pair',
                                       out_lib=f_lalign_lib,
                                       quiet=f_coffee_log+'_2', convert='' ),

                ## combine alignments to one
                self.coffee_align_inp( [ f_fast_lib, f_lalign_lib],
                                       clean_aln=0, newtree=f_fast_tree,
                                       output=['clustalw','phylip',
                                               'score_html', 'pir_aln'],
##                                     run_name=T.stripFileName(f_final_aln),
outfile=f_final_aln,
quiet=f_coffee_log+'_4')
            ]

        ## Normal alignment run with structural alignment
        ## (more than one template)
        else:
            r = [
                ## fast global pair-wise alignment
                ## why not use slow pair? better for distant sequences
                self.coffee_align_inp( [ f_sequences, f_templates, f_target ],
                                       method='fast_pair',
                                       out_lib=f_fast_lib,
                                       quiet=f_coffee_log+'_1', convert='' ),

                ## uses LALIGN to find the best pair-wise local alignments
                self.coffee_align_inp( [ f_sequences, f_templates, f_target ],
                                       method='lalign_id_pair',
                                       out_lib=f_lalign_lib,
                                       quiet=f_coffee_log+'_2', convert='' ),

                ## uses SAP to do pairwise structural alignments
                self.coffee_align_inp( pdbFiles, newtree=f_sequence_tree,
                                       method='sap_pair',
                                       weight=1000,
                                       out_lib=f_sap_lib,
                                       output='fasta_aln',
                                       outfile=f_struct_aln,
                                       quiet=f_coffee_log+'_3' ),

                ## fix lib file
                (f_sap_lib, f_struct_aln),

                ## combine alignments to one
                self.coffee_align_inp( [ f_fast_lib, f_lalign_lib, f_sap_lib],
                                       clean_aln=0, newtree=f_fast_tree,
                                       output=['clustalw','phylip',
                                               'score_html',
                                               'pir_aln', 'score_ascii'],
                                       ##                                    run_name=T.stripFileName(f_final_aln),
                                       outfile=f_final_aln,
                                       quiet=f_coffee_log+'_4')
            ]

        ## add to internal command 'queue'
        self.commands += r

        if self.verbose:
            f = open( self.outFolder + self.F_COFFEE_INP, 'w' )
            for cmd in r:
                f.write( str(cmd) + '\n\n' )
            f.close()


    def repair_lib_file(self, fname ):
        """
        Msap_pair alignments mess up sequence ids of certain
        (single-chain?) structures.
        @note: Dirty hack..

        @param fname: Msap_pair file to repair
        @type  fname: str
        """
        ex_sapfix = re.compile('_[A-Z]{0,1}_[A-Z]')
        ex_alpha = re.compile('.alpha') ## new from T-Coffe ver. 3.32
        bak_fname = fname+'_original'

        os.rename( fname,  bak_fname)

        fout = open( fname, 'w' )
        for l in open( bak_fname ):
            l = re.sub( ex_alpha, '', l )
            fout.write( re.sub( ex_sapfix, '_', l ) )

        fout.close()


    def repair_target_fasta( self, fname ):
        """
        Change target ID to 'target' and remove '|' etc. that mess up
        things in the alignments.

        @param fname: fasta file to repair
        @type  fname: str

        @raise AlignerError: if cannot fix target sequence file
        """
        bak_fname = fname+'_original'
        os.rename( fname,  bak_fname)

        try:
            fout = open( fname, 'w' )
            fin = open( bak_fname )

            l0 = fin.readline() #skipp header line

            fout.write( '>target\n' )

            for l in fin:
                fout.write( l )

            fin.close()
            fout.close()

        except Exception, why:
            print T.lastError()
            raise AlignerError("Cannot fix target sequence file: "+str(why))


    def go( self, host=None ):
        """
        Run the previously added commands, delete internal command list.

        @param host: host name for remote execution
                     (default: None=local execution)
        @type  host: str

        @raise AlignerError: if T_Coffee execution failed
        """
        try:
            if self.verbose:
                self.logWrite('\nALIGNING...\n')

            for cmd in self.commands:

                if type( cmd ) == tuple:
                    for f in cmd:
                        self.repair_lib_file( f )
                else:

                    if host:
                        tc = TCoffee( cmd, host, verbose=self.verbose,
                                      debug=self.debug, log=self.log)
                    else:
                        tc = TCoffee( cmd, verbose=self.verbose,
                                      debug=self.debug, log=self.log)

                    ## run t_coffee
                    output, error, returncod = tc.run()


        except EnvironmentError, why:
            self.logWrite("ERROR: Can't run t_coffee: %s Error: %s"\
                          %( why, error ) )
            raise AlignerError( "Can't run t_coffee: %s Error: %"\
                                %( why, error ) )


#############
##  TESTING        
#############
import Biskit.test as BT    

class Test(BT.BiskitTest):
    """
    Base Test class (does not actually run any test)
    """
    def prepare(self):
        import tempfile
        import shutil, stat
        from Biskit import LogFile, StdLog

        ## collect the input files needed

        self.outfolder = tempfile.mkdtemp( '_test_Aligner' )

        os.mkdir( self.outfolder +'/templates' )
        os.mkdir( self.outfolder +'/sequences/' )

        shutil.copytree( T.testRoot() + '/Mod/project/templates/t_coffee',
                         self.outfolder + '/templates/t_coffee' )

        os.chmod( self.outfolder + '/templates/t_coffee', stat.S_IWRITE )

        shutil.copy( T.testRoot() + '/Mod/project/templates/templates.fasta',
                     self.outfolder + '/templates' )

        shutil.copy( T.testRoot() + '/Mod/project/sequences/nr.fasta',
                     self.outfolder + '/sequences/' )

        shutil.copy( T.testRoot() + '/Mod/project/target.fasta',
                     self.outfolder  )

        ## when copying from a write-protected biskit, the temp files are also
        ## becomming write-protected
        try:
            os.system( 'chmod -R +xrw %s' % self.outfolder )
        except Exception, error:
            T.lastError( str( error ) )

        if not self.local:
            self.a_log = LogFile( self.outfolder + '/Aligner.log' )
        else:
            self.a_log = StdLog()


    def t_Aligner(self, run=True ):

        self.a = Aligner( outFolder=self.outfolder, verbose=self.local,
                          log=self.a_log)

        self.a.align_for_modeller_inp()

        if run:
            self.a.go()
            if self.local:
                self.log.add(
                    'The alignment result is in %s/t_coffee' % self.outfolder)

    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )


class TestDry( Test ):
    """Test case dry run -- only initialize without running t-coffee"""

    TAGS = [BT.EXE]

    def test_AlignerDryRun(self):
        """Mod.Aligner dry run test"""
        self.t_Aligner(run=False)


class TestLong( Test ):
    """Complete test case for Aligner"""

    TAGS = [BT.EXE, BT.LONG]

    def test_Aligner(self):
        """Mod.Aligner test"""
        self.t_Aligner(run=True)

if __name__ == '__main__':

    BT.localTest(debug=0)
