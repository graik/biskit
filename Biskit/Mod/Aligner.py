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
## last $Author$
## last $Date$

import re, os
import commands

import settings
import modUtils
from TemplateCleaner import TemplateCleaner as TC
from SequenceSearcher import SequenceSearcher as SS

import Biskit.tools as tools
from Biskit.Errors import *

class AlignerError( BiskitError ):
    pass

class Aligner:
    """
    Take template CA traces and template sequences plus non-template
    sequences.
    -> T-Coffee script or T-Coffee alignment
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

    def __init__( self, outFolder='.', log=None, verbose=1 ):
        """
        outFolder - str, base folder for t_coffee output folder
        log - LogFile - deviate all messages to file [None=STDOUT]
        """
        self.log = log
        self.outFolder = tools.absfile( outFolder )
        self.verbose = verbose

        ## recognize file types for adding t_coffee type code
        self.ex_inp_types = { 'P' : re.compile('.+\.[Pp][Dd][Bb]$|'+
                                               '.+\.[Aa][Ll][Pp][Hh][Aa]$'),
                              'L' : re.compile('.+\.[Ll][Ii][Bb]$') }

        self.commands = []  ## list of t_coffee commands to run

        fn = self.outFolder + self.F_RESULT_FOLDER
        if not os.path.exists( fn ):
            self.logWrite( 'Creating new output folder '+ fn)
            os.mkdir( fn)


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg


    def __add_type_code( self, inp ):
        """
        Try guessing a t_coffee type code for file name and add type code.
        inp - str, file name
        -> type_code + file name OR file name unchanged
        """
        for code, pattern in self.ex_inp_types.items():
            if pattern.match( inp ):
                return code + inp
        return inp


    def coffee_align_inp( self, input, method=None, out_lib=None,
                          output=None, outfile=None, **kw ):
        """
        Create a single t_coffee command.
        input  - [ str ], list of input files
        method - str, t_coffee method ('M' is added automatically) [None]
        output - [ str ], output format(s) [None]
        out_lib- str, output library file if needed [None]
        outfile- str, output file (base) if needed [None]
        ..     - additional param=value pairs
        -> str, t_coffee command
        """
        r = settings.t_coffee_bin

        input = [ self.__add_type_code( i ) for i in tools.toList(input) ]

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
            for o in tools.toList( output ):
                r += ' ' + o

        if outfile:
            r += ' -outfile ' + outfile

        for param, value in kw.items():
            r += ' -'+param + ' ' + str( value )

        return r


    def __default_input_files( self, pdbFiles=None, fasta_templates=None,
                                fasta_sequences=None, fasta_target=None ):
        """
        Fetch missing fasta and pdb files from standard locations.
        -> { .. } 
        """
        r = {}

        ## collect CA trace PDBs from default location
        if not pdbFiles:
            folder = self.outFolder + TC.F_COFFEE
            fs = os.listdir( folder )
            pdbFiles= [ folder + f for f in fs if f[-6:].upper()=='.ALPHA' ]

        r['pdbs'] = tools.toList( pdbFiles )
        r['templates'] = fasta_templates or self.outFolder + TC.F_FASTA
        r['sequences'] = fasta_sequences or self.outFolder + SS.F_FASTA_NR
        r['target']    = fasta_target or self.outFolder + SS.F_FASTA_TARGET

        return r

        
    def align_for_modeller_inp( self, pdbFiles=None, fasta_templates=None,
                                fasta_sequences=None, fasta_target=None,
                                f_fast_tree=None, f_sequence_tree=None ):
        """
        Prepare alignment commands for homology modelling.
        If verbose==1, the commands are mirrored to t_coffee.inp .
        pdbFiles  - [ str ], template PDBs [from templates/nr/chain_index.txt ]
        fasta_templates  - str, template sequences [templates/templates.fasta ]
        fasta_sequences  - str, more sequences           [ sequences/nr.fasta ]
        """
        ## Fetch default files where needed
        d = self.__default_input_files( pdbFiles, fasta_templates,
                                        fasta_sequences, fasta_target )

        pdbFiles = d['pdbs']
        f_templates = d['templates']
        f_sequences = d['sequences']
        f_target    = d['target']
        
        pdbFiles = [ tools.absfile( f ) for f in pdbFiles ]

        f_templates = tools.absfile( f_templates )
        f_sequences = tools.absfile( f_sequences )

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
        if len( pdbFiles ) == 1:
            print '\n WARNING! Only one template avaliable:' + str(pdbFiles)
            print '          Structural alignment (SAP) will not be performed'
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
                                       output=['clustalw','phylip','score_html',
                                               'pir_aln'],
##                                     run_name=tools.stripFileName(f_final_aln),
                                       outfile=f_final_aln,
                                       quiet=f_coffee_log+'_4')
                ]
            
        ## Normal alignment run with structural alignment (more than one template)
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
                ## TODO: skip this step if there is only one template sequence
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
                                       output=['clustalw','phylip','score_html',
                                               'pir_aln', 'score_ascii'],
 ##                                    run_name=tools.stripFileName(f_final_aln),
                                       outfile=f_final_aln,
                                       quiet=f_coffee_log+'_4')
                ]

        ## add to internal command 'queue'
        self.commands += r

        if self.verbose:
            f = open( self.outFolder + self.F_COFFEE_INP, 'a' )
            for cmd in r:
                f.write( str(cmd) + '\n\n' )
            f.close()
        

    def repair_lib_file(self, fname ):
        """
        Msap_pair alignments mess up sequence ids of certain
        (single-chain?) structures.
        Dirty hack..
        """
        ex_sapfix = re.compile('_[A-Z]{0,1}_[A-Z]')
        bak_fname = fname+'_original'

        os.rename( fname,  bak_fname)

        fout = open( fname, 'w' )
        for l in open( bak_fname ):
            fout.write( re.sub( ex_sapfix, '_', l ) )

        fout.close()


    def repair_target_fasta( self, fname ):
        """
        Change target ID to 'target' and remove '|' etc. that mess up
        things in the alignments.
        !! AlignerError
        """
        bak_fname = fname+'_original'
        os.rename( fname,  bak_fname)

        try:
            fout = open( fname, 'w' )
            fin = open( bak_fname )

            l0 = fin.readline()
            
            fout.write( '>target\n' )

            for l in fin:
                fout.write( l )

            fin.close()
            fout.close()

        except Exception, why:
            print tools.lastError()
            raise AlignerError("Cannot fix target sequence file: "+str(why))
        

    def go( self, host=None ):
        """
        Run the previously added commands, delete internal command list.
        host - str, host name for remote execution [None=local execution]
        """
        try:
            self.logWrite('\nALIGNING...\n')
            for cmd in self.commands:

                if type( cmd ) == tuple:
                    for f in cmd:
                        self.repair_lib_file( f )
                else:

                    if host:
                        cmd = "ssh %s '%s'" % (host, cmd)

                    self.logWrite( 'Running t_coffee .. ')
                    self.logWrite( cmd )

                    status, output = commands.getstatusoutput( cmd )

                    self.logWrite( '..done: ' + str( output ) + '\n' )

        except EnvironmentError, why:
            self.logWrite("ERROR: Can't run t_coffee: "+ str( why ) )
            raise AlignerError( "Can't run t_coffee: " + str( why ) )
        

##########
## TEST ##
##########
if __name__ == '__main__':

    a = Aligner( outFolder=tools.projectRoot()+'/test/Mod/project/' )

    a.align_for_modeller_inp()

    a.go()
