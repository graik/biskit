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

## Contributions: Olivier PERIN
##
## last $Author$
## $Date$
## $Revision$
"""
Modeling benchmark
"""

import os
from Biskit.PDBModel import PDBModel
from Biskit.ModelList import ModelList
from Biskit.IcmCad import IcmCad as CAD

import Biskit.tools as T
import numpy as N
from Biskit import StdLog, EHandler

from Biskit.Mod import Modeller

class Benchmark:
    """
    This class regroups the different methods to analyse
    the output files of interest from modeller
    """
    # standard directory for input
    F_INPUT_FOLDER = Modeller.F_RESULT_FOLDER

    # standard directory for output
    F_RESULT_FOLDER = '/benchmark'

    # Standard Input Files
    # PDB of the Known Structure
    F_INPUT_REFERENCE= '/reference.pdb'

    # PDBModels list Pickle
    F_PDBModels = F_INPUT_FOLDER + '/PDBModels.list'

    # Standard Output files for RMSD
    F_RMSD_AA = F_RESULT_FOLDER +'/rmsd_aa.out'
    F_RMSD_CA = F_RESULT_FOLDER +'/rmsd_ca.out'
    F_RMSD_RES = '/rmsd_res'

    # standard Ouput file for PDBModels
    F_PDBModels_OUT = F_RESULT_FOLDER +'/PDBModels.list'

    # standard Ouput file for final reference
    F_FINAL_REFERENCE = F_RESULT_FOLDER +'/reference_final.pdb'


    def __init__(self, outFolder='.', verbose=1):
        """
        @param outFolder: folder for output (default: .)
        @type  outFolder: str
        @param verbose: write intermediary files (default: 1)
        @type  verbose: 1|0
        """
        self.outFolder = T.absfile( outFolder )

        self.verbose = verbose

        self.prepareFolders()


    def prepareFolders( self ):
        """
        Create folders needed by this class.
        """
        if not os.path.exists(self.outFolder + self.F_RESULT_FOLDER):
            os.mkdir(self.outFolder + self.F_RESULT_FOLDER)


    def output_fittedStructures(self, pdb, reference, index, atom_mask,
                                output_folder = None):
        """
        Takes a model and a reference structure, performs both a
        normal fillting to the reference and an itterative fitting.
        Then returns the two fitted models.

        @param pdb: model
        @type  pdb: PDBModel
        @param reference: reference model
        @type  reference: PDBModel
        @param index: index number
        @type  index: int
        @param atom_mask: atom mask for discharding certain atoms
        @type  atom_mask: [int]
        @param output_folder: output folder
                       (default: None S{->} outFolder/L{F_RESULT_FOLDER}
        @type  output_folder: str

        @return: pdb_if, pdb -  model, fitted iteratively (n=10),
                                model, without iterative fitting
        @rtype: PDBModel, PDBModel
        """
        output_folder = output_folder or self.outFolder + self.F_RESULT_FOLDER

        tmp_model = pdb.compress( atom_mask )

        ## iterative fit, max 10 itterations
        pdb_if = pdb.transform( *tmp_model.transformation\
                                ( reference, n_it=10,
                                  profname="rms_outliers") )
        ## normal fit
        pdb = pdb.transform( *tmp_model.transformation(reference, n_it=0) )

        ## info about discarded residues in iterative fit
        pdb_if.atoms.set( "rms_outliers",
                          tmp_model.profile("rms_outliers"),
                          mask=atom_mask )

        ## wirte iteratively fitted structure to disc
        pdb_if.writePdb( '%s/Fitted_%02i.pdb'%( output_folder, index ) )

        return pdb_if, pdb


    def calc_rmsd(self, fitted_model_if, fitted_model_wo_if, reference, model):
        """
        Takes the two fitted structures (with and without iterative fitting),
        the known structure (reference), and the associated model inside the
        pdb_list. Calculates the different RMSD and set the profiles

        @param fitted_model_if: itteratively fitted model
        @type  fitted_model_if: PDBModel
        @param fitted_model_wo_if: normaly fitted model
        @type  fitted_model_wo_if: PDBModel
        @param reference: reference model
        @type  reference: PDBModel
        @param model: model
        @type  model: PDBModel
        """
        ## first calculate rmsd for heavy atoms and CA without
        ## removing any residues from the model
        mask_CA = fitted_model_wo_if.maskCA()

        rmsd_aa = fitted_model_wo_if.rms( reference, fit=0 )
        rmsd_ca = fitted_model_wo_if.rms( reference, mask=mask_CA, fit=1 )

        model.info["rmsd2ref_aa_wo_if"] = rmsd_aa
        model.info["rmsd2ref_ca_wo_if"] = rmsd_ca

        outliers_mask = N.logical_not(fitted_model_if.profile("rms_outliers"))

        ## Now remove the residues that were outliers in the iterative fit
        ## and calculate the rmsd again
        fitted_model_if = fitted_model_if.compress( outliers_mask )
        reference = reference.compress( outliers_mask )

        mask_CA = fitted_model_if.maskCA()

        rmsd_aa_if = fitted_model_if.rms( reference, fit=0 )
        rmsd_ca_if = fitted_model_if.rms( reference, mask=mask_CA, fit=1 )

        model.info["rmsd2ref_aa_if"] = rmsd_aa_if
        model.info["rmsd2ref_ca_if"] = rmsd_ca_if
        model.info["rmsd2ref_aa_outliers"] = 1.*(len(outliers_mask) \
                                                 - N.sum(outliers_mask)) / len(outliers_mask)
        model.info["rmsd2ref_ca_outliers"] = 1.*(N.sum(mask_CA) \
                                                 - N.sum(N.compress(mask_CA, outliers_mask))) \
             / N.sum(mask_CA)


    def output_rmsd_aa(self, pdb_list, output_file = None):
        """
        Write a file rmsd_aa.dat that contains
        all the global heavy atom rmsd values.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param output_file: output file
                            (default: None S{->} outFolder/L{F_RMSD_AA})
        @type  output_file: str
        """
        output_file = output_file or self.outFolder + self.F_RMSD_AA
        rmsd_aa_out = open(output_file, 'w')

        rmsd_aa_out.write("# Heavy atom RMSD.\n")
        rmsd_aa_out.write("#  column 1. Normal fitting\n")
        rmsd_aa_out.write("#  column 2. Itterative fitting\n")
        rmsd_aa_out.write("#  column 3. Fraction of discarded residues in 2\n")
        for i in range(len(pdb_list)):
            rmsd_aa_out.write("model_%02i\t%f\t%f\t%f\n"\
                              %(i, pdb_list[i].info["rmsd2ref_aa_wo_if"],
                                pdb_list[i].info["rmsd2ref_aa_if"],
                                pdb_list[i].info["rmsd2ref_aa_outliers"]))
        rmsd_aa_out.close()


    def output_rmsd_ca(self, pdb_list, output_file = None):
        """
        Write a file rmsd_ca.dat that contains
        all the global c_alpha rmsd values.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param output_file: output file
                            (default: None S{->} outFolder/L{F_RMSD_CA})
        @type  output_file: str        
        """
        output_file = output_file or self.outFolder + self.F_RMSD_CA
        rmsd_ca_out = open('%s'%output_file, 'w')

        rmsd_ca_out.write("# Carbon alpha RMSD.\n")
        rmsd_ca_out.write("#  column 1. Normal fitting\n")
        rmsd_ca_out.write("#  column 2. Itterative fitting\n")
        rmsd_ca_out.write("#  column 3. Fraction of discarded residues in 2\n")
        for i in range(len(pdb_list)):
            rmsd_ca_out.write("model_%02i\t%7.3f\t%7.3f\t%7.3f\n"\
                              %(i, pdb_list[i].info["rmsd2ref_ca_wo_if"],
                                pdb_list[i].info["rmsd2ref_ca_if"],
                                pdb_list[i].info["rmsd2ref_ca_outliers"] ))
        rmsd_ca_out.close()


    def rmsd_res(self, coord1, coord2):
        """
        Calculate the rsmd on residue level for c-alpha between a
        model and its reference.

        @param coord1: first set of coordinates
        @type  coord1: array
        @param coord2: second set of coordinates
        @type  coord2: array

        @return: rmsd_res: rmsd per c-alpha
        @rtype: [float]
        """
        rmsd_res = []

        for i in range( len(coord1) ):
            rmsd = N.sqrt( (N.power(coord1[i][0]-coord2[i][0],2) +  \
                            N.power(coord1[i][1]-coord2[i][1],2 )+ \
                            N.power(coord1[i][2]-coord2[i][2],2 )))
            rmsd_res.append(rmsd)

        return rmsd_res


    def output_rmsd_res(self, pdb_list, output_folder = None):
        """
        Write a file that contains the rmsd profile on a residue
        level for c_alpha.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param output_folder: output folder (default: None S{->} outFolder/
                              L{F_RESULT_FOLDER}/L{F_RMSD_RES})
        @type  output_folder: str      
        """
        for m, t in zip(pdb_list, range(len(pdb_list))):
            output_file = output_folder or self.outFolder + \
                        self.F_RESULT_FOLDER + self.F_RMSD_RES + '_%02i'%t

            file = open(output_file, 'w')
            rmsd_res_list = m.compress(m.maskCA()).atoms["rmsd2ref_if"]

            file.write("# Carbon alpha rmsd profile.\n")
            for i in range(len(rmsd_res_list)):
                file.write("%4i\t%5.3f\n"%(i+1,rmsd_res_list[i]))

            file.close()


    def cad(self, reference, model):
        """
        Calculates the CAD Contact Area Difference between
        the model and its reference structure and set a profile
        for the model

        @param reference: reference model
        @type  reference: PDBModel
        @param model: model
        @type  model: PDBModel
        """
        reference = reference.compress(reference.maskProtein())
        model = model.compress(model.maskProtein())

        model_list = []
        model_list.append(model)

        x = CAD(reference, model_list, debug=0, verbose=1)
        r = x.run()

        model.info["CAD"] = r


    def write_PDBModels(self, pdb_list, output_file = None):
        """
        Pickles the list of PDBModels to disc.

        @param pdb_list: list of models
        @type  pdb_list: ModelList
        @param output_file: output file
                       (default: None S{->} outFolder/L{F_PDBModels_OUT})
        @type  output_file: str         
        """
        output_file = output_file or self.outFolder + self.F_PDBModels_OUT
        T.dump(pdb_list, '%s'%(output_file))


    def go(self, model_list = None, reference = None):
        """
        Run benchmarking.

        @param model_list: list of models
                           (default: None S{->} outFolder/L{F_PDBModels})
        @type  model_list: ModelList
        @param reference: reference model
                        (default: None S{->} outFolder/L{F_INPUT_REFERENCE})
        @type  reference: PDBModel
        """
        model_list = model_list or self.outFolder + self.F_PDBModels
        reference = reference or self.outFolder + self.F_INPUT_REFERENCE

        pdb_list = T.load('%s'%model_list)
        reference = PDBModel(reference)

        # check with python 2.4
        iref, imodel = reference.compareAtoms(pdb_list[0])

        mask_casting = N.zeros(len(pdb_list[0]))
        N.put(mask_casting, imodel, 1)

        reference = reference.take(iref)
        #reference_mask_CA = reference_rmsd.maskCA()

        atom_mask = N.zeros(len(pdb_list[0]))
        N.put(atom_mask,imodel,1)

        rmask = pdb_list[0].profile2mask("n_templates", 1,1000)
        amask = pdb_list[0].res2atomMask(rmask)

        mask_final_ref = N.compress(mask_casting, amask)
        mask_final = mask_casting * amask

        reference = reference.compress(mask_final_ref)

        for i in range(len(pdb_list)):

            #self.cad(reference, pdb_list[i])

            pdb_list[i], pdb_wo_if = self.output_fittedStructures(\
                pdb_list[i], reference, i, mask_final)

            fitted_model_if = pdb_list[i].compress(mask_final)
            fitted_model_wo_if = pdb_wo_if.compress(mask_final)

            coord1 = reference.getXyz()
            coord2 = fitted_model_if.getXyz()

            aprofile = self.rmsd_res(coord1,coord2)

            self.calc_rmsd(fitted_model_if, fitted_model_wo_if,
                           reference, pdb_list[i])

            pdb_list[i].atoms.set('rmsd2ref_if', aprofile,
                                  mask=mask_final, default = -1,
                                  comment="rmsd to known reference structure")

        self.output_rmsd_aa(pdb_list)
        self.output_rmsd_ca(pdb_list)
        self.output_rmsd_res(pdb_list)

        self.write_PDBModels(pdb_list)



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
        from Biskit import PDBModel

        ## collect the input files needed
        self.outfolder = tempfile.mkdtemp( '_test_Benchmark' )
        os.mkdir( self.outfolder +'/modeller' )

        dir = '/Mod/project/validation/1DT7'
        shutil.copy( T.testRoot() + dir + '/modeller/PDBModels.list',
                     self.outfolder + '/modeller' )    
        shutil.copy( T.testRoot() + dir + '/reference.pdb',
                     self.outfolder )


    def cleanUp(self):
        T.tryRemove( self.outfolder, tree=1 )

    def test_Benchmark(self):
        """Mod.Benchmark test"""
        from Biskit import Pymoler

        self.b = Benchmark( self.outfolder )

        self.b.go()

        pdb = T.load( self.outfolder + "/modeller/PDBModels.list" )[0]

        reference = PDBModel(self.outfolder  + "/reference.pdb" )
        tmp_model = pdb.clone()

        reference = reference.compress( reference.maskCA() )
        pdb       = pdb.compress( pdb.maskCA() )
        tmp_model = tmp_model.compress(tmp_model.maskCA())

        tm = tmp_model.transformation( reference, n_it=0,
                                       profname="rms_outliers")
        pdb = pdb.transform( tm )

        if self.local:
            pm = Pymoler()
            pm.addPdb( pdb, "m" )
            pm.addPdb( reference, "r" )
            pm.colorAtoms( "m", tmp_model.profile("rms_outliers") )
            pm.add('set ribbon_trace,1')
            pm.add('show ribbon')
            pm.show()

            if self.DEBUG:
                self.log.add(
                    'The result from the benchmarking is in %s/benchmark'%\
                    self.outfolder)

            globals().update( locals() )


if __name__ == '__main__':

    BT.localTest()
