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
##
## last $Author$
## last $Date$
## $Revision$

"""
Distribute calculation of contact matrices for many complexes
over many nodes.
"""

from Biskit.PVM import TrackingJobMaster
from Biskit.ReduceCoordinates import ReduceCoordinates
from Biskit.PVM.hosts import cpus_all, nice_dic
import Biskit.tools as t
import Biskit.mathUtils as MU
import Biskit.settings as settings
from Biskit.PDBDope import PDBDope
from Biskit.PDBModel import PDBProfiles
from Biskit.Errors import BiskitError
from Biskit.LogFile import StdLog
from Biskit import EHandler

from Complex import Complex
from ComplexList import ComplexList
from ComplexEvolvingList import ComplexEvolvingList

import numpy as N
import tempfile
import os.path

slave_path = t.projectRoot()+"/Biskit/Dock/ContactSlave.py"

class ContactMaster(TrackingJobMaster):
    """
    Distribute calculations done on all complexes of a ComplexList.
    """

    def __init__(self, complexLst, chunks=5, hosts=cpus_all, refComplex=None,
                 updateOnly=0, niceness = nice_dic, force = [],
                 outFile = 'complexes_cont.cl', com_version=-1,
                 show_output = 0, add_hosts=0, verbose=1, log=StdLog() ):
        """
        @param complexLst: input list
        @type  complexLst: ComplexList
        @param chunks: number of complexes to hand over to each slave
                       at each call
        @type  chunks: int
        @param hosts: list of hostnames
        @type  hosts: [str]
        @param niceness: dictionary mapping hosts to nice values::
                         {'hostname':nice_value, ..., 'default': nice_value}
                          nice_value: int or str
        @type  niceness: dict
        @param outFile: file name for output list of complexes with
                        calculated contacts (default: 'complexes_cont.cl')
        @type  outFile: str
        @param com_version: contact the given version of a ComplexEvolving,
                            only valid if input is ComplexEvolvingList
                            (default: -1)
        @type  com_version: int
        @param show_output: show x-window for each slave or not (default: 0)
        @type  show_output: 1|0
        @param add_hosts: add hosts to PVM (can take some time) (default: 0)
        @type  add_hosts: 1|0
        @param force: force update of given info keys
        @type  force: [str]
        @param verbose: print progress infos (default: 1)
        @type  verbose: 0|1

        @raise BiskitError: if attempting to extract version from list
                            that is not of type ComplexEvolvingList.
        """
        self.outFile = outFile

        self.verbose = verbose
        self.log = log

        ## set name for error log file
        self.ferror = os.path.dirname(
            t.absfile(self.outFile)) + '/contacting_errors.log'

        ## reference for fraction of native contacts
        self.c_ref_res = None
        self.c_ref_atom_4_5 = self.c_ref_atom_7_5 = self.c_ref_atom_10 = None

        ## masks to apply before comparing to reference atom matrix (casting)
        self.mask_rec, self.mask_lig = None, None

        ## for interface rms calculation
        self.ref_interface_model = None
        self.mask_rec_interface = self.mask_lig_interface = None

        ## update None dictonary entries
        self.updateOnly = updateOnly

        ## force update of given info keys
        self.force = force

        ## extract given version of each Complex from ComplexEvolvingList
        self.complexLst_original = None
        self.com_version = None

        if isinstance( complexLst, ComplexEvolvingList ):
            complexLst = self.__extractVersion( complexLst, com_version)
        else:
            if com_version != -1: raise BiskitError(
                'Input list must be ComplexEvolvingList to extract version.')

        ## make sure models have a surfaceMask
        self.__addSurfaceMasks( complexLst )

        ## store class type so that same class can be returned
        self.list_class = complexLst.__class__

        complexDic, self.remainLst = self.__complexDic( complexLst )

        if verbose:
            self.log.write( '%i complexes total in list\n'%len(complexLst) )
            self.log.write( '\tof which %i need updating,\n'%len(complexDic) )
            self.log.write( '\tand %i are complete\n'%len(self.remainLst) )

        self.__check_models( complexLst )

        ## coarse models for reduced atom contacts
        self.reduced_recs, self.reduced_ligs = {}, {}
        self.reduced_refCom = None
        self.c_ref_ratom_10 = None

        if not self.force or 'fnarc_10' in self.force or \
           'c_ratom_10' in self.force:

            self.__create_reduced_models( complexLst, refCom=refComplex )
            if refComplex:
                self.__reduced_refContacts()

        ## reference complex data
        self.c_ref_res_4_5 = self.c_ref_atom_4_5 = None
        self.c_ref_atom_10 = self.mask_rec = self.mask_lig = None
        self.mask_interface = self.mask_interface_bb = None
        self.ref_interface_model = self.ref_interface_model_bb = None

        ## get reference residue and atom contact matrices
        if refComplex:
            self.__refMasks( refComplex, complexLst[0] )

        self.finished = 0

        TrackingJobMaster.__init__( self, complexDic, chunks,
                                    hosts, niceness, slave_path, show_output=show_output,
                                    add_hosts=add_hosts, verbose=verbose )

        if verbose: print "JobMaster initialized."


    def __check_models( self, cl ):
        """
        Perform some checks for models assocciated with ComplexList.
        """
        ## stray models are not in the ComplexLists.models because the don't
        ## have a identical pickled source
        if cl.strayModels() != ({},{}):
            print """Warning: some receptor or ligand models have in-memory
            coordinates or profiles. This will slow down the PVM traffic
            and might have other bad consequences...
            Make sure that all models referenced in the ComplexList are
            PDBModel.isChangedFromDisc() == (0,0)!
            If that's the case but you still get this message, the model
            registry of the list may be messed up. Re-create the ComplexList
            with:
            cl = ComplexList( cl.toList() )
            """

        ## For some reason profiles are sometimes marked as "changed" and
        ## their content saved in the model even though the identical profile
        ## exists in the source model. This can blow up pvm traffic.
        ## BUT:!!
        ## If the profiles are used by the ContactSlave you should do the
        ## opposite, so that the slave doesn't need to unpickle the source
        ## model
##         self.log.write("synchronizing models with source ...")
##         for m in cl.recModels() + cl.ligModels():
##             m.update()
##             m.slim()
##         self.log.write("done\n")


    def __extractVersion( self, cel, com_version=-1 ):
        """
        Get a ComplexList by extracting a specified version (generation) of
        a ComplexEvolvingList.

        @param cel: ComplexEvolvingList
        @type  cel: ComplexEvolvingList
        @param com_version: version of ComplexEvolvingList to get
                            (default: -1, last version)
        @type  com_version: int

        @return: ComplexList 
        @rtype: ComplexList
        """
        self.complexLst_original = cel
        self.com_version = com_version

        ## mark complexes to detect mix-ups later on
        for i in range( len( cel ) ):
            cel[i][ com_version ]['$temp$'] = i

        return cel.toComplexList( version=com_version )


    def __addSurfaceMasks( self, cl ):
        """
        Add surface area to rec and lig model if not already there.

        @param cl: ComplexList
        @type  cl: ComplexList
        """
        for m in cl.recModels():
            if not m.profile('surfMask', 0) != 0:
                d = PDBDope( m )
                d.addSurfaceMask()

        for m in cl.ligModels():
            if not m.profile('surfMask', 0) != 0:
                d = PDBDope( m )
                d.addSurfaceMask()


    def __addSurface( self, m ):
        """
        Add surface area to model.

        @param m: model
        @type  m: PDBModel       
        """
        if m.profile('relAS', 0) is not 0:
            return 

        d = PDBDope( m )
        d.addSurfaceRacer()


    def __refMasks( self, refCom, normalCom ):
        """
        Create residue and atom contact masks for reference complex.

        @param refCom: reference complex
        @type  refCom: Complex
        @param normalCom: complex with different atom/residue content than
                          the reference complex to be used for determining
                          the sections that are identical in both complexes
        @type  normalCom: Complex
        """
        if self.verbose: self.log.write( 'Analyzing reference complex ... ')

        NC = normalCom; RC = refCom

        ## indices to apply for casting, cast reference complex, del. Hydrogens
        i_rec_ref, i_lig_ref, i_rec, i_lig = RC.compareAtoms( NC )
        h_rec = N.nonzero( RC.rec_model.maskH() )[0]
        h_lig = N.nonzero( RC.lig_model.maskH() )[0]
        i_rec_ref = [ i for i in i_rec_ref if i not in h_rec ]
        i_lig_ref = [ i for i in i_lig_ref if i not in h_lig ]

        RC = RC.take( i_rec_ref, i_lig_ref )

        ## convert casting indices for normalCom to mask
        m_rec = N.zeros( len( NC.rec_model ), N.int )
        N.put( m_rec, i_rec, 1 )
        m_lig = N.zeros( len( NC.lig_model ), N.int )
        N.put( m_lig, i_lig, 1 )

        self.mask_rec = m_rec * NC.rec_model.maskHeavy()
        self.mask_lig = m_lig * NC.lig_model.maskHeavy()

        ## reference residue contacts
        cont_4_5 = RC.resContacts( cutoff=4.5, refComplex = NC, cache_pw=1 )
        self.c_ref_res_4_5 = MU.packBinaryMatrix( cont_4_5 )

        cont_10 = RC.resContacts( cutoff=10., refComplex= NC, cache_pw=1 )

        ## reference atom contacts
        m = RC.atomContacts( 4.5, cache=1, map_back=0)
        self.c_ref_atom_4_5 = MU.packBinaryMatrix( m )

        m = RC.atomContacts( 10., cache=1, map_back=0)
        self.c_ref_atom_10  = MU.packBinaryMatrix( m )

        ## reference structure of contacting residues all atoms, cutoff 4.5
        x, y = self.__ref_interface( RC, NC, cont_4_5,
                                     self.mask_rec, self.mask_lig )
        self.mask_interface      = MU.packBinaryMatrix( x )
        self.ref_interface_model = y

        ## reference structure of contacting residues backbone only, cutoff 10
        x, y = self.__ref_interface( RC, NC, cont_10,
                                     self.mask_rec, self.mask_lig, bb=1 )
        self.mask_interface_bb      = MU.packBinaryMatrix( x )
        self.ref_interface_model_bb = y

        ## 
        self.mask_rec = MU.packBinaryMatrix( self.mask_rec )
        self.mask_lig = MU.packBinaryMatrix( self.mask_lig )

        if self.verbose: self.log.write('done\n')


    def __ref_interface( self, RC, NC, res_contacts, mask_rec, mask_lig, bb=0):
        """
        @param RC: reference complex, casted to normal complex
        @type  RC: Complex
        @param NC: normal complex, not casted to RC
        @type  NC: Complex
        @param res_contacts: residue contact matrix for defining
                             interface residues
        @type  res_contacts: matrix
        @param mask_rec: atom mask to cast NC.rec() to RC.rec()
        @type  mask_rec: [1|0]
        @param mask_lig: atom mask to cast NC.lig() to RC.lig()
        @type  mask_lig: [1|0]
        @param bb: apply a backbone mask (default: 0)
        @type  bb: 1|0

        @return: atom mask for getting interface from NC.model(),
                 reference interface
        @rtype: [int], PDBModel
        """
        try:
            ## calculate atom mask that extracts all interface residues from NC
            if_rec = NC.rec_model.res2atomMask( N.sum( res_contacts, 1) )
            if_rec = if_rec * NC.rec_model.maskHeavy() * mask_rec
            if bb:
                if_rec = if_rec * NC.rec_model.maskBB()

            if_lig = NC.lig_model.res2atomMask( N.sum( res_contacts, 0) )
            if_lig = if_lig * NC.lig_model.maskHeavy() * mask_lig
            if bb:
                if_lig = if_lig * NC.lig_model.maskBB()

            mask_interface = N.concatenate( (if_rec, if_lig) )

            ## extract all interface residues of ref complex RC into a PDBModel
            if_rec = N.compress( mask_rec, if_rec ) ## reduce to mask for RC
            if_lig = N.compress( mask_lig, if_lig )
            mask_interface_ref = N.concatenate( (if_rec, if_lig) )

            ref_interface_model = RC.model().compress(mask_interface_ref)
            ref_interface_model.residues = PDBProfiles()
            ref_interface_model.atoms = PDBProfiles()
            ref_interface_model.xyz = ref_interface_model.xyz.tolist()

            return mask_interface, ref_interface_model

        except ValueError, why:
            EHandler.error("ValueError (rescontacts: " + str(res_contacts)+")")


    def __reduced_refContacts( self ):
        """
        calculate contact mask of atom-reduced reference complex.
        """
        if self.verbose: self.log.write("get reduced reference contacts...")

        c = self.reduced_refCom.atomContacts( 10, cache=1, map_back=0)
        self.c_ref_ratom_10 = MU.packBinaryMatrix( c )

        if self.verbose: self.log.write(' done\n')


    def __toBeUpdated( self, com ):
        """
        Check if complex is to be updated.

        @param com: Complex
        @type  com: Complex

        @return: 1, if complex needs to be contacted
        @rtype: 1|0
        """
        if not self.updateOnly:
            return 1

        for ikey in com.keys():
            if com[ikey] == None:

                return 1

        return 0


    def __complexDic( self, cl ):
        """
        Collect info about complexes in ComplexList that needs to be updated.

        @param cl: ComplexList
        @type  cl: ComplexList

        @return: dictionary mapping solution number to Complex,
                 remaining complexes as a ComplexList
        @rtype: {int:Complex}, ComplexList
        """
        update = {}
        remain = ComplexList()

        if self.verbose: self.log.write('setting up task list ')
        for i in range( len( cl ) ):

            c = cl[i]

            if self.__toBeUpdated( c ):
                update[i] = c
            else:
                remain += [c]

        if self.verbose: self.log.write(' done\n')

        return update, remain


    def __create_reduced_models( self, cl, refCom=None ):
        """
        Create rec and lig models with pooled atoms (ca. every 3 atoms) but
        keep only those united atoms where the singular atoms have an average
        rel surf acc area > 30%. Do the same to the reference complex rec and
        lig. The result are PDBModels with actually less atoms than residues
        but whose contact matrix (and fnac) still resembles the one calculated
        from all atoms.

        @param cl: ComplexList
        @type  cl: ComplexList
        @param refCom: reference complex
        @type  refCom: Complex        
        """
        if self.verbose:
            self.log.write('preparing/saving coarse models of receptors&ligands...')

        rec_models = cl.recModels()
        lig_models = cl.ligModels()

        ## cast with reference complex, if available
        if refCom:
            ref_rec = refCom.rec_model
            ref_lig = refCom.lig_model

            i_rec, i_ref = rec_models[0].compareAtoms( ref_rec )
            ref_rec = ref_rec.take( i_ref )

            i_lig, i_ref = lig_models[0].compareAtoms( ref_lig )
            ref_lig = ref_lig.take( i_ref )

            rec_models = [ m.take( i_rec ) for m in rec_models ]
            lig_models = [ m.take( i_lig ) for m in lig_models ]

        ## add surface profile so that reducer will create averaged profile
        self.__addSurface( rec_models[0] )
        self.__addSurface( lig_models[0] )

        ## reduce all ligands and receptors
        r_reducer = ReduceCoordinates( rec_models[0] )
        l_reducer = ReduceCoordinates( lig_models[0] )

        r,l = {}, {}
        for m in rec_models:
            r[ m.source ] = r_reducer.reduceToModel( m.getXyz() )

        for m in lig_models:
            l[ m.source ] = l_reducer.reduceToModel( m.getXyz() )

        ## reduce refComplex
        if refCom:
            ref_rec = r_reducer.reduceToModel( ref_rec.getXyz() )
            ref_lig = l_reducer.reduceToModel( ref_lig.getXyz() )
            self.reduced_refCom = Complex( ref_rec, ref_lig,
                                           refCom.ligandMatrix )

        ## keep only (average) surface atoms
        ## play around with cutoff: lower- more atoms, better correl. to fnac
        i_r_surf = N.nonzero( r[ rec_models[0].source ].\
                              profile2mask( 'relAS', 30 ) )[0]
        i_l_surf = N.nonzero( l[ lig_models[0].source ].\
                              profile2mask( 'relAS', 30 ) )[0]

        for m in r.values(): m.keep( i_r_surf )
        for m in l.values(): m.keep( i_l_surf )

        if refCom:
            ref_rec.keep( i_r_surf )
            ref_lig.keep( i_l_surf )

        ## save changed models
        self.reduced_recs, self.reduced_ligs = r,l

        for m in l.values() + r.values():
            f = tempfile.mktemp( '_reduced.model', dir=settings.tempDirShared )
            m.saveAs(f)

        if self.verbose: self.log.write(' done\n')


    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.

        @param slave_tid: slave task id
        @type  slave_tid: int

        @return: dictionaty with slave parameters
        @rtype: dict        
        """
        r = {'ferror': self.ferror, 'force' : self.force,
             'c_ref_res_4_5' : self.c_ref_res_4_5,
             'c_ref_atom_4_5': self.c_ref_atom_4_5,
             'c_ref_atom_10' : self.c_ref_atom_10,
             'ref_interface_model'   : self.ref_interface_model,
             'ref_interface_model_bb': self.ref_interface_model_bb,
             'mask_rec'       : self.mask_rec,
             'mask_lig'       : self.mask_lig,
             'mask_interface'   : self.mask_interface,
             'mask_interface_bb': self.mask_interface_bb,
             'reduced_recs' : self.reduced_recs,
             'reduced_ligs' : self.reduced_ligs,
             'c_ref_ratom_10'  : self.c_ref_ratom_10,
             'log':self.log.fname,
             'verbose':self.verbose
             }

        return r


    def isFinished( self ):
        return self.finished


    def cleanup( self ):
        """
        Remove temporary files.
        Overrides TrackingJobMaster method
        """
        if self.verbose: print "Deleting temporary coarse rec/lig models..."

        for m in self.reduced_recs.values() + self.reduced_ligs.values():
            if not t.tryRemove( str( m.source ) ):
                print "Cannot remove ", str( m.source )


    def getResult( self, **arg ):
        """
        Return result ComplexList, if it is available.

        @return: resulting ComplexList
        @rtype: ComplexList
        """
        return self.complexLst


    def done(self ):
        """
        Collect the last complexes and write result ComplexList to file.

        @raise BiskitError: if Complex version mixup
        """
        if self.verbose: print "Assembling new ComplexList...",
        self.complexLst = self.list_class( self.result.values() )
        self.result.clear()

        self.complexLst += self.remainLst

        ## update complexes in ComplexEvolvingList 
        if self.complexLst_original is not None:

            if self.verbose:
                print "Copying values into version %i ..." % self.com_version,

            for cev, c in zip( self.complexLst_original, self.complexLst ):

                check_1 = c['$temp$']
                check_2 = cev[ self.com_version ]['$temp$']
                del c['$temp$']
                del cev[self.com_version]['$temp$']
                if not check_1 == check_2:
                    raise BiskitError('Complex version mixup: %i != %i' \
                                      % (check_1, check_2) )

                cev[ self.com_version ].info.update( c.info )

            self.complexLst = self.complexLst_original

        if self.verbose: print "\nSaving result to %s..." % self.outFile
        t.dump( self.complexLst, self.outFile )

        self.finished = 1


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case

    Requires PVM installed and running and at least one (non-master) node.
    """

    TAGS = [ BT.PVM ]

    def prepare(self):
        self.cl_out = tempfile.mktemp('_test.cl')
        self.master = None


    def cleanUp(self):
        if self.master is not None:
            t.tryRemove( self.master.outFile )
            t.tryRemove( self.master.ferror )


    def test_ContactMaster(self):
        """Dock.ContactMaster test"""
        niceness = {'default': 0}
        self.hosts = cpus_all[:4]

        lst = t.load( t.testRoot() + "/dock/hex/complexes.cl")
        lst = lst[:9]

        refcom = t.load( t.testRoot() + "/com/ref.complex")

        self.master = ContactMaster( lst, chunks = 3, hosts = self.hosts,
                                     niceness = niceness,
                                     show_output = self.local,
                                     verbose = self.local,
                                     refComplex = refcom,
                                     outFile = self.cl_out )

        self.assert_( len(self.hosts) > 0,
                      'master needs at least one pvm node for calculations' )

        if len(self.hosts) > 0:
        ## wait for calculation to finish, then load contacted list
            self.cl_cont = self.master.calculateResult()

            if self.local:
                print 'Any error are written to: %s'%master.ferror
                ## plot atom and residue contacts vs. rmsd
                self.p = self.cl_cont.plot( 'rms', 'fnac_10','fnarc_10' )
                self.p.show()

            self.assertAlmostEqual(N.sum(self.cl_cont.valuesOf('fnac_10')),
                                   0.50811038550663579, 5 )


if __name__ == '__main__':

    BT.localTest()

