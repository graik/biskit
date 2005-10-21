##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
"""
Distribute calculation of contact matrices for many complexes
over many nodes.
"""

from Biskit.PVM.TrackingJobMaster import TrackingJobMaster
from Biskit.ReduceCoordinates import ReduceCoordinates
from Biskit.hosts import cpus_all, nice_dic
import Biskit.tools as t
import Biskit.mathUtils as MU
import Biskit.settings as settings
from Biskit.PDBDope import PDBDope
from Biskit.PDBModel import PDBProfiles
from Biskit.Errors import BiskitError

from Complex import Complex
from ComplexList import ComplexList
from ComplexEvolvingList import ComplexEvolvingList

import Numeric as N
import tempfile
import os.path


slave_path = t.projectRoot()+"/Biskit/Dock/ContactSlave.py"

class ContactMaster(TrackingJobMaster):
    """
    Distribute calculations done on all complexes of a ComplexList.
    """

    def __init__(self, complexLst, chunks=5, hosts=cpus_all, refComplex=None,
                 updateOnly=0, niceness = nice_dic, force = [],
                 outFile = 'complexes_cont.cl',
                 com_version=-1,
                 show_output = 0, add_hosts=0,
                 xplor_restraints=None, xplor_scale=None ):
        """
        cmplxLst       - ComplexList, input list
        chunks         - int, number of complexes to hand over to each slave
                         at once
        hosts          - [str], hostnames
        niceness       - {'hostname':nice_value, ..., 'default': nice_value}
                         nice_value: int or str
        cutoff         - float, distance (A) cutoff for atom-atom contacts
        outFile        - str, file name for list of complexes with contacts
        com_version    - int, contact the given version of a ComplexEvolving,
                         only valid if input is ComplexEvolvingList [-1]
        show_output    - 1||0, show x-window for each slave or not
        add_hosts      - 1||0, add hosts to PVM (can take some time) [0]
        force          - [str], force update of given info keys

        for xplorEnergy only:
          xplor_restraints - lst
          xplor_scale      - lst
        """
        self.outFile = outFile

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

        ## for xplorEnergy only
        self.xplor_restraints = xplor_restraints
        self.xplor_scale = xplor_scale
             
        ## extract given version of each Complex from ComplexEvolvingList
        self.complexLst_original = None
        self.com_version = None
        
        if isinstance( complexLst, ComplexEvolvingList ):
            complexLst = self.__extractVersion( complexLst, com_version)
        else:
            if com_version != -1: raise BiskitError(
                'Input list must be ComplexEvolvingList to extract version.')

##         ## make sure models have a surfaceMask
##         self.__addSurfaceMasks( complexLst )
 
        ## store class type so that same class can be returned
        self.list_class = complexLst.__class__

        complexDic, self.remainLst = self.__complexDic( complexLst )

        t.flushPrint( '%i complexes total in list\n'%len(complexLst) )
        t.flushPrint( '\tof which %i need updating,\n'%len(complexDic) )
        t.flushPrint( '\tand %i are complete\n'%len(self.remainLst) )
         
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
               add_hosts=add_hosts )

        print "JobMaster initialized."


    def __check_models( self, cl ):

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
##         t.flushPrint("synchronizing models with source ...")
##         for m in cl.recModels() + cl.ligModels():
##             m.update()
##             m.slim()
##         t.flushPrint("done\n")

                        
    def __extractVersion( self, cel, com_version=-1 ):
        """
        cel          - ComplexEvolvingList
        com_version  - int
        -> ComplexList
        """

        self.complexLst_original = cel
        self.com_version = com_version

        ## mark complexes to detect mix-ups later on
        for i in range( len( cel ) ):
            cel[i][ com_version ]['$temp$'] = i
        
        return cel.toComplexList( version=com_version )


    def __addSurfaceMasks( self, cl ):
        """
        add surface area to model if not already there
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
        add surface area to model
        """
        if m.profile('relASA', 0) != 0:
            return 
        
        d = PDBDope( m )
        d.addASA()


    def __refMasks( self, refCom, normalCom ):
        """
        create reference residue and atom contact masks.
        """
        t.flushPrint( 'Analyzing reference complex ... ')

        NC = normalCom; RC = refCom

        ## indices to apply for casting, cast reference complex, del. Hydrogens
        i_rec_ref, i_lig_ref, i_rec, i_lig = RC.compareAtoms( NC )
        h_rec = N.nonzero( RC.rec_model.maskH() )
        h_lig = N.nonzero( RC.lig_model.maskH() )
        i_rec_ref = [ i for i in i_rec_ref if i not in h_rec ]
        i_lig_ref = [ i for i in i_lig_ref if i not in h_lig ]
        
        RC = RC.take( i_rec_ref, i_lig_ref )

        ## convert casting indices for normalCom to mask
        m_rec = N.zeros( len( NC.rec_model ), 'i' )
        N.put( m_rec, i_rec, 1 )
        m_lig = N.zeros( len( NC.lig_model ), 'i')
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

        t.flushPrint('done\n')


    def __ref_interface( self, RC, NC, res_contacts, mask_rec, mask_lig, bb=0):
        """
        RC - Complex, reference complex, casted to normal complex
        NC - Complex, normal complex, not casted to RC
        res_contacts - residue contact matrix for defining interface residues
        mask_rec - [ int ], atom mask to cast NC.rec() to RC.rec()
        mask_lig - [ int ], atom mask to cast NC.lig() to RC.lig()
        -> ([ int ], PDBModel),
           atom mask for getting interface from NC.model(), reference interface
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
            ref_interface_model.rProfiles = PDBProfiles()
            ref_interface_model.aProfiles = PDBProfiles()
            ref_interface_model.xyz = ref_interface_model.xyz.tolist()

            return mask_interface, ref_interface_model

        except ValueError, why:
            EHandler.error("ValueError (rescontacts: " + str(res_contacts)+")")
        

    def __reduced_refContacts( self ):
        """
        calculate contact mask of atom-reduced reference complex.
        """
        t.flushPrint("get reduced reference contacts...")

        c = self.reduced_refCom.atomContacts( 10, cache=1, map_back=0)
        self.c_ref_ratom_10 = MU.packBinaryMatrix( c )

        t.flushPrint(' done\n')


    def __toBeUpdated( self, com ):
        """-> 1, if complex needs to be contacted"""

        if not self.updateOnly:
            return 1

        for ikey in com.keys():
            if com[ikey] == None:
                
                return 1
        
        return 0


    def __complexDic( self, cl ):
        """
        -> dict {soln:Complex,..}
        """
        update = {}
        remain = ComplexList()

        t.flushPrint('setting up task list ')
        for i in range( len( cl ) ):

            c = cl[i]

            if self.__toBeUpdated( c ):
                update[i] = c
            else:
                remain += [c]

        t.flushPrint(' done\n')

        return update, remain


    def __create_reduced_models( self, cl, refCom=None ):
        """
        Create rec and lig models with pooled atoms (ca. every 3 atoms) but
        keep only those united atoms where the singular atoms have an average
        rel surf acc area > 30%. Do the same to the reference complex rec and
        lig. The result are PDBModels with actually less atoms than residues
        but whose contact matrix (and fnac) still resembles the one calculated
        from all atoms. 
        """
        t.flushPrint('preparing/saving coarse models of receptors&ligands...')

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
                              profile2mask( 'relASA', 30 ) )
        i_l_surf = N.nonzero( l[ lig_models[0].source ].\
                              profile2mask( 'relASA', 30 ) )

        for m in r.values(): m.keep( i_r_surf )
        for m in l.values(): m.keep( i_l_surf )

        if refCom:
            ref_rec.keep( i_r_surf )
            ref_lig.keep( i_l_surf )

        ## save changed models
        self.reduced_recs, self.reduced_ligs = r,l

        oldTemp, tempfile.tempdir = tempfile.tempdir, settings.tempDirShared

        for m in l.values() + r.values():
            f = tempfile.mktemp( '_reduced.model' )
            m.saveAs(f)

        tempfile.tempdir = oldTemp

        t.flushPrint(' done\n')
    

    def getInitParameters(self, slave_tid):
        """
        hand over parameters to slave once.
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
             'xplor_restraints' : self.xplor_restraints,
             'xplor_scale'  : self.xplor_scale
             }
        
        return r

    def isFinished( self ):
        return self.finished


    def cleanup( self ):
        """Overrides TrackingJobMaster method"""

        print "Deleting temporary coarse rec/lig models..."
        
        for m in self.reduced_recs.values() + self.reduced_ligs.values():
            if not t.tryRemove( str( m.source ) ):
                print "Cannot remove ", str( m.source )


    def getResult( self, **arg ):
        """
        Return result ComplexList, if it is available.
        -> ComplexList
        """
        return self.complexLst


    def done(self ):
        """
        Collect the last complexes and write result ComplexList to file.
        """
        print "Assembling new ComplexList...",
        self.complexLst = self.list_class( self.result.values() )
        self.result.clear()

        self.complexLst += self.remainLst

        ## update complexes in ComplexEvolvingList 
        if self.complexLst_original is not None:

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

        print "\nSaving result to %s..." % self.outFile
        t.Dump( self.complexLst, self.outFile )

        self.finished = 1


## TESTING ###
##############

if __name__ == '__main__':
    niceness = {'default': 0}
    hosts = cpus_all[:4]

##    lst = t.Load( t.testRoot() + "/dock/hex01/complexes.cl")
    lst = t.Load( t.absfile('~/interfaces/c11/dock_xray/hex01/complexes.cl'))
    lst = lst[:15]

    refcom = t.Load( t.absfile('~/interfaces/c11/com_wet/ref.complex'))

    master = ContactMaster(lst, chunks=3, hosts=hosts, niceness=niceness,
                           show_output = 1, refComplex=refcom,
                           outFile=t.testRoot()+'/test.cl')

    master.start()

