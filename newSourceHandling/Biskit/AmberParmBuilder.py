## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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
## $Revision$
## last $Date$
## last $Author$

"""
Create Amber topology and coordinate file from PDB.
"""

import os, tempfile, copy
import numpy.oldnumeric as N

import Biskit.tools as t
import Biskit.settings as s
import Biskit.mathUtils as MU
from Biskit.LogFile import LogFile, StdLog
from Biskit.PDBModel import PDBModel
from Biskit.Errors import BiskitError
from Biskit.PDBCleaner import PDBCleaner
from Biskit.AmberLeap import AmberLeap
from Biskit import Executor

class AmberError( BiskitError ):
    pass

class AmberParmBuilder:
    """
    AmberParmBuilder
    ================
    Create Amber topology and coordinate file from PDB.

      - parmMirror():
         ...builds a fake parm that exactly mirrors a given PDB file.
         This parm can be used for ptraj but not for simulations.
         Currently, parmMirror only accepts amber-formatted PDBs as
         input. It should be possible to create topologies that have
         the same content and order of atoms as an xplor PDB but
         some atoms will have different names.

      - parmSolvated():
         ...builds a solvated system for PME simulations (incl. closing
         of S-S bonds, capping of chain breaks). parmSolvated accepts
         both xplor and amber-formatted PDBs as input.

    Requires the amber programs C{tleap} and C{ambpdb}.
    Requires leap template files in C{biskit/external/amber/leap/}.

    @note: The design of AmberParmBuilder is less than elegant. It
           would make more sense to split it into two classes that
           are both derrived from Executor.
    """

    ## script to create a parm that exactly mirrors a given PDB
    script_mirror_pdb = """
    logFile %(f_out)s
    source %(leaprc)s
    %(fmod)s
    %(fprep)s
    p = loadPdb %(in_pdb)s
    %(delete_atoms)s
    saveAmberParm p %(out_parm)s %(out_crd)s
    quit
    """

    ## tleap command to close a single S-S bond
    ss_bond = "bond p.%i.SG p.%i.SG\n"

    ## leap script for solvated topology
    F_leap_in = t.projectRoot() + '/external/amber/leap/solvate_box.leap'
    ## PDB with ACE capping residue
    F_ace_cap = t.projectRoot() + '/external/amber/leap/ace_cap.pdb'
    ## PDB with NME capping residue
    F_nme_cap = t.projectRoot() + '/external/amber/leap/nme_cap.pdb'

    def __init__( self, model,
                  leap_template=F_leap_in,
                  leaprc=None,
                  leap_out=None, leap_in=None,
                  leap_pdb=None,
                  log=None,
                  debug=0,
                  verbose=0,
                  **kw ):
        """
        @param model: model
        @type  model: PDBModel or str
        @param leap_template: path to template file for leap input
        @type  leap_template: str
        @param leaprc: forcefield parameter file or code (e.g. ff99)
        @type  leaprc: str
        @param leap_out: target file for leap.log (default: discard)
        @type  leap_out: str
        @param leap_in: target file for leap.in script (default: discard)
        @type  leap_in: str
        @param kw: kw=value pairs for additional options in the leap_template
        @type  kw: key=value
        """
        self.m = PDBModel( model )

        self.leap_template = leap_template
        self.leaprc  = leaprc

        self.leap_pdb = leap_pdb or tempfile.mktemp( '_leap_pdb' )
        self.keep_leap_pdb = leap_pdb is not None

        self.leap_in = leap_in
        self.leap_out= leap_out

        self.log = log or StdLog()

        self.debug = debug
        self.verbose = verbose

        self.__dict__.update( kw )


    def __runLeap( self, in_script, in_pdb, norun=0, **kw ):
        """
        Create script file and run Leap.

        @param in_script: content of ptraj script with place holders
        @type  in_script: str
        @param in_pdb: PDB file to load into tleap
        @type  in_pdb: str
        @param norun: 1 - only create leap scrip (default: 0)
        @type  norun: 1|0
        @param kw: key=value pairs for filling place holders in script
        @type  kw: key=value

        @raise AmberError: if missing option for leap input file or
                           if could not create leap input file
        """
        x = AmberLeap( in_script,
                       in_pdb=in_pdb,
                       log=self.log, verbose=self.verbose, debug=self.debug,
                       catch_out=True,
                       f_in=self.leap_in,
                       f_out=self.leap_out,
                       **kw )
        if norun:
            x.generateInp()
        else:
            x.run()
        
##         ## create leap script
##         try:
##             ## use own fields and given kw as parameters for leap script
##             d = copy.copy( self.__dict__ )
##             d.update( kw )

##             in_script = in_script % d
##             f = open( self.leap_in, 'w')
##             f.write( in_script )
##             f.close()

##             if self.verbose:
##                 self.log.add('leap-script: ')
##                 self.log.add( in_script )

##         except IOError:
##             raise AmberError('Could not create leap input file')
##         except:
##             raise AmberError('missing option for leap input file\n'+\
##                              'available: %s' % (str( d.keys() ) ))

##         ## run tleap
##         args = '-f %s' % self.leap_in

##         if not norun:
##             self.exe = Executor('tleap', args, log=self.log,verbose=1,
##                                 catch_out=0)
##             self.output, self.error, self.status = self.exe.run()

##             if not os.path.exists( kw['out_parm'] ):
##                 raise AmberError, "tleap failed"

##         ## clean up

##         if not self.keep_leap_in and not self.debug:
##             t.tryRemove( self.leap_in )
##         if not self.keep_leap_out and not self.debug:
##             t.tryRemove( self.leap_out)


    def parm2pdb( self, f_parm, f_crd, f_out, aatm=0 ):
        """
        Use ambpdb to build PDB from parm and crd.

        @param f_parm: existing parm file
        @type  f_parm: str
        @param f_crd: existing crd file
        @type  f_crd: str
        @param f_out: target file name for PDB
        @type  f_out: str

        @return: f_out, target file name for PDB
        @rtype: str

        @raise AmberError: if ambpdb fail
        """
##         cmd = '%s -p %s -aatm < %s > %s' % \
        args = '-p %s %s' % (f_parm, '-aatm'*aatm )

        x = Executor('ambpdb', args, f_in=f_crd, f_out=f_out,
                     log=self.log, verbose=1, catch_err=1)

        output,error,status = x.run()

        if not os.path.exists( f_out ):
            raise AmberError, 'ambpdb failed.'

        return f_out


    def __ssBonds( self, model, cutoff=4. ):
        """
        Identify disulfide bonds.

        @param model: model
        @type  model: PDBModel        
        @param cutoff: distance cutoff for S-S distance (default: 4.0)
        @type  cutoff: float
        
        @return: list with numbers of residue pairs forming S-S
        @rtype: [(int, int)]
        """
        m = model.compress( model.mask( ['SG'] ) )

        if len( m ) < 2:
            return []

        pw = MU.pairwiseDistances( m.xyz, m.xyz )

        pw = N.less( pw, cutoff )

        r = []
        for i in range( len( pw ) ):
            for j in range( i+1, len(pw) ):
                if pw[i,j]:
                    r += [ (m.atoms['residue_number'][i],
                            m.atoms['residue_number'][j]) ]
        return r


    def __cys2cyx( self, model, ss_residues ):
        """
        Rename all S-S bonded CYS into CYX.

        @param model: model
        @type  model: PDBModel
        @param ss_residues: original residue numbers of S-S pairs
        @type  ss_residues: [(int, int)]
        """
        ss = []
        for a,b in ss_residues:
            ss += [a,b]

        for a in model:
            if a['residue_number'] in ss:
                a['residue_name'] = 'CYX'


    def capACE( self, model, chain ):
        """
        Cap N-terminal of given chain.

        @param model: model
        @type  model: PDBMode
        @param chain: index of chain to be capped
        @type  chain: int
        """
        self.log.add('Capping N-terminal of chain %i.' % chain )
        m_ace = PDBModel( self.F_ace_cap )

        chains_before = model.takeChains( range(chain), breaks=1 )
        m_chain       = model.takeChains( [chain], breaks=1 )
        chains_after  = model.takeChains( range(chain+1, model.lenChains(1)),
                                          breaks=1 )

        m_term  = m_chain.resModels()[0]

        ## we need 3 atoms for superposition, CB might mess things up but
        ## could help if there is no HN
        if 'HN' in m_term.atomNames():
            m_ace.remove( ['CB'] )

        ## rename overhanging residue in cap PDB
        for a in m_ace:
            if a['residue_name'] != 'ACE':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]-1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## fit cap onto first residue of chain
        m_ace = m_ace.magicFit( m_term )

        ## concat cap on chain
        m_chain = m_ace.resModels()[0].concat( m_chain )

        ## re-assemble whole model
        return chains_before.concat( m_chain, chains_after )


    def capNME( self, model, chain ):
        """
        Cap C-terminal of given chain.

        @param model: model
        @type  model: PDBMode
        @param chain: index of chain to be capped
        @type  chain: int        
        """
        self.log.add('Capping C-terminal of chain %i.' % chain )
        m_nme   = PDBModel( self.F_nme_cap )

        chains_before = model.takeChains( range(chain), breaks=1 )
        m_chain       = model.takeChains( [chain], breaks=1 )
        chains_after  = model.takeChains( range(chain+1, model.lenChains(1)),
                                          breaks=1 )

        m_term  = m_chain.resModels()[-1]

        ## rename overhanging residue in cap PDB, renumber cap residue
        for a in m_nme:
            if a['residue_name'] != 'NME':
                a['residue_name'] = m_term.atoms['residue_name'][0]
            else:
                a['residue_number'] = m_term.atoms['residue_number'][0]+1
                a['chain_id']       = m_term.atoms['chain_id'][0]
                a['segment_id']     = m_term.atoms['segment_id'][0]

        ## chain should not have any terminal O after capping
        m_chain.remove( ['OXT'] )            

        ## fit cap onto last residue of chain
        m_nme = m_nme.magicFit( m_term )

        ## concat cap on chain
        m_chain = m_chain.concat( m_nme.resModels()[-1] )

        if m_chain._PDBModel__terAtoms != []:
            m_chain._PDBModel__terAtoms = [ len( m_chain ) - 1 ]

        ## re-assemble whole model
        return chains_before.concat( m_chain, chains_after )


    def centerModel( self, model ):
        """
        Geometric centar of model.
        
        @param model: model
        @type  model: PDBMode
        """
        center = N.average( model.getXyz() )
        model.setXyz( model.xyz - center )


    def leapModel( self, hetatm=0 ):
        """
        Get a clean PDBModel for input into leap.

        @param hetatm: keep HETATM records (default: 0)
        @type  hetatm: 1|0

        @return: model
        @rtype: PDBMod
        """
        m = self.m.clone()
        m.xplor2amber()

        cleaner = PDBCleaner( m, log=self.log )
        m = cleaner.process( keep_hetatoms=hetatm, amber=1 )

        m.renumberResidues( addChainId=1 )

        self.centerModel( m )

        return m


    def __fLines( self, template, values ):
        if not type( values ) is list:
            values = [ values ]

        return ''.join( [ template % v for v in values ] )


    def parmSolvated( self, f_out, f_out_crd=None, f_out_pdb=None,
                      hetatm=0, norun=0,
                      cap=0, capN=[], capC=[],
                      fmod=[], fprep=[],
                      box=10.0, **kw ):
        """
        @param f_out: target file for parm (topology)
        @type  f_out: str
        @param f_out_crd: target file for crd (coordinates)
                          (default:|f_out_base|.crd)
        @type  f_out_crd: str
        @param f_out_pdb: target file for pdb (default:|f_out_base|.pdb)
        @type  f_out_pdb: str
        @param hetatm: keep hetero atoms (default: 0)
        @type  hetatm: 1|0
        @param cap: put ACE and NME capping residue on chain breaks 
                    (default: 0)
        @type  cap: 1|0
        @param capN: indices of chains that should get ACE cap (default: [])
        @type  capN: [int]
        @param capC: indices of chains that should get NME cap (default: [])
        @type  capC: [int]
        @param box: minimal distance of solute from box edge (default: 10.0)
        @type  box: float
        @param fmod: list of files with amber parameter modifications
                    (to be loaded into leap with loadAmberParams) (default:[])
        @type  fmod: [str]
        @param fprep: list of files with amber residue definitions
                    (to be loaded into leap with loadAmberPrep) (default: [])
        @type  fprep: [str]
        @param kw: additional key=value pairs for leap input template
        @type  kw: key=value

        @raise IOError:
        """
        f_out = t.absfile( f_out )
        f_out_crd = t.absfile( f_out_crd ) or t.stripSuffix( f_out ) + '.crd'
        f_out_pdb = t.absfile( f_out_pdb ) or t.stripSuffix( f_out ) +\
                    '_leap.pdb'

        fmod  = [ t.absfile( f ) for f in t.toList( fmod )  ]
        fprep = [ t.absfile( f ) for f in t.toList( fprep ) ]

        try:
            if self.verbose: self.log.add( '\nCleaning PDB file for Amber:' )
            m = self.leapModel( hetatm=hetatm )

            if cap:
                end_broken = m.atom2chainIndices( m.chainBreaks() )
                capC = MU.union( capC, end_broken )
                capN = MU.union( capN, N.array( end_broken ) + 1 )

            for i in capN:
                if self.verbose:
                    self.log.add( 'Adding ACE cap to chain %i' % i )
                m = self.capACE( m, i )

            for i in capC:
                if self.verbose:
                    self.log.add( 'Adding NME cap to chain %i' % i )
                m = self.capNME( m, i )

            m.renumberResidues( addChainId=1 )  ## again, to accomodate capping

            template = open( self.leap_template ).read()

            leap_mod = self.__fLines( 'm = loadAmberParams %s\n', fmod )
            leap_prep= self.__fLines( 'loadAmberPrep %s\n', fprep )

            ss = self.__ssBonds( m, cutoff=4. )
            self.__cys2cyx( m, ss )
            leap_ss  = self.__fLines( self.ss_bond, ss )
            if self.verbose:
                self.log.add('Found %i disulfide bonds: %s' % (len(ss),str(ss)))

            if self.verbose:
                self.log.add( 'writing cleaned PDB to %s'  % self.leap_pdb )
            m.writePdb( self.leap_pdb, ter=3 )

            self.__runLeap( template, in_pdb=self.leap_pdb,
                            out_parm=f_out, out_crd=f_out_crd,
                            ss_bonds=leap_ss, fmod=leap_mod,
                            fprep=leap_prep, norun=norun,
                            box=box, **kw )

            if not norun:
                parm_pdb = self.parm2pdb( f_out, f_out_crd, f_out_pdb )

            if not self.keep_leap_pdb and not self.debug:
                t.tryRemove( self.leap_pdb )

        except IOError, why:
            raise IOError, why


    def __deleteAtoms( self, m, i_atoms ):
        """
        Delet atoms

        @param m: model
        @type  m: PDBMode
        @param i_atoms: atom index
        @type  i_atoms: [int]
        
        @return: leap statements for deleting given atoms
        @rtype: [str]
        """
        cmd = 'remove p.%(res)i p.%(res)i.%(atom)i'

        rM = m.resMap()
        rI = m.resIndex()

        s = []
        for i in i_atoms:
            res = m.atoms['residue_number'][i]
            ## substract atom index of first atom in of this residue
            atm = i - rI[ rM[ i ] ] + 1  
            s += [ cmd % {'res':res, 'atom':atm} ]

        return '\n'.join( s )


    def __inverseIndices( self, model, i_atoms ):
        """
        @param model: model
        @type  model: PDBMode
        @param i_atoms: atom index
        @type  i_atoms: [int]
    
        @return: remaining atom indices of m that are NOT in i_atoms
        @rtype: [int]
        """
        mask = N.zeros( len( model ),N.Int )
        N.put( mask, i_atoms, 1 )
        return N.nonzero( N.logical_not( mask ) )


    def parmMirror( self, f_out, f_out_crd=None, fmod=[], fprep=[], **kw ):
        """
        Create a parm7 file whose atom content (and order) exactly mirrors
        the given PDBModel. This requires two leap runs. First we get a
        temporary topology, then we identify all atoms added by leap and
        build a final topology where these atoms are deleted.
        This parm is hence NOT suited for simulations but can be used to parse
        e.g. a trajectory or PDB into ptraj.

        @param f_out: target parm file
        @type  f_out: str
        @param f_out_crd: target crd file (default: f_out but ending .crd)
        @type  f_out_crd: str
        """
        f_out = t.absfile( f_out )
        f_out_crd = t.absfile( f_out_crd ) or t.stripSuffix( f_out ) + '.crd'

        ## if there are hydrogens, recast them to standard amber names
        aatm = 'HA' in self.m.atomNames() ## 'HB2' in self.m.atomNames()

        ## First leap round ##
        m_ref = self.m.clone()
        m_ref.xplor2amber( aatm=aatm )
        tmp_in = tempfile.mktemp( 'leap_in0.pdb' )
        m_ref.writePdb( tmp_in, ter=3 )

        tmp_parm = tempfile.mktemp( '_parm0' )
        tmp_crd  = tempfile.mktemp( '_crd0' )

        leap_mod = self.__fLines( 'm = loadAmberParams %s\n', fmod )
        leap_prep= self.__fLines( 'loadAmberPrep %s\n', fprep )

        self.__runLeap( self.script_mirror_pdb,
                        leaprc=self.leaprc, fmod=leap_mod, fprep=leap_prep,
                        in_pdb=tmp_in, out_parm=tmp_parm, out_crd=tmp_crd,
                        delete_atoms='' )

        tmp_pdb = self.parm2pdb( tmp_parm, tmp_crd,
                                 tempfile.mktemp( 'leap_out.pdb' ), aatm=aatm )

        if not self.debug:
            t.tryRemove( tmp_parm )
            t.tryRemove( tmp_crd )
            t.tryRemove( tmp_in )

        ## load model with missing atoms added by leap
        m_leap = PDBModel( tmp_pdb  )

        ## compare atom content
        iLeap, iRef = m_leap.compareAtoms( m_ref )

        ## check that ref model doesn't need any change
        if iRef != range( len( m_ref ) ):
            raise AmberError, "Cannot create exact mirror of %s.\n" % tmp_in +\
                  "Leap has renamed/deleted original atoms in %s."% tmp_pdb

        ## indices of atoms that were added by leap
        delStr = self.__deleteAtoms( m_leap,
                                     self.__inverseIndices( m_leap, iLeap ) )

        ## Second leap round ##
        self.__runLeap( self.script_mirror_pdb, leaprc=self.leaprc,
                        in_pdb=tmp_pdb, fmod=leap_mod, fprep=leap_prep,
                        out_parm=f_out, out_crd=f_out_crd,
                        delete_atoms=delStr )

        if not self.debug:
            t.tryRemove( tmp_pdb )


#############
## TESTING ##
import Biskit.test as BT
import tempfile
import Biskit.tools as T

class Test( BT.BiskitTest ):
    """Test AmberParmBuilder"""

    TAGS = [ BT.EXE ]

    def prepare(self):
        root = T.testRoot() + '/amber/'
        self.ref = PDBModel( T.testRoot() + '/amber/1HPT_0.pdb')
        self.refdry = root + '1HPT_0dry.pdb'

        self.dryparm = tempfile.mktemp('.parm', 'dry_')
        self.drycrd  = tempfile.mktemp('.crd', 'dry_')
        self.drypdb  = tempfile.mktemp('.pdb', 'dry_')
        self.wetparm = tempfile.mktemp('.parm', 'wet_')
        self.wetcrd  = tempfile.mktemp('.crd', 'wet_')
        self.wetpdb  = tempfile.mktemp('.pdb', 'wet_')

    def cleanUp(self):
        T.tryRemove( self.dryparm )
        T.tryRemove( self.drycrd )
        T.tryRemove( self.drypdb )
        T.tryRemove( self.wetparm )
        T.tryRemove( self.wetcrd )
        T.tryRemove( self.wetpdb )
        

    def test_AmberParmMirror(self):
        """AmberParmBuilder.parmMirror test"""
        ref = self.ref
        mask = N.logical_not( ref.maskH2O() ) ## keep protein and Na+ ion
        self.mdry = ref.compress( mask )

        self.a = AmberParmBuilder( self.mdry, verbose=self.local,
                                   debug=self.DEBUG )

        self.a.parmMirror(f_out=self.dryparm,
                          f_out_crd=self.drycrd )

        self.a.parm2pdb( self.dryparm, self.drycrd, self.drypdb )

        self.m1 = PDBModel(self.drypdb)
        self.m2 = PDBModel(self.refdry)

        eq = N.array( self.m1.xyz == self.m2.xyz )
        self.assert_( eq.all() )


    def test_AmberParmSolvated( self ):
        """AmberParmBuilder.parmSolvated test"""
        ## remove waters and hydrogens
        self.mdry = self.ref.compress( self.ref.maskProtein() )
        self.mdry = self.mdry.compress( self.mdry.maskHeavy() )

        self.a = AmberParmBuilder( self.mdry,
                                   verbose=self.local, debug=self.DEBUG)

        self.a.parmSolvated( self.wetparm, f_out_crd=self.wetcrd,
                             f_out_pdb=self.wetpdb,
                             box=2.5 )

        self.m3 = PDBModel( self.wetpdb )

        m3prot = self.m3.compress( self.m3.maskProtein() )
        refprot= self.ref.compress( self.ref.maskProtein() )
        
        self.assertEqual( self.ref.lenChains(), self.m3.lenChains() )
        self.assertEqual( refprot.atomNames(), m3prot.atomNames() )



if __name__ == '__main__':

    BT.localTest()
