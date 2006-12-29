##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

import os, tempfile, commands, copy
import Numeric as N

import Biskit.tools as t
import Biskit.settings as s
import Biskit.mathUtils as MU
from Biskit.LogFile import LogFile, StdLog
from Biskit.PDBModel import PDBModel
from Biskit.Errors import BiskitError
from Biskit.PDBCleaner import PDBCleaner
from Biskit import Executor

class AmberError( BiskitError ):
    pass

class AmberParmBuilder:
    """
    Create Amber topology and coordinate file from PDB.

    parmMirror():
    ...builds a fake parm that exactly mirrors a given PDB file.
      This parm can be used for ptraj but not for simulations.
      Currently, parmMirror only accepts amber-formatted PDBs as input.
      It should be possible to create topologies that have the same content
      and order of atoms as an xplor PDB but some atoms will have different
      names.

    parmSolvated():
    ...builds a solvated system for PME simulations (incl. closing of
      S-S bonds, capping of chain breaks).
      parmSolvated accepts both xplor and amber-formatted PDBs as input.

    Requires the amber programs tleap and ambpdb.
    Requires leap template files in /external/amber/leap/.

    Note: The design of AmberParmBuilder is less than elegant. It would make
    more sense to split it into two classes that are both derrived from
    Executor.
    """

    ## script to create a parm that exactly mirrors a given PDB
    script_mirror_pdb = """
    logFile %(leap_out)s
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
                  leaprc=s.leaprc,
                  leap_bin=s.tleap_bin,
                  leap_out=None, leap_in=None,
                  leap_pdb=None,
                  log=None,
                  debug=0,
                  verbose=0,
                  **kw ):
        """
        model   - PDBModel or str
        leap_template - str, path to template file for leap input [..]
        leaprc        - str, path to parameter file for leap [..]
        leap_out      - str, target file for leap.log [default: discard]
        leap_in       - str, target file for leap.in script [default: discard]
        leap_bin      - str, path to tleap binary [..]
        kw            - kw=value pairs for additional options in the
                        leap_template
        """
        self.m = PDBModel( model )

        self.leap_template = t.absfile( leap_template )
        self.leap_bin= leap_bin
        self.leaprc  = leaprc

        self.leap_in = leap_in or tempfile.mktemp( '_leap_in' )
        self.keep_leap_in = leap_in is not None

        self.leap_pdb = leap_pdb or tempfile.mktemp( '_leap_pdb' )
        self.keep_leap_pdb = leap_pdb is not None

        self.leap_out= leap_out or tempfile.mktemp( '_leap_log' )
        self.keep_leap_out = leap_out is not None

        self.log = log or StdLog()

        self.debug = debug
        self.verbose = verbose

        self.__dict__.update( kw )


    def __runLeap( self, in_script, norun=0, **kw ):
        """
        Create script file and run Leap.
        in_script - str, content of ptraj script with place holders
        norun     - 1|0, only create leap script                [0]
        **kw      - key=value pairs for filling place holders in script
        """
        ## create leap script
        try:
            ## use own fields and given kw as parameters for leap script
            d = copy.copy( self.__dict__ )
            d.update( kw )
            
            in_script = in_script % d
            f = open( self.leap_in, 'w')
            f.write( in_script )
            f.close()

            if self.verbose:
                self.log.add('leap-script: ')
                self.log.add( in_script )

        except IOError:
            raise AmberError('Could not create leap input file')
        except:
            raise AmberError('missing option for leap input file\n'+\
                             'available: %s' % (str( d.keys() ) ))

        ## run tleap
        args = '-f %s' % self.leap_in

        if not norun:
            self.exe = Executor('tleap', args, log=self.log,verbose=1,
                                catch_out=0)
            self.output, self.error, self.status = self.exe.run()

            if not os.path.exists( kw['out_parm'] ):
                raise AmberError, "tleap failed"

        ## clean up
            
        if not self.keep_leap_in and not self.debug:
            t.tryRemove( self.leap_in )
        if not self.keep_leap_out and not self.debug:
            t.tryRemove( self.leap_out)


    def parm2pdb( self, f_parm, f_crd, f_out, aatm=0 ):
        """
        Use ambpdb to build PDB from parm and crd.
        f_parm - str, existing parm file
        f_crd  - str, existing crd file
        f_out  - str, target file name for PDB
        """
##         cmd = '%s -p %s -aatm < %s > %s' % \
        args = '-p %s %s' % (f_parm, '-aatm'*aatm )

        x = Executor('ambpdb', args, f_in=f_crd, f_out=f_out,
                     log=self.log, verbose=1)

        output,error,status = x.run()

        if not os.path.exists( f_out ):
            raise AmberError, 'ambpdb failed.'

        return f_out
    

    def __ssBonds( self, model, cutoff=4. ):
        """
        -> [ (int, int) ], list with numbers of residue pairs forming S-S
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
                    r += [ (m.atoms[i]['residue_number'],
                            m.atoms[j]['residue_number']) ]
        return r


    def __cys2cyx( self, model, ss_residues ):
        """
        Rename all S-S bonded CYS into CYX.
        model - PDBModel
        ss_residues - [ ( int, int) ], original residue numbers of S-S pairs
        """
        ss = []
        for a,b in ss_residues:
            ss += [a,b]

        for a in model.atoms:
            if a['residue_number'] in ss:
                a['residue_name'] = 'CYX'


    def capACE( self, model, chain ):
        """
        Cap N-terminal of given chain.
        model - PDBModel
        chain - int, index of chain to be capped
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
        for a in m_ace.atoms:
            if a['residue_name'] != 'ACE':
                a['residue_name'] = m_term.atoms[0]['residue_name']
            else:
                a['residue_number'] = m_term.atoms[0]['residue_number']-1
                a['chain_id']       = m_term.atoms[0]['chain_id']
                a['segment_id']     = m_term.atoms[0]['segment_id']

        ## fit cap onto first residue of chain
        m_ace = m_ace.magicFit( m_term )

        ## concat cap on chain
        m_chain = m_ace.resModels()[0].concat( m_chain )

        ## re-assemble whole model
        return chains_before.concat( m_chain, chains_after )

        
    def capNME( self, model, chain ):
        """
        Cap C-terminal of given chain.
        model - PDBModel
        chain - int, index of chain to be capped
        """
        self.log.add('Capping C-terminal of chain %i.' % chain )
        m_nme   = PDBModel( self.F_nme_cap )

        chains_before = model.takeChains( range(chain), breaks=1 )
        m_chain       = model.takeChains( [chain], breaks=1 )
        chains_after  = model.takeChains( range(chain+1, model.lenChains(1)),
                                          breaks=1 )

        m_term  = m_chain.resModels()[-1]

        ## rename overhanging residue in cap PDB, renumber cap residue
        for a in m_nme.atoms:
            if a['residue_name'] != 'NME':
                a['residue_name'] = m_term.atoms[0]['residue_name']
            else:
                a['residue_number'] = m_term.atoms[0]['residue_number']+1
                a['chain_id']       = m_term.atoms[0]['chain_id']
                a['segment_id']     = m_term.atoms[0]['segment_id']

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
        center = N.average( model.getXyz() )
        model.setXyz( model.xyz - center )

    def leapModel( self, hetatm=0 ):
        """Get a clean PDBModel for input into leap."""
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
        f_out     - str, target file for parm (topology)
        f_out_crd - str, target file for crd (coordinates) [|f_out_base|.crd]
        f_out_pdb - str, target file for pdb               [|f_out_base|.pdb]
        hetatm    - 1|0, keep hetero atoms                                [0]
        cap       - 1|0, put ACE and NME capping residue on chain breaks  [0]
        capN      - [ int ], indices of chains that should get ACE cap   [[]]
        capC      - [ int ], indices of chains that should get NME cap   [[]]
        box       - float, minimal distance of solute from box edge    [10.0]
        fmod      - [ str ], list of files with amber parameter modifications
                    (to be loaded into leap with loadAmberParams)        [[]]
        fprep     - [ str ], list of files with amber residue definitions
                    (to be loaded into leap with loadAmberPrep)          [[]]
        **kw .. additional key=value pairs for leap input template
        """
        f_out = t.absfile( f_out )
        f_out_crd = t.absfile( f_out_crd ) or t.stripSuffix( f_out ) + '.crd'
        f_out_pdb = t.absfile( f_out_pdb ) or t.stripSuffix( f_out ) +\
                    '_leap.pdb'

        fmod  = [ t.absfile( f ) for f in t.toList( fmod )  ]
        fprep = [ t.absfile( f ) for f in t.toList( fprep ) ]

        try:
            self.log.add( 'Cleaning PDB file for Amber:' )
            m = self.leapModel( hetatm=hetatm )

            if cap:
                end_broken = m.atom2chainIndices( m.chainBreaks() )
                capC = MU.union( capC, end_broken )
                capN = MU.union( capN, N.array( end_broken ) + 1 )

            for i in capN:
                m = self.capACE( m, i )
            for i in capC:
                m = self.capNME( m, i )
            
            m.renumberResidues( addChainId=1 )  ## again, to accomodate capping

            template = open( self.leap_template ).read()

            leap_mod = self.__fLines( 'm = loadAmberParams %s\n', fmod )
            leap_prep= self.__fLines( 'loadAmberPrep %s\n', fprep )

            ss = self.__ssBonds( m, cutoff=4. )
            self.__cys2cyx( m, ss )
            leap_ss  = self.__fLines( self.ss_bond, ss )
            self.log.add('Found %i disulfide bonds: %s' % (len(ss),str(ss)))

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
        """ -> [ str ], leap statements for deleting given atoms """

        cmd = 'remove p.%(res)i p.%(res)i.%(atom)i'

        rM = m.resMap()
        rI = m.resIndex()

        s = []
        for i in i_atoms:
            res = m.atoms[i]['residue_number']
            ## substract atom index of first atom in of this residue
            atm = i - rI[ rM[ i ] ] + 1  
            s += [ cmd % {'res':res, 'atom':atm} ]

        return '\n'.join( s )


    def __inverseIndices( self, model, i_atoms ):
        """ -> [ int ], remaining atom indices of m that are NOT in i_atoms"""

        mask = N.zeros( len( model ),'i' )
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
        f_out     - str, target parm file
        f_out_crd - str, target crd file (default: f_out but ending .crd)
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
if __name__ == '__main__':

    f = t.testRoot() + '/Amber/'

    traj = t.Load( f + 'lig.etraj')
    m = traj.ref

    ## create solvated topology from xplor PDB
    a = AmberParmBuilder( m, verbose=1, debug=1 )

    a.parmSolvated( f + 'lig_solvated.parm', box=2.5,
                    f_out_pdb = f + 'lig_solvated.pdb',
                    f_out_crd = f + 'lig_solvated.crd')


    ## create amber pdb without water and hydrogen
    m = m.compress( m.maskProtein() )
    m = m.compress( m.maskHeavy() )

    f_in = f + '/AmberParmBuilder/lig_stripped.pdb'

    m.writePdb( f_in )

    ## create mirror parm for this stripped PDB

    a = AmberParmBuilder( f_in, verbose=1, debug=1 )

    a.parmMirror( f + '/AmberParmBuilder/lig_stripped.parm' )