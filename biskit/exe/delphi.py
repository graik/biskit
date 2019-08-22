##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg
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

"""
Wrapper for Delphi -- Poisson-Boltzman electrostatic potential calculation.
"""

import tempfile, os
import numpy as N
import re
import operator

import biskit.tools as T
import biskit.mathUtils as U

from biskit import PDBModel, AmberResidueLibrary, AmberPrepParser, AtomCharger

from biskit.exe.executor import Executor
from biskit.exe.reduce import Reduce


class DelphiError( Exception ):
    pass

class PDB2DelphiCharges( object ):
    """
    Generate a Delphi Charge file from the atomic partial charges assigned
    to a Input structure. Charge records are built as generically as possible.
    That means without residue number or chain ID, in order to minimize the
    number of records. 
    
    Cases where the same name is used for residues with different atom
    content, e.g. C-terminal ALA and normal ALA, require special treatment:
    The most common residue is described by a generic record (ALA, no residue
    number) and any additional versions with different atom content are
    described with explicit residue number and chain ID.
    
    """
    
    def __init__( self, model ):
        """
        """
        self.model = model
        self.resmap = None
        
    def prepare( self ):
        if not 'partial_charge' in self.model.atoms:
            ac = AtomCharger()
            ac.charge( self.model )
        
        self.resmap = self.mapResidues()
        
    def mapResidues( self ):
        """
        """
        resmap = {}
        resmodels = self.model.resModels()
        
        for res in resmodels:
            resname = res['residue_name'][0]
            akey = res.atomkey(compress=False)
            ## initialize dict of dict if not existing
            resmap[ resname ] = resmap.get(resname,{})
            resmap[ resname ][akey] = resmap[ resname ].get(akey, [])
            resmap[ resname ][ akey ]+= [ res ]

        return resmap

    def res2delphi(self, restype, resnumber=None, chain=None, comment='' ):
        """
        Generate the Delphi charge record for a single residue (residue type).
        @type restype: AmberResidueType
        @param resnumber: residue number to use in the Delphi record
        @type  resnumber: int
        @return: several lines formatted as in Delphi charge files
        @rtype: str
        """
        resnumber = resnumber or ''
        chain = chain or ''

        r = ''
        for atom in restype.atoms.iterDicts():
 
            atom['resnumber'] = str(resnumber)
            atom['chain'] = chain
            r += '%(name)-5s %(residue_name)3s %(resnumber)-3s%(chain)1s  %(partial_charge)6.3f'\
                  % atom

            if comment:
                r += '  ! %s' % comment
                comment = ''  # only print once
            r += '\n'
        return r
    
    def manyres2delphi(self, akeydict ):
        """
        Generate one or several charge records for one or more residues with
        the same name (e.g. normal ALA and C-terminal ALA).
        
        @param akeydict: dict of residues indexed by atom keys
        @type  akeydict: { str : PDBModel }
        @return: Delphi Charge record for one or more residues with the same
                 name
        """
        if len( akeydict ) == 1:
            return self.res2delphi( list(akeydict.values())[0][0] )
        
        ## determine which of several residues is used most
        keyres = [ (len(reslist), reslist) for reslist in akeydict.values() ]
        keyres.sort(key=operator.itemgetter(0))  ## force sorting by 1st item only
        keyres.reverse()

        ## one record w/o number and chain id for all residues with these atoms
        r = ''
        key, res = keyres[0]
        r += self.res2delphi( res[0] )
        
        ## now add less used residues with number and chainID
        for key, residues in keyres[1:]:
            for res in residues:
                r += self.res2delphi( res, resnumber=res['residue_number'][0],
                                      chain=res['chain_id'][0] )
            
        return r
    
    
    def tofile( self, fname ):
        if not self.resmap:
            self.prepare()
        
        fname = T.absfile( fname )
        f = open( fname, 'w' )
        f.write('! charge file generated by Biskit.delphi.DelphiCharges\n')
        f.write('atom__resnumbc_charge_\n')

        try:
            for resname, akeys in self.resmap.items():
                f.write( self.manyres2delphi( akeys ) )
        finally:
            f.close()
   
    
class Delphi( Executor ):
    """
    Calculate electrostatic potentials and potential maps with Delphi.
    
    The current default workflow of this wrapper is:
    1. take input model/structure, remove hydrogens
    2. optionally: cap chain breaks and probable false N- and C-termini with 
       NME and ACE residues to cancel artificial charges 
       (see L{Biskit.PDBCleaner})
    3. add and optimize hydrogens with the reduce program
       (see L{Biskit.Reduce} )
    4. adapt residue and atom names to Amber conventions
    5. match each residue (by atom content) to a residue from a list of 
       Amber residue topology files and assign partial charge to each atom
       (see L{Biskit.AtomCharger})
    6. Create custom delphi charge file with Amber partial charges
    7. Determine center and dimensions of the grid used in Delphi calculation
    8. Run Delphi in temporary folder 
    9. parse result energies into result dictionary
    
    Usage
    =====

    >>> D = Delphi( inputmodel )
    >>> result = D.run()
    >>> result
    {'scharge' :  1.4266   # surface charge
     'egrid' :  9105.51    # total grid energy
     'ecoul' :  -9849.664  # couloumb energy
     'erxn'  :  -664.7469   # corrected reaction field energy
     'erxnt' :  -21048.13  # total reaction field energy
     'eself' :  -20383.39 }  # self reaction field energy
    
    Note
    ====
    
    All energy values are in units of kT. The important terms for the
    calculation of (free) energy differences are 'egrid', 'erxn' and 'ecoul'.
    
    Grid dimensions
    ===============
    
    If no further options are given, the Delphi wrapper will determine grid
    dimensions that are centered on the molecule and lead to at most 60%
    "filling" of the grid in any of the x, y or z dimensions. In other words,
    the longest dimension of the solute along x, y or z axis will still have a
    20% margin to the grid boundary. This fill-factor can be overriden with
    the 'perfil' parameter.
    
    The density of the grid is, by default, 2.3 points per Angstroem
    and can be adjusted with the 'scale' parameter. The number of grid points 
    (gsize) will be calculated accordingly. 
    
    For a given combination of scale and perfil, this calculation will give
    approximately (but not exactly) the same grid dimensions as the one the
    delphi program would determine by itself. Minor differences may lead to
    two grid points more or less. This can also lead to some changes in
    resulting energies.
    
    For more control over grid dimensions, you should call the method
    Delphi.setGrid() *before* Delphi.run(). For example:
    
    >>> D = Delphi( inputmodel )
    >>> D.setGrid( scale=1.2, perfil=80 )
    
    ... will fix a grid centered on the molecule, with 1/1.2 A grid spacing and 
    at least 10% distance margin to the boundary. Another example:
    
    >>> D = Delphi( inputmodel )
    >>> D.setGrid( acenter=(0,0,0), scale=1.0, gsize=100 )
    
    ... circumvents any calculation and fixes a 100 x 100 x 100 grid centered 
    on 0 and with 1/1.0 A grid spacing -- a box of 100 A x 100 A x 100 A. If
    the scale parameter is not given, it will default to whatever was specified
    at the Delphi() constructor (and from there default to 2.3). 
    
    Delphi.setGrid() returns the three grid parameters as a dictionary for 
    later re-use. So if you want to use the same grid dimensions for two 
    different structures you can first calculate the dimensions based on one
    structure and then apply the same dimensions to another Delphi run:
    
    >>> D1 = Delphi( complex )
    >>> grid1 = D1.setGrid( scale=2.4, perfil=60 )
    >>> print grid1
    {'acenter': (0.1, 10., -0.5), 'scale': 2.400, 'gsize': 99 }
    
    >>> ligand = complex.takeChains( [1] ) ## extract part of structure
    >>> D2 = Delphi( ligand )
    >>> D2.setGrid( **grid1 )
    
    Note: The '**' syntax unpacks the dictionary as keyword parameters into the 
    method.
    
    
    Customization
    =============
    
    The default delphi parameter file can be replaced (parameter template).
    Note though that you should keep the place holders for input and output
    files. The default parameter file is taken from:

        Biskit/data/delphi/delphi_simple.prm
        
    As always, the place holders (e.g. %(salt)f ) in this file are replaced
    by the value of a variable of the same name (e.g. salt) within the name
    space of the Executor instance. In other words:
    
    >>> D = Delphi( inputmodel, salt=0.2, ionrad=2.5 )
    
    ...will override the default values for salt and ionradius for which there
    are the place hoders %(salt)f and %(ionrad)f in the template file.
    
    
    The default handling and matching of atomic partial charges can be 
    modified by:
    
        * providing a ready-made Delphi charge file (parameter: f_charges)
        * or providing an alternative list of Amber topology files from which 
          residues are looked up by their atom content (parameter: topologies)
        * or circumvent charge matching (parameter: addcharge=False) and 
          provide a input PDBModel with an atom profile 'partial_charge' that 
          is used instead

    @note: Command configuration: biskit/Biskit/data/defaults/exe_delphi.dat
    """

    F_RADII = 'radii_amber_fromdelphi.siz'  ## default atom radius file
    F_PARAMS = 'delphi_simple.prm'    ## default delphi parameter file
    
    ## list of Amber topology files in decending priority
    ## will be mined for matching residues to assign charges
    F_RESTYPES = ['all_amino03.in',
                  'all_aminoct03.in',
                  'all_aminont03.in',
                  'all_nuc02.in' ]
    
    RE_E_GRID = r'total grid energy\s+:\s+(?P<egrid>[0-9\-\.]+)\s+kt'
    RE_E_COUL = r'coulombic energy\s+:\s+(?P<ecoul>[0-9\-\.]+)\s+kt'
    RE_E_SELF = r'self-reaction field energy\s+:\s+(?P<eself>[0-9\-\.]+)\s+kt'
    RE_E_RXN  = r'corrected reaction field energy\s*:\s+(?P<erxn>[0-9\-\.]+)\s+kt'
    RE_E_RXNT = r'total reaction field energy\s*:\s+(?P<erxnt>[0-9\-\.]+)\s+kt'
    RE_SURFCH = r'total s\.charge\,no epsin carrying\s*:\s+(?P<scharge>[0-9\-\.]+)'
    
    
    def __init__( self, 
                  model, 
                  template=None, 
                  topologies=None,
                  f_charges=None,
                  f_radii=None,
                  f_map=None,
                  addcharge=True,
                  protonate=True,
                  autocap=False,
                  indi=4.0, exdi=80.0, salt=0.15, ionrad=2, prbrad=1.4, 
                  bndcon=4, scale=2.3, perfil=60, 
                  **kw ):
        """
        @param model: structure for which potential should be calculated
        @type  model: PDBModel
        @param template: delphi command file template [None=use default]
        @type  template: str
        @param f_radii: alternative delphi atom radii file [None=use default]
        @type  f_radii: str
        @param topologies: alternative list of residue charge/topology files
                           [default: amber/residues/all*]
        @type  topologies: [ str ]
        @param f_charges: alternative delphi charge file 
                          [default: create custom]
        @type  f_charges: str
        @param f_map   : output file name for potential map [None= discard]
        @type  f_map   : str
        @param addcharge: build atomic partial charges with AtomCharger
                          [default: True]
        @type  addcharge: bool
        
        @param protonate: (re-)build hydrogen atoms with reduce program (True)
                          see L{Biskit.Reduce}
        @type  protonate: bool
        @param autocap: add capping NME and ACE residues to any (auto-detected)
                        false N- or C-terminal and chain breaks (default: False)
                        see L{Biskit.Reduce} and L{Biskit.PDBCleaner}
        @type  autocap: bool

        @param indi: interior dilectric (4.0)
        @param exdi: exterior dielectric (80.0)
        @param salt: salt conc. in M (0.15)
        @param ionrad: ion radius (2)
        @param prbrad: probe radius (1.4) 
        @param bndcon: boundary condition (4, delphi default is 2)
        @param scale:  grid spacing (2.3)
        @param perfil: grid fill factor in % (for automatic grid, 60) 
        
        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        template = template or os.path.join(T.dataRoot(),'delphi',self.F_PARAMS)
        
        tempdir = self.newtempfolder( tempdir=True )  ## create new temp folder
        f_in = tempfile.mktemp( '.inp', 'delphi_', dir=tempdir )
        
        self.f_pdb = tempfile.mktemp( '.pdb', 'delphi_', dir=tempdir)

        self.keep_map = f_map is not None
        self.f_map = f_map or \
            tempfile.mktemp( '_mapout.phi', 'delphi_', dir=tempdir )

##        self.f_map = None
        self.f_radiisrc = os.path.join(T.dataRoot(),'delphi',self.F_RADII)
        self.f_radii = f_radii or os.path.join(tempdir, 'radii.siz')
        
        self.topologies = topologies or self.F_RESTYPES
        self.f_charges = f_charges or tempfile.mktemp( '.crg', 'delphi_',
                                                       dir=tempdir )
        
        self.protonate = protonate
        self.autocap = autocap
        self.addcharge = addcharge
        
        ## DELPHI run parameters
        self.indi=indi  # interior dilectric(4.0)
        self.exdi=exdi  # exterior dielectric(80.0)
        self.salt=salt  # salt conc. in M (0.15)
        self.ionrad=ionrad # ion radius (2)
        self.prbrad=prbrad # probe radius (1.4) 
        self.bndcon=bndcon # boundary condition (4, delphi default is 2)
        
        ## DELPHI parameters for custom grid
        self.scale=scale   # grid spacing (2.3)
        self.perfil=perfil # grid fill factor in % (for automatic grid, 60)
        self.gsize = None
        self.acenter = None
        self.strcenter = '(0.0,0.0,0.0)'
        
        kw['tempdir'] = tempdir
        kw['cwd']     = tempdir
        
        Executor.__init__( self, 'delphi', 
                           template=template,
                           f_in=f_in,
                           args=f_in,
                           catch_err=True,
                           **kw )
        
        self.model = model
        self.delphimodel = None
    

    def delphiDimensions( self, model ):
        """
        Calculate "geometric" center and molecular dimensions as defined by 
        Delphi (the delphi geometric center is NOT exactly what a geometric
        center is commonly defined as). 
        @param model: PDBModel for which center and dimensions should be
                      calculated
        @type model:  PDBModel
        @return: center and dimensions
        @rtype : N.array([x,y,z] of float), N.array([x,y,z] of float)
        """
        m = model.compress( model.maskHeavy() )
        xyz = m.getXyz()
        
        ## largest protein length in x, y, z direction
        dimensions = N.max( xyz, 0 ) - N.min( xyz, 0 )
        
        center = N.min( xyz, 0) + dimensions / 2
        
        ## + 1 C diameter
        dimensions += 2*1.7
        
        return center, dimensions
        

    def setGrid( self, acenter=None, gsize=None, scale=None, perfil=None ):
        """
        Specify or calculate Delphi grid. There are two options:
        (1) specify the actual grid by giving acenter, scale and gsize.
        (2) calculate new grid dimensions from acenter, scale and perfil.
        If not given, acenter defaults to the geometric center of the structure
        model.
        @param acenter: center coordinates for the grid
        @type  acenter: [float, float, float]
        @param gsize: number of grid points in x, y, and z direction
        @type  gsize: int
        @param scale: distance between grid points
        @type scale : float
        @param perfil: percent fill factor 
        """
        center, dimensions = self.delphiDimensions( self.model )

        if acenter is not None:
            self.acenter = acenter
        else:
            if self.acenter is None:
                self.acenter = center
        
        self.strcenter = str( tuple( self.acenter) )
        self.scale= scale or self.scale
        
        if not gsize:
            self.perfil = perfil or self.perfil
              
            ## grid size in number of points at self.scale per Angstrom density
            gsize = (N.max( dimensions ) * 100. / self.perfil) * self.scale
            gsize = int( round( gsize ) )
            
            ## grid size must be an uneven number
            if not gsize % 2:
                gsize += 1
                
        self.gsize = gsize
        
        return self.getGrid()
    
    def getGrid( self ):
        return {'acenter':self.acenter, 'scale':self.scale, 'gsize':self.gsize}
        

    def __prepareRadii( self ):
        """
        Link default atom radius file into working directory.
        """
        try:
            os.symlink( self.f_radiisrc, self.f_radii)

        except OSError as error:
            raise DelphiError('Error preparing atom radius file for Delphi\n'+\
                                  'Error: %r\n' % error +\
                                  'radii file: %r\n' % self.f_radiisrc +\
                                  'target file: %r\n' % self.f_radii)
        

    def __prepareCharges(self, f_out ):

        try:
            if self.addcharge:
                if self.verbose:
                    self.log.add(
                        '\nAssigning atomic charges with AtomCharger...')
    
                if self.topologies is None:
                    reslib = None
                else:
                    reslib = AmberResidueLibrary(self.topologies,log=self.log,
                                                 verbose=self.verbose)
                ac = AtomCharger( reslibrary=reslib, 
                                  log=self.log, verbose=self.verbose )
                ac.charge( self.delphimodel )
            
            if self.verbose:
                self.log.add('Creating Delphi charge file %s'%f_out)

            dc = PDB2DelphiCharges( self.delphimodel )
            dc.prepare()
            dc.tofile( f_out )
            
            if self.verbose:
                qmissing = self.delphimodel['partial_charge']==0
                self.log.add('\nAtoms without charges: %i' % N.sum(qmissing))
                if N.sum(qmissing) > 0:
                    self.log.add('Warning: there are atoms without charge:')
                    m = self.delphimodel.compress( qmissing )
                    for a in m:
                        self.log.add(
                            '%(serial_number)4i %(name)-4s %(residue_name)3s %(residue_number)3i %(chain_id)s'\
                            % a)

        except IOError as why: 
            raise IOError('Error creating custom delphi charge file '+f_out+\
                  '( '+str(why)+' )')

    
    def prepare( self ):
        """
        Overrides Executor method.
        """
        Executor.prepare( self )
                
        ## if setGrid hasn't been called yet, create automatic grid
        if not self.gsize:
            self.setGrid()
        
        if self.protonate:
            reducer = Reduce( self.model, verbose=self.verbose,
                              autocap=self.autocap,
                              tempdir=self.tempdir, cwd=self.cwd,
                              log=self.log, debug=self.debug )
            if self.verbose: 
                self.log.add('adding hydrogen atoms to input structure\n')
    
            self.delphimodel = reducer.run()
        else:
            self.delphimodel = self.model.clone()
            
        self.delphimodel.xplor2amber()

        if not os.path.exists( self.f_radii ):
            self.__prepareRadii()
        
        if not os.path.exists( self.f_charges ):
            self.__prepareCharges( self.f_charges )
        
        self.delphimodel.writePdb( self.f_pdb )

    def cleanup( self ):
        """
        Tidy up the mess you created.
        """        
        if not self.debug:
            T.tryRemove( self.f_pdb )
            if not self.keep_map:
                T.tryRemove( self.f_map )

        Executor.cleanup( self )

    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.output is None or \
               not 'energy calculations done' in self.output

    def fail( self ):
        """
        Overrides Executor method. Called when execution fails.
        """
        s = 'Delphi failed. Please check the program output in the '+\
          'field `output` of this Delphi instance (e.g. `print x.output`)!\n'
        self.log.add( s )
        if self.output:
            s = 'The last message from DelPhi reads as follows:\n'
            s += '\n'.join( self.output.split('\n')[-3:] )
            self.log.add( s )
        else:
            self.log.add( 'There does not seem to be any DelPhi output.')

        raise DelphiError(s)

    def postProcess( self ):
        """
        Called directly after execution. Read delphi output.
        """
        try:
            f = open( self.f_out, 'r')
            self.output = f.read()
            f.close()
        except IOError:
            self.output = None
    

    def parseOutput( self ):
        """
        Assumes output file has been parsed into self.output
        """
        r = {}
        for pattern in [self.RE_E_COUL, self.RE_E_GRID, self.RE_E_RXN, 
                        self.RE_E_SELF, self.RE_E_RXNT, self.RE_SURFCH]:
            ex = re.compile( pattern )
            hit = ex.search( self.output )
            try:
                r.update( hit.groupdict() )
            except:
                if self.verbose:
                    self.log.writeln('Warning, no match for: ' + pattern)
        
        for k, v in r.items():
            r[k] = float( v )
            
        return r
    
    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parseOutput()
            


#############
##  TESTING        
#############
import biskit.test as BT
import tempfile

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE, BT.LONG, BT.FAILS ]
    MODEL= None

    def prepare( self ):
        self.fcrg = tempfile.mktemp( '.crg', 'delphicharges_' )
        self.fmap = tempfile.mktemp( '.phi','delphimap_' )

        
    def cleanUp( self ):
        if not self.DEBUG:
            T.tryRemove( self.fcrg )
            T.tryRemove( self.fmap )

    
    def test_delphi( self ):
        """Delphi test"""
        if self.local: print('Loading PDB...')

        self.m1 = self.MODEL or PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        Test.MODEL = self.m1
        self.m1.addChainFromSegid()

        if self.local: print('Starting Delphi')
        self.x = Delphi( self.m1, scale=1.2, debug=self.DEBUG,
                         verbose=self.local, f_map=self.fmap )

        if self.local:
            print('Running')

        self.r = self.x.run()

        if self.local:
            print("Result: ")
            print(self.r)
            
        expect_delphi_v5 = {'scharge': 1.427, 'egrid': 9091., 
                            'ecoul': -9870, 'eself': -20420, 'erxn': -666.7}

        expect_delphi_v6 = {'scharge': 1.426, 'egrid': 8765., 
                            'ecoul': -10335., 'eself': -20420., 'erxn': -664.} 
##                            'erxnt': -21083. }
        expect = expect_delphi_v6

        if self.local:
            print("verifying results... ")
            print("Note: numeric values can differ on different hardware.")
        
        for k, v in list(expect.items()):
            self.assertAlmostEqual( expect[k], self.r[k], -1 )
        
        self.assertTrue(os.path.exists( self.fmap ), 'Potential map not found' )
        self.assertTrue(os.path.getsize( self.fmap) > 1000, 'empty potential map')


    def test_delphiCharges2( self ):
        """
        PDB2DelphiCharges test
        """
        if self.local:
            T.errWrite( 'loading PDB...' )

        self.m1 = self.MODEL or PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        Test.MODEL = self.m1
        if self.local:
            T.errWriteln( 'Done.' )
        
        if self.local:
            T.errWrite( 'Adding hydrogens to model (reduce)...' )

        self.rmodel = Reduce( self.m1, verbose=self.local ).run()
        self.rmodel.xplor2amber()
        if self.local:
            T.errWriteln( 'Done.' )

        ac = AtomCharger()
        ac.charge(self.rmodel)
        self.rmodel.addChainFromSegid()
        
        self.dc = PDB2DelphiCharges( self.rmodel )
        self.dc.prepare()
        
        self.assertEqual( len(self.dc.resmap['LYS']), 2 )  # normal and N'
        self.assertEqual( len(self.dc.resmap['SER']), 2 )  # normal and C'
        
        if self.local:
            T.errWriteln( 'writing delphi charge file to %s' % self.fcrg )
        self.dc.tofile( self.fcrg )
        
        self.assertTrue( os.path.exists( self.fcrg ) )
        
        

if __name__ == '__main__':

    BT.localTest(debug=False)
    
