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
Delphi-based protocol for electrostatic component of free energy of binding.
"""

import numpy as N
import copy

from biskit import PDBModel
from biskit import StdLog
from biskit import AtomCharger
from biskit.exe import Delphi, DelphiError, Reduce
from biskit.dock import Complex

import biskit.tools as T
import biskit.mathUtils as U

class DelphiBindingEnergy( object ):
    """
    Determine the electrostatic component of the free energy of binding using
    several rounds of Delphi calculations. DelphiBindingEnergy accepts a 
    binary complex (L{Biskit.Dock.Complex}) as input, performs several house
    keeping tasks (optional capping of free terminals, h-bond optimization and
    protonation with the reduce program, assignment of Amber partial atomic 
    charges) and then calls Delphi six times:

    1 - 3: one run each for complex, receptor, and ligand with zero ionic 
           strength
    4 - 6: one run each for complex, receptor, and ligand with custom salt
           concentration (default: physiological 0.15 M).
           
    The grid position and dimensions are determined once for the complex and
    then kept constant between all runs so that receptor, ligand and complex
    calculations can indeed be compared to each other.
           
    The free energy of binding is calculated as described in the Delphi
    documentation -- see: delphi2004/examples/binding.html -- according to the
    energy partioning scheme:
    
    dG = dG(coulomb) + ddG(solvation) + ddG(ions)
    Please have a look at the source code of DelphiBindingEnergy.bindingEnergy() 
    for a detailed breakup of these terms.
    
    Use:
    ====
    
    >>> com = Biskit.Complex( 'myrec.pdb', 'mylig.pdb' )
    >>> dg = DelphiBindingEnergy( com, verbose=True )
    >>> r = dg.run()

    Variable r will then receive a dictionary with the free energy of binding
    and its three components:
    
         {'dG_kt': free E. in units of kT, 
          'dG_kcal': free E. in units of kcal / mol, 
          'dG_kJ': free E. in units of kJ / mol,
          'coul': coulomb contribution in kT , 
          'solv': solvation contribution in kT, 
          'ions': salt or ionic contribution in kT }
          
    The modified complex used for the calculation (with hydrogens added and
    many other small changes) can be recovered from the DelphiBindingEnergy
    instance:
    
    >>> modified_com = dg.delphicom
    
    The run() method assigns the result dictionary to the info record
    of both the original and the modified complex. That means:
    
    >>> r1 = com.info['dG_delphi']
    >>> r2 = modified_com.info['dG_delphi']
    >>> r1 == r2
    True
    
    The result of the original DelPhi calculations (with and without salt for
    receptor, ligand and complex) are assigned to the info dictionaries of the
    complex' receptor and ligand model as well as to the complex itself:
    
    com.info['delphi_0.15salt'] ...delphi run with 0.15M salt on complex
    com.info['delphi_0salt']    ...delphi run without salt on complex

    com.lig().info['delphi_0.15salt'] ...delphi run with 0.15M salt on ligand
    com.lig().info['delphi_0salt']    ...delphi run without salt on ligand
    com.rec().info['delphi_0.15salt'] ...delphi run with 0.15M salt on receptor
    com.rec().info['delphi_0salt']    ...delphi run without salt on receptor
    
    From there the individual values for solvation, ionic and couloumb
    contributions can be recovered. See L{Biskit.Delphi} for a description of
    this result dictionary.
     
    Customization:
    ==============
    
    The most important Delphi parameters (dielectric, salt concentration, 
    grid scales, ion and probe radius) can be adjusted by passing parameters
    to the constructor (see documentation of __init__). The default parameters
    should be reasonable. By default,
    we create a grid that covers every linear dimension to at least 60% 
    (perfil=60) and has a density of 2.3 points per Angstroem (scale=2.3). 
    Such high densities come at much larger computational cost. It is
    recommended to test different values and average results.
    
    Note: Any parameters that are not recognized by 
    DelphiBindingEnergy() will be passed on to the Biskit.Delphi instance
    of each Delphi run and, from there, passed on to Biskit.Executor.
    
    The internal pipeline of DelphiBindingEnergy consists of:
    
    * adding hydrogens and capping of protein chain breaks and terminals with
      Biskit.Reduce( autocap=True ). The capping is performed by 
      L{Biskit.PDBCleaner} called from within Reduce.
    
    * mapping Amber partial atomic charges into the structure with 
      Biskit.AtomCharger()
      
    * setting up the various delphi runs with L{Biskit.Delphi}.
    
    A point worth looking at is the automatic capping of protein chain
    breaks and premature terminals with ACE and NME residues. This should
    generally be a good idea but the automatic discovery of premature C-termini
    or N-termini is guess work at best. See L{Biskit.PDBCleaner} for a more
    detailed discussion. You can override the default behaviour by setting
    autocap=False (no capping at all, this is now the default) and you can 
    then provide a complex structure that is already pre-treated by 
    L{Biskit.Reduce}.
    For example:
    
    >>> m = PDBModel('mycomplex.pdb')
    >>> m = Reduce( m, capN=[0], capC=[2] )
    >>> com = Biskit.Complex( m.takeChains([0]), m.takeChains([1,2]) )
    >>> dg = DelphiBindingEnergy( com, protonate=False )
    
    In this case, the original structure would receive a ACE N-terminal
    capping of the first chain and the second chain would receive a C-terminal
    NME capping residue. The first chain is considered receptor and the second
    and third chain are considered ligand.
    
    @see: L{Biskit.PDBCleaner}
    @see: L{Biskit.Reduce}
    @see: L{Biskit.AtomCharger}
    @see: L{Biskit.Delphi}
    """
    
    def __init__(self, com, 
                 protonate=True, addcharge=True,
                 indi=4.0, exdi=80.0, salt=0.15, ionrad=2, prbrad=1.4, 
                 bndcon=4, scale=2.3, perfil=60, 
                 template=None, topologies=None,
                 f_charges=None,
                 verbose=True, debug=False, log=None, tempdir=None, cwd=None,
                 **kw ):
        """
        @param com: complex to analyze
        @type  com: Biskit.Complex
        @param protonate: (re-)build hydrogen atoms with reduce program (True)
                          see L{Biskit.Reduce}
        @type  protonate: bool
        @param addcharge: Assign partial charges from Amber topologies (True)
        @type  addcharge: bool
        @param indi: interior dilectric (4.0)
        @param exdi: exterior dielectric (80.0)
        @param salt: salt conc. in M (0.15)
        @param ionrad: ion radius (2)
        @param prbrad: probe radius (1.4) 
        @param bndcon: boundary condition (4, delphi default is 2)
        @param scale:  grid spacing (2.3)
        @param perfil: grid fill factor in % (for automatic grid, 60) 
        
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

        @param kw: additional key=value parameters for Delphi or Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        self.com = com
        self.delphicom = None
        
        self.protonate = protonate
        self.addcharge = addcharge
        ## DELPHI run parameters
        self.indi=indi  # interior dilectric(4.0)
        self.exdi=exdi  # exterior dielectric(80.0)
        self.salt=salt  # salt conc. in M (0.15)
        self.ionrad=ionrad # ion radius (2)
        self.prbrad=prbrad # probe radius (1.4) 
        self.bndcon=bndcon # boundary condition (4, delphi default is 2)
        
        ## DELPHI parameters for custom grid
        self.scale=scale   # grid spacing (1.2 / A)
        self.perfil=perfil # grid fill factor in % (for automatic grid, 60)
        
        ## DELPHI parameter file and Amber residue definitions or charge file
        self.template = template
        self.topologies = topologies
        self.f_charges = f_charges
        
        ## pump everything else into name space, too
        self.__dict__.update( kw )
        
        ## prepare globally valid grid
        self.grid = None
        
        self.verbose = verbose
        self.log = log or StdLog()
        self.debug = debug
        self.tempdir = tempdir
        self.cwd = cwd
        
        self.ezero = self.esalt = None # raw results assigned by run()

    def prepare( self ):
        """
        """
        mrec = self.com.rec()
        mlig = self.com.lig()
        
        m = mrec.concat( mlig )

        if self.protonate:
            if self.verbose:
                self.log.add( '\nRe-building hydrogen atoms...' )
            
            tempdir = self.tempdir
            if tempdir:
                tempdir += '/0_reduce'
                
            r = Reduce( m, tempdir=tempdir, log=self.log,
                        autocap=False,
                        debug=self.debug, verbose=self.verbose )
            m = r.run()
            
        if self.addcharge:
            if self.verbose:
                self.log.add( '\nAssigning charges from Amber topologies...')
                    
            ac = AtomCharger( log=self.log, verbose=self.verbose )
            ac.charge( m )
            
        mrec = m.takeChains( range( mrec.lenChains() ) )
        mlig = m.takeChains( range( mrec.lenChains(), m.lenChains() ) )
            
        self.delphicom = Complex( mrec, mlig )


    def setupDelphi( self, model, grid={}, **kw ):
        """
        """
        params = copy.copy( self.__dict__ )
        params.update( kw )
        params['protonate'] = False   # run reduce once for complex only
        params['addcharge'] = False   # run AtomCharger once for complex only
        
        d = Delphi( model, **params )
        d.setGrid( **grid )
        
        return d
    
    def processThreesome( self, **kw ):
        """
        Calculate Delphi energies for rec, lig and complex
        @return: result dictionaries for com, rec, lig
        @rtype: ( dict, dict, dict )
        """
        dcom = self.setupDelphi( self.delphicom.model(), **kw )
        grid = dcom.getGrid()
        
        drec = self.setupDelphi( self.delphicom.rec(), grid=grid, **kw )
        dlig = self.setupDelphi( self.delphicom.lig(), grid=grid, **kw )
        
        if self.verbose:
            self.log.add('\nrunning Delphi for complex...')
        r_com = dcom.run()

        if self.verbose:
            self.log.add('\nrunning Delphi for receptor...')
        r_rec = drec.run()
        
        if self.verbose:
            self.log.add('\nrunning Delphi for ligand...')
        r_lig = dlig.run()
        
        return r_com, r_rec, r_lig
    
    def processSixsome( self, **kw ):
        """
        Calculate Delphi energies for rec, lig and complex with and without
        salt.
        @return: com, rec, lig delphi calculations with and without ions
        @rtype: ( {}, {} )
        """
        ## with ions
        if self.verbose:
            self.log.add('\nDelphi calculations with %1.3f M salt' % self.salt)
            self.log.add('=======================================')
        ri_com, ri_rec, ri_lig = self.processThreesome( salt=self.salt, **kw )

        ## w/o ions
        if self.verbose:
            self.log.add('\nDelphi calculations with zero ionic strenght')
            self.log.add('==============================================')
        r_com, r_rec, r_lig = self.processThreesome( salt=0.0, **kw )
        
        rsalt = { 'com':ri_com, 'rec':ri_rec, 'lig':ri_lig }
        rzero = { 'com':r_com, 'rec':r_rec, 'lig':r_lig }
        
        if self.verbose:
            self.log.add('\nRaw results')
            self.log.add('============')
            self.log.add('zero ionic strength: \n' + str(rzero) )
            self.log.add('\nwith salt: \n' + str(rsalt) )
    
        return rzero, rsalt
    
    def bindingEnergy( self, ezero, esalt ):
        """
        Calculate electrostatic component of the free energy of binding 
        according to energy partitioning formula.
        @param ezero: delphi energies calculated at zero ionic strength
        @type  ezero: {'com':{}, 'rec':{}, 'lig':{} }
        @param esalt: delphi energies calculated with ions (see salt=)
        @type  esalt: {'com':{}, 'rec':{}, 'lig':{} }
        @return: {'dG_kt': free E. in units of kT, 
                  'dG_kcal': free E. in units of kcal / mol, 
                  'dG_kJ': free E. in units of kJ / mol,
                  'coul': coulomb contribution in kT , 
                  'solv': solvation contribution in kT, 
                  'ions': salt or ionic contribution in kT }
        @rtype: dict { str : float }
        """
        delta_0 = {}
        delta_i = {}
        for k, v in ezero['com'].items():
            delta_0[ k ] = ezero['com'][k] - ezero['rec'][k] - ezero['lig'][k]
            delta_i[ k ] = esalt['com'][k] - esalt['rec'][k] - esalt['lig'][k]
        
        dG_coul = delta_0['ecoul']
        ddG_rxn = delta_0['erxn']
        ddG_ions= delta_i['egrid'] - delta_0['egrid']
        
        dG = dG_coul + ddG_rxn + ddG_ions  # in units of kT
        dG_kcal = round( dG * 0.593, 3)  # kcal / mol
        dG_kJ = round( dG * 2.479, 3)    # kJ / mol
        
        r = {'dG_kt':dG, 'dG_kcal':dG_kcal, 'dG_kJ':dG_kJ,
             'coul':dG_coul, 'solv':ddG_rxn, 'ions':ddG_ions}
        return r
    
    def run( self ):
        """
        Prepare structures and perform Delphi calculations for complex, 
        receptor and ligand, each with and without ions (salt). Calculate
        free energy of binding according to energy partitioning.
        Detailed delphi results are saved into object fields 'ezero' and
        'esalt' for calculations without and with ions, respectively.
        
        @return: {'dG_kt': free E. in units of kT, 
                  'dG_kcal': free E. in units of kcal / mol, 
                  'dG_kJ': free E. in units of kJ / mol,
                  'coul': coulomb contribution in kT , 
                  'solv': solvation contribution in kT, 
                  'ions': salt or ionic contribution in kT }
        @rtype: dict { str : float }
        """
        self.prepare()
        self.ezero, self.esalt = self.processSixsome()
        self.result = self.bindingEnergy( self.ezero, self.esalt )
        
        self.delphicom.info['dG_delphi'] = self.result
        self.com.info['dG_delphi'] = self.result
        
        for com in [self.delphicom, self.com]:
            key = 'delphi_%4.2fsalt' % self.salt
            com.rec_model.info[key] = self.esalt['rec']
            com.lig_model.info[key] = self.esalt['lig']
            com.lig_model.lig_transformed = None
            key = 'delphi_0salt'
            com.rec_model.info[key] = self.ezero['rec']
            com.lig_model.info[key] = self.ezero['lig']
            com.lig_model.lig_transformed = None  # reset Complex.lig() cache
            com.info[key] = self.ezero['com']
        
        return self.result

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE, BT.LONG ]
    
    def test_bindingE( self ):
        """bindingEnergyDelphi test (Barnase:Barstar)"""
        self.com = T.load( T.testRoot() + '/com/ref.complex' )
        self.dG = DelphiBindingEnergy( self.com, log=self.log, scale=1.2,
                                       verbose=self.local, debug=self.debug )
        self.r = self.dG.run()

##        self.assertAlmostEqual( self.r['dG_kt'], 21., 0 )
        self.assertTrue( abs(self.r['dG_kt'] - 24.6) < 4 )

        if self.local:
            self.log.add(
                '\nFinal result: dG = %3.2f kcal/mol'%self.r['dG_kcal'])
            
    def test_errorcase1( self ):
        """bindinEnergyDelphi test (error case 01)"""
        self.m = PDBModel( T.testRoot() + '/delphi/case01.pdb' )
        rec = self.m.takeChains( [0,1] )
        lig = self.m.takeChains( [2] )
        self.com = Complex( rec, lig )
        
        self.dG = DelphiBindingEnergy( self.com, log = self.log, scale=0.5,
                                       verbose=self.local, debug=self.debug )
        self.r = self.dG.run()

        if self.local:
            self.log.add(
                '\nFinal result: dG = %3.2f kcal/mol'%self.r['dG_kcal'])
        

if __name__ == '__main__':

    BT.localTest(debug=False)
