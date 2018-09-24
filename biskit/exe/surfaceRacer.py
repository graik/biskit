## numpy-oldnumeric calls replaced by custom script; 09/06/2016
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
Calculates accessible, molecular surface areas and average curvature.
"""
## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import tempfile
import os.path
import string
import biskit.core.oldnumeric as N0

from biskit import EHandler
import biskit.settings as S
import biskit.tools as T

from .executor import Executor, TemplateError
from . import surfaceRacerTools as SRT

class SurfaceRacer_Error( Exception ):
    pass


class SurfaceRacer( Executor ):
    """
    Run SurfaceRacer
    ================

    Runs surface_racer_5 with given PDBModel, and returns a dictionary with
    accessible, molecular surface areas (MS and AS) and average curvature.
    Hydrogens should not be present in the pdb-file during the calculation.
    If a probe radius of 1.4 A and the Richards vdw radii set is used the
    relative exposure is also calculated.
    
    Note, surfaceRacer does not tolerate most non-protein or hydrogen atoms.
    Hydrogen andn solvent atoms are filtered out by default before calling
    the program. The mask used is accessible as SurfaceRacer.mask. This means,
    the arrays of values comming back from the calculation do not neccessarily
    have the same lenght as the PDBModel going in. 
    You will probably want to use PDBDope.addSurfaceRacer which is dealing 
    with this problem automatically. Better still, filter your model manually
    before giving it to surfaceRacer. E.g.:
    
        >>> model = model.compress( model.maskProtein() )
        
    ... should be the most robust approach.

    Options
    -------
        probe    - float, probe radii
        vdw_set  - int, Van der Waals radii set
                      1. Richards (1977)
                      2. Chothia  (1976)
        mode     - int, calculation mode
                      1. Accessible surface area only
                      2. Accessible and molecular surface areas
                      3. Accessible, molecular surface areas and average
                         curvature

    Example usage
    -------------
        >>> x = SurfaceRacer( model, 1.4, verbose=1 )
        >>> result = x.run()

    References
    ----------
       - U{http://monte.biochem.wisc.edu/~tsodikov/surface.html}
       - Tsodikov, O. V., Record, M. T. Jr. and Sergeev, Y. V. (2002).
         A novel computer program for fast exact calculation of accessible
         and molecular surface areas and average surface curvature.
         J. Comput. Chem., 23, 600-609. 
    """

    inp = \
"""
%(vdw_set)i

%(f_pdb_name)s
%(probe).2f
%(mode)i


"""

    def __init__( self, model, probe, vdw_set=1, mode=3, mask=None, **kw ):
        """
        SurfaceRacer creates three output files::
          result.txt - contains breakdown of surface areas and is writen
                         to the directory where the program resides. This
                         file is discarded here.
          <file>.txt - contains the accessible, molecular surface areas
                         and average curvature information parsed here.
                         The filename is that of the input pdb file but
                         with a .txt extension.
          <file>_residue.txt - new in version 5.0 and not used by this wrapper
          stdout     - some general information about the calculation.
                         Redirected to /dev/null

        @param model: model analyze
        @type  model: PDBModel
        @param probe: probe radii, Angstrom
        @type  probe: float
        @param vdw_set: Van del Waals radii set (default: 1)::
                          1 - Richards (1977)
                          2 - Chothia  (1976)
        @type  vdw_set: 1|2
        @param mode: calculation mode (default: 3)::
                      1- Accessible surface area only
                      2- Accessible and molecular surface areas
                      3- Accessible, molecular surface areas and
                         average curvature
        @type  mode: 1|2|3
        @param mask: optional atom mask to apply before calling surface racer
                     (default: heavy atoms AND NOT solvent)
        @type mask: [ bool ]

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

        Executor.__init__( self, 'surfaceracer', template=self.inp,\
                           **kw )

        self.model = model.clone()
        self.mask = mask if mask is not None else \
            model.maskHeavy() * N0.logical_not( model.maskSolvent())
        self.model = self.model.compress( self.mask )

        ## will be filled in by self.prepare() after the temp folder is ready
        self.f_pdb = None
        self.f_pdb_name = None
        self.f_out_name = None

        ## parameters that can be changed
        self.probe = probe
        self.vdw_set = vdw_set
        self.mode = mode

        ## random data dictionaries
        self.ranMS = SRT.ranMS
        self.ranAS = SRT.ranAS
        self.ranMS_Nter = SRT.ranMS_N
        self.ranAS_Nter = SRT.ranAS_N
        self.ranMS_Cter = SRT.ranMS_C
        self.ranAS_Cter = SRT.ranAS_C

        ## count failures
        self.i_failed = 0


    def prepare( self ):
        """
        Overrides Executor method.
        - link surfaceracer binary to temporary folder
        - put cleaned up pdb into same folder
        - adapt working directory and output names to same folder
        (changes self.cwd, f_pdb, f_pdb_names, f_out_name)
        """
        self.cwd = self.__prepareFolder()

        self.f_pdb = tempfile.mktemp( '_surfaceracer.pdb', self.cwd )
        self.f_pdb_name = os.path.split(self.f_pdb)[1]

        ## The SurfRace output file has the same name as the input
        ## pdb, but with a txt extension.
        self.f_out_name = self.f_pdb[:-3]+'txt'

        self.__prepareModel( self.model, self.f_pdb )


    def __prepareFolder( self ):
        """
        The surfrace binary, radii files and the PDB have to be in the
        current working directory. Otherwise we get a pretty silent
        segmentation fault.

        @return temp folder
        @rtype str
        """
        try:
            folder = T.tempDir()
            binlnk = os.path.join( folder, self.exe.name )

            if not os.path.exists( binlnk ):
                os.symlink( self.exe.bin, binlnk )

            radii_txt = T.dataRoot() + '/surface_racer_3/radii.txt'

            target = os.path.join(folder, 'radii.txt')
            if not os.path.exists( target ):
                os.symlink( radii_txt, target )

            self.exe.bin = binlnk

            return folder + '/'

        except OSError as error:
            raise SurfaceRacer_Error('Error preparing temporary folder for SurfaceRacer\n'+\
                  'Error: %r\n' % error +\
                  'folder: %r\n' % folder +\
                  'binary: %r\n' % self.exe.bin +\
                  'binary link: %r' % binlnk) 



    def __prepareModel( self, model, f_pdb_out ):
        """
        Prepare a model that SurfaceRacer likes.
         - Surface curvature should only be calculated on heavy atoms.
         - Delete Hydrogens, sequential numbering ...

        @param model: model 
        @type  model: PDBModel
        @param f_pdb_out: name of pdb file to write
        @type  f_pdb_out: str

        @raise SurfaceRacer_Error: if contains other than heavy atoms
        """
        ## Hydrogens should not be present in calculation of curvature
        model = model.compress( model.maskHeavy() )
        model.renumberResidues()
        if sum(model.maskHeavy()) == model.lenAtoms():
            model.writePdb( f_pdb_out, wrap=1, left=0 )
        else:
            raise SurfaceRacer_Error('The pdb file that was to be written as input for SurfaceRacer contains none heavy atoms.')


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb, verbose=self.verbose )
            T.tryRemove( self.f_out_name, verbose=self.verbose )
            T.tryRemove( os.path.join(self.cwd, 'result.txt'), 
                         verbose=self.verbose )
            T.tryRemove( self.f_pdb[:-4]+'_residue.txt', verbose=self.verbose)


    def parse_result( self):
        """
        Parse the SurfaceRacer output file which has the same nawe as the input
        pdb, but with a txt extension. The output ends up un the same folder
        as the input. In addition a file called result.txt is created in the
        same directory as the binary.

        @return: dictionary with curvature and surface data
        @rtype: dict
        """
        curv = [] ## average curvature
        ms   = [] ## molecular surface area
        asa  = [] ## accessible surface area

        try:
            out_file = open( self.f_out_name )
            lines = out_file.readlines()
            out_file.close()
        except:
            raise SurfaceRacer_Error('SurfaceRacer result file %s does not exist. You have probably encountered a very rare SurfaceRacer round off error that have caused the program to terminate. The simplest remedy to this problem is to increase the probe radii with a very small number, for example from %.3f to %.3f.'%(self.f_out_name, self.probe,self.probe+0.001  ))

        if  len(lines) == 0:
            raise SurfaceRacer_Error('SurfaceRacer result file %s empty'%self.f_out_name)

        ## don't parse cavity information, find first occurance or 'CAVITY'
        end = len(lines)
        for i in range( len(lines)-1, 0, -1 ):
            if lines[i][:6]=='CAVITY':
                end = i

        for i in range( end ):
            curv += [ float( str.strip( lines[i][-11:-1] ) ) ]
            ms   += [ float( str.strip( lines[i][-17:-11] ) ) ]
            asa  += [ float( str.strip( lines[i][-24:-17] ) ) ]

        result = {'curvature':N0.array(curv),
                  'MS':N0.array(ms),
                  'AS':N0.array(asa),
                  'surfaceRacerInfo':{'probe_radius':self.probe,
                                      'vdw_set':self.vdw_set}
                  }

        ## check curvature profile integrity
        result['curvature'] = \
              self.__checkProfileIntegrity( result['curvature'], 1.0, -1.0 )

        return result


    def __checkProfileIntegrity( self, profile, upperLimit=1.0,
                                 lowerLimit=-1.0):
        """
        In some cases SurfaceRacer generates incorrect curvature
        values for some atoms. This function sets values outside
        a given range to 0

        @param profile: profile name
        @type  profile: str
        @param upperLimit: upper limit for a valid value (default: 1.0)
        @type  upperLimit: float
        @param lowerLimit: lower limit for a valid value (default: -1.0)
        @type  lowerLimit: float

        @return: profile with inspected values
        @rtype: [float]
        """
        mask = N0.greater( profile, upperLimit )
        mask += N0.less( profile, lowerLimit )

        for i in  N0.nonzero(mask):
            print('WARNING! Profile value %.2f set to O\n'%profile[i])
            profile[i] = 0

        return profile


    def __relExposure( self, key='AS', clip=1 ):
        """
        Calculate how exposed an atom is relative to the same
        atom in a GLY-XXX-GLY tripeptide, an approximation of
        the unfolded state. See L{Biskit.surfaceRacerTools.relExposure()}

        @param key: caclulate relative Molecular Surface or
                    Acessible Surface (default: AS)
        @type  key: MS|AS
        @param clip: clip MS ans AS values above 100% to 100% (default: 1)
        @type  clip: 1|0

        @return: relAS - list of relative solvent accessible surfaces OR
                 relMS - list of relative molecular surface exposure
        @rtype: [float]
        """
        if not key=='MS' and not key=='AS':
            raise SurfaceRacer_Error('Incorrect key for relative exposiure: %s '%key)

        self.result['rel'+key ] = SRT.relExposure( self.model,
                                                   self.result[key],
                                                   key )

    def isFailed( self ):
        """
        Overrides Executor method
        """
        if not os.path.exists(self.f_out_name):
            T.errWrite( '\nSurfaceRacer result file %s does not exist. You have probably encountered a very rare SurfaceRacer round off error that have caused the program to terminate. Will now try to recalculate the surface with a slightly increased surface probe radii: increasing radii from %.3f to %.3f.\n'%(self.f_out_name, self.probe,self.probe+0.001))
            return 1

        return self.error is None 


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )

        self.result = self.parse_result()

        ## if probe radius other than 1.4 A the relative surface exposure
        ## cannot be calculated, but allow this check to be a little flexible
        ## if we ate forced to slightly increase the radii to excape round off
        ## SurfaceRacer errors
        try:
            if round(self.probe, 1) == 1.4 and self.vdw_set == 1:
                self.__relExposure('MS')
                self.__relExposure('AS')
            else:
                EHandler.warning("No relative accessabilities calculated "+\
                                 "when using a prob radius other than 1.4 A"+\
                                 " or not using the Richards vdw radii set.")
        except KeyError as what:
            EHandler.warning("Missing standard accessibilities for some "+\
                             "atoms. No relative accesibilities calculated.")
            if 'relMS' in self.result: del self.result['relMS']
            if 'relAS' in self.result: del self.result['relAS']


    def fail( self ):
        """
        Called if external program failed, Overrides Executor method.

        In some very rare cases SurfaceRacer round off error cause the program
        to terminate. The simplest remedy to this problem is to increase the
        probe radii with a very small number and rerun the calculation.
        """
        self.i_failed += 1

        if self.i_failed < 2:
            self.probe = self.probe + 0.001
            self.run()

        Executor.fail( self )

#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.EXE ]

    def test_SurfaceRacer(self):
        """SurfaceRacer test"""

        from biskit import PDBModel
        import biskit.mathUtils as MA
        import numpy as N

        if self.local: print('Loading PDB...')
        f = T.testRoot('/lig/1A19.pdb')
        m = PDBModel(f)
        m = m.compress( m.maskProtein() )

        if self.local: print('Starting SurfaceRacer')

        self.x = SurfaceRacer( m, 1.4, vdw_set=1, debug=self.DEBUG, verbose=0 )

        if self.local:
            print('Running ...')

        self.r = self.x.run()

        c= self.r['curvature']
        ms= self.r['MS']

        if self.local:
            print("Curvature: weighted mean %.6f and standard deviation %.3f"\
                  %(MA.wMean(c,ms), MA.wSD(c,ms)))

            print('Relative MS of atoms 10 to 20:',self.r['relMS'][10:20])

            print('Relative AS of atoms 10 to 20:',self.r['relAS'][10:20])

        self.e = ( N0.sum(self.r['relMS'][10:20]), N0.sum(self.r['relAS'][10:20]),
              N0.sum(self.r['curvature'][10:20]) )

        delta = N.sum(N.abs(N.array(self.e) - N.array(self.EXPECT)))
        self.assertAlmostEqual( delta, 0.0, 5 )


    EXPECT = (570.29631242113669, 356.72813304721222, 0.80000000000000004)


if __name__ == '__main__':

    BT.localTest(debug=False)
