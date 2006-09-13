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
## $Revision$

"""
Calculates accessible, molecular surface areas and average curvature.
"""

import tempfile, re
import os.path
import string
import Numeric as N

from Biskit import Executor, TemplateError
import Biskit.settings as S
import Biskit.tools as T

import Biskit.surfaceRacerTools as SRT

class SurfaceRacer_Error( Exception ):
    pass


class SurfaceRacer( Executor ):
    """
    Run SurfaceRacer
    ================
    
    Runs surface_racer_3 with given PDBModel, and returns a dictionary with
    accessible, molecular surface areas (MS and AS) and average curvature.
    Hydrogens should not be present in the pdb-file during the calculation.
    If a probe radius of 1.4 A and the Richards vdw radii set is used the
    relative exposure is also calculated.

    Options
    -------
        probe    - float, probe radii
        vdw_set  - int, Van del Waals radii set
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

    def __init__( self, model, probe, vdw_set=1, mode=3, **kw ):
        """
        SurfaceRacer creates three output files::
          result.txt - contains breakdown of surface areas and is writen
                         to the directory where the program resides. This
                         file is discharded here.
          file *.txt - contains the accessible, molecular surface areas
                         and average curvature information paresd here.
                         The filename is that of the input pdb file but
                         with a .txt extension.
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

        ## FastSurf have to be run i the local directory where
        ##   the application resides. Look for the application folder.
        try:
            if os.path.exists( S.surfaceracer_bin ):
                dir = os.path.split(S.surfaceracer_bin)[0] +'/'
        except:
            if os.path.exists( T.projectRoot() +'/external/surface_racer_3'):
                dir =  T.projectRoot() +'/external/surface_racer_3/'
            else:
                raise SurfaceRacer_Error, 'Cannot find SurfaceRacer directory. Set your path in ~/.biskit/settings.dat as surfaceracer_bin'

        Executor.__init__( self, 'surfaceracer', template=self.inp,\
                           cwd=dir, **kw )

        self.model = model.clone( deepcopy=1 )
        self.model = self.model.compress( self.model.maskHeavy() )

        ## temporary pdb-file
        self.f_pdb = tempfile.mktemp( '_surfaceracer.pdb', self.cwd )
        self.f_pdb_name = os.path.split(self.f_pdb)[1]

        ## The SurfRace output file has the same name as the input
        ## pdb, but with a txt extension.
        self.f_out_name = self.f_pdb[:-3]+'txt'

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

    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.__prepareModel( self.model, self.f_pdb ) 


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
            raise SurfaceRacer_Error, \
                  'The pdb file that was to be written as input for SurfaceRacer contains none heavy atoms.'


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )
            T.tryRemove( self.f_out_name )


    def parse_result( self, output ):
        """
        Parse the SurfaceRacer output file which has the same nawe as the input
        pdb, but with a txt extension. The output ends up un the same folder
        as the input. In addition a file called result.txt is created in the
        same directory as the binary.
        
        @param output: full path to input pdb-file
        @type  output: str

        @return: dictionary with curvature and surface data
        @rtype: dict
        """
        curv = [] ## average curvature
        ms   = [] ## molecular surface area
        as   = [] ## accessible surface area

        try:
            out_file = open( self.f_out_name )
            lines = out_file.readlines()
            out_file.close()
        except:
            raise SurfaceRacer_Error,\
                  'SurfaceRacer result file %s does not exist. You have probably encountered a very rare SurfaceRacer round off error that have caused the program to terminate. The simplest remedy to this problem is to increase the probe radii with a very small number, for example from %.3f to %.3f.'%(self.f_out_name, self.probe,self.probe+0.001  )
        
        if  len(lines) == 0:
            raise SurfaceRacer_Error,\
                  'SurfaceRacer result file %s empty'%self.f_out_name

        ## don't parse cavity information, find first occurance or 'CAVITY'
        end = len(lines)
        for i in range( len(lines)-1, 0, -1 ):
            if lines[i][:6]=='CAVITY':
                end = i

        for i in range( end ):
            curv += [ float( string.strip( lines[i][-11:-1] ) ) ]
            ms   += [ float( string.strip( lines[i][-17:-11] ) ) ]
            as   += [ float( string.strip( lines[i][-24:-17] ) ) ]

        result = {'curvature':N.array(curv),
                  'MS':N.array(ms),
                  'AS':N.array(as),
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
        mask = N.greater( profile, upperLimit )
        mask += N.less( profile, lowerLimit )

        for i in  N.nonzero(mask):
            print 'WARNING! Profile value %.2f set to O\n'%profile[i]
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
            raise SurfaceRacer_Error,\
                  'Incorrect key for relative exposiure: %s '%key

        self.result['rel'+key ] = SRT.relExposure( self.model,
                                                   self.result[key],
                                                   key )

    def isFailed( self ):
        """
        Overrides Executor method
        """
        if not os.path.exists(self.f_out_name):
            T.flushPrint( '\nSurfaceRacer result file %s does not exist. You have probably encountered a very rare SurfaceRacer round off error that have caused the program to terminate. Will now try to recalculate the surface with a slightly increased surface probe radii: increasing radii from %.3f to %.3f.\n'%(self.f_out_name, self.probe,self.probe+0.001))
            return 1
        
        return self.error is None 


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )

        self.result = self.parse_result( self.output )

        ## if probe radius other than 1.4 A the relative surface exposure
        ## cannot be calculated, but allow this check to be a little flexible
        ## if we ate forced to slightly increase the radii to excape round off
        ## SurfaceRacer errors
        if round(self.probe, 1) == 1.4 and self.vdw_set == 1:
            self.__relExposure('MS')
            self.__relExposure('AS')
        else:
            T.flushPrint("\nNo relative accessabilities calculated when using a prob radius other than 1.4 A or not using the Richards vdw radii set.")


    def failed( self ):
        """
        Called if external program failed, Overrides Executor method.

        In some very rare cases SurfaceRacer round off error cause the program
        to terminate. The simplest remedy to this problem is to increase the
        probe radii with a very small number and rerun the calculation.
        """
        self.probe = self.probe + 0.001
        self.run()
        
    
#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: sum of relMS, relAS and curvature for atoms 10 to 20
        @rtype:  float, float, float
        """
        from Biskit import PDBModel
        import Biskit.mathUtils as MA
        
        if local: print 'Loading PDB...'
        f = T.testRoot()+'/lig/1A19.pdb'
        m = PDBModel(f)
        m = m.compress( m.maskProtein() )

        if local: print 'Starting SurfaceRacer'
        x = SurfaceRacer( m, 1.4, vdw_set=1, debug=0, verbose=0 )

        if local: print 'Running'
        r = x.run()

        c= r['curvature']
        ms= r['MS']
    
        if local:
            print "Curvature: weighted mean %.6f and standard deviation %.3f"\
              %(MA.wMean(c,ms), MA.wSD(c,ms))

            print 'Relative MS of atoms 10 to 20 atoms:', r['relMS'][10:20]

            print 'Relative AS of atoms 10 to 20 atoms:', r['relAS'][10:20]
            
            globals().update( locals() )
        

        return N.sum(r['relMS'][10:20]), N.sum(r['relAS'][10:20]), N.sum(r['curvature'][10:20])


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: sum of relMS, relAS and curvature for atoms 10 to 20
        @rtype:  float, float, float
        """
        
        return 570.47829086283639, 356.81939295543083, 0.80000000000000004

    
if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()



