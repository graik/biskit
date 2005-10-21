##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
## 
## last $Author$
## last $Date$
## $Revision$

import tempfile, re
import os.path
import string
import Numeric

from Biskit import Executor, TemplateError
import Biskit.settings as S
import Biskit.tools as T

class SurfaceRacer_Error( Exception ):
    pass


class SurfaceRacer( Executor ):
    """
    Runs surface_racer_3 with given PDBModel, and returns a dictionary with
    accessible, molecular surface areas (MS and AS) and average curvature.
    Hydrogens should not be present in the pdb-file during the calculation.

    Options:
        probe    - float, probe radii
        vdw_set  - int, Van del Waals radii set
                      1 - Richards (1977)
                      2 - Chothia  (1976)
        mode     - int, calculation mode
                      1- Accessible surface area only
                      2- Accessible and molecular surface areas
                      3- Accessible, molecular surface areas and average curvature

    Example usage:
        x = SurfaceRacer( model, 1.4, verbose=1 )
        result = x.run()

    References:
       http://monte.biochem.wisc.edu/~tsodikov/surface.html
       Tsodikov, O. V., Record, M. T. Jr. and Sergeev, Y. V. (2002).
       A novel computer program for fast exact calculation of accessible
       and molecular surface areas and average surface curvature.
       J. Comput. Chem., 23, 600-609. 
    """
    
    inp = \
"""%(vdw_set)i
%(f_pdb_name)s
%(probe).2f
%(mode)i
         
        
"""
    
    def __init__( self, model, probe, vdw_set=1, mode=3, **kw ):
        """
        SurfaceRacer creates three output files:
          result.txt - contains breakdown of surface areas and is writen
                         to the directory where the program resides. This
                         file is discharded here.
          file *.txt - contains the accessible, molecular surface areas
                         and average curvature information paresd here.
                         The filename is that of the input pdb file but
                         with a .txt extension.
          stdout     - some general information about the calculation.
                         Redirected to /dev/null
                  
        model    - PDBModel
        probe    - float, probe radii
        vdw_set  - int, Van del Waals radii set
                      1 - Richards (1977)
                      2 - Chothia  (1976)
        mode     - int, calculation mode
                      1- Accessible surface area only
                      2- Accessible and molecular surface areas
                      3- Accessible, molecular surface areas and average curvature

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
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
                raise FastSurf_Error, 'Cannot find SurfaceRacer directory. Set your path in ~/.biskit/settings.dat as surfaceracer_bin'

        Executor.__init__( self, 'surfaceracer', template=self.inp,\
                           f_out='/dev/null', cwd=dir, **kw )

        self.model = model.clone( deepcopy=1 )
                    
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


    def prepare( self ):
        """Overrides Executor method."""
        self.__prepareModel( self.model, self.f_pdb ) 


    def __prepareModel( self, model, f_pdb_out ):
        """
        Surface curvature should only be calculated on heavy atoms.
        Delete Hydrogens, sequential numbering ...
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
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )
            T.tryRemove( self.f_out_name )
            

    def parse_result( self, output ):
        """
        Parse the SurfaceRacer output file which has the same mane as the input
        pdb, but with a txt extension. The output ends up un the same folser
        as the input. In addition a file called result.txt is created in the
        same directory as the binary.
        pdbFile - str, full path to input pdb-file
        """
        curv = [] ## average curvature
        ms   = [] ## molecular surface area
        as   = [] ## accessible surface area
        
        lines = open( self.f_out_name ).readlines()
    
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

        result = {'curvature':Numeric.array(curv),
                  'MS':Numeric.array(ms),
                  'AS':Numeric.array(as),
                  'surfaceRacerInfo':{'probe_radius':self.probe,
                                  'vdw_set':self.vdw_set}
                  }

        ## check curvature profile integrity
        result['curvature'] = \
             self.__checkProfileIntegrity( result['curvature'], 1.0, -1.0 )
        
        return result


    def __checkProfileIntegrity( self, profile, upperLimit=1.0, lowerLimit=-1.0):
        """
        In some cases SurfaceRacer generates incorrect curvature
        values for some atoms. This function sets values outside
        a given range to 0
        """
        mask = Numeric.greater( profile, upperLimit )
        mask += Numeric.less( profile, lowerLimit )

        for i in  Numeric.nonzero(mask):
            print 'WARNING! Profile value %.2f set to O\n'%profile[i]
            profile[i] = 0
            
        return profile


    def isFailed( self ):
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_result( self.output )
        


#######
## test
if __name__ == '__main__':

    from Biskit import PDBModel
    import Biskit.tools as T
    import glob
    import Biskit.mathUtils as MA

    print "Loading PDB..."

    f = glob.glob( T.testRoot()+'/rec_pc2_00/pdb/*_1_*pdb.gz' )[1]
    m = PDBModel(f)
    m = m.compress( m.maskProtein() )

    print "Starting SurfaceRAcer"
    x = SurfaceRacer( m, 1.4, vdw_set=2, debug=1, verbose=1 )

    print "Running"
    r = x.run()

    print "Result: ",
    
    c= r['curvature']
    ms= r['MS']
    print "weighted mean %.6f and standard deviation %.3f"%(MA.wMean(c,ms), MA.wSD(c,ms))
