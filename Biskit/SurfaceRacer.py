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

import tempfile, re
import os.path
import string
import Numeric

from Biskit import Executor, TemplateError
import Biskit.settings as S
import Biskit.tools as T


## dictionary with average accessabilities in 500 random
## GLY-XX-GLY peptides calculated with SurfaceRacer using
## a probe radius of 1.4 A

ranMS = {'CYS': {'C'  :  1.969, 'CB' : 21.339, 'CA' :  8.077,
                 'O'  : 14.038, 'N'  :  5.713, 'SG' : 34.959},
         'GLN': {'C'  :  1.736, 'CB' : 18.742, 'CA' :  6.384,
                 'CG' : 20.672, 'O'  : 14.136, 'CD' :  2.825,
                 'N'  :  7.409, 'NE2': 23.326, 'OE1': 14.894},
         'ILE': {'C'  :  1.151, 'CB' :  7.933, 'CA' :  4.369,
                 'O'  : 14.185,  'N' :  5.428, 'CD1': 31.332,
                 'CG1': 17.690, 'CG2': 28.741},
         'SER': {'C'  :  2.336, 'OG' : 17.413, 'CB' : 26.667,
                 'CA' :  8.803, 'O'  : 13.894, 'N'  :  6.849},
         'VAL': {'C'  :  1.332, 'CB' :  7.877, 'CA' :  6.509,
                 'O'  : 13.741, 'N'  :  4.786, 'CG1': 31.049,
                 'CG2': 30.663},
         'LYS': {'C'  :  1.722, 'CB' : 18.077, 'CA' :  6.360,
                 'CG' : 16.832, 'CE' : 18.940, 'CD' : 16.329,
                 'NZ' : 32.635, 'O'  : 14.175, 'N'  :  7.007},
         'ASN': {'C'  :  2.050, 'CB' : 13.362, 'CA' : 15.158,
                 'CG' : 20.928, 'O'  :  6.352, 'N'  :  1.645,
                 'OD1': 14.146, 'ND2':  6.379 },
         'PRO': {'C'  : 15.067, 'CB' : 14.374, 'CA' :  5.090,
                 'CG' :  4.201, 'O'  : 15.123, 'CD' : 15.837,
                 'N'  :  1.779},
         'GLY': {'CA' : 20.259, 'C'  :  5.731, 'O'  : 21.192,
                 'N'  : 15.386},
         'THR': {'C'  : 19.678, 'CB' :  0.712, 'CA' :  1.439,
                 'OG1': 11.571, 'O'  :  8.170, 'N'  : 15.934,
                 'CG2': 14.310},
         'PHE': {'C'  :  6.163, 'CE1': 32.630, 'CB' :  1.976,
                 'CA' : 12.422, 'CG' : 21.216, 'O'  :  7.309,
                 'N'  :  1.825, 'CZ' : 14.120, 'CD1':  6.826,
                 'CD2': 15.219, 'CE2': 12.588},
         'ALA': {'CB' : 15.389, 'CA' : 15.277, 'C'  : 34.475,
                 'O'  :  9.672, 'N'  :  1.697},
         'HIS': {'C'  : 14.138, 'CE1':  7.174, 'CB' : 26.080,
                 'CA' :  3.136, 'CG' : 13.978, 'O'  :  8.065,
                 'N'  :  1.666, 'CD2': 14.372, 'ND1': 21.724,
                 'NE2':  8.415},
         'MET': {'C'  :  1.564, 'CB' : 14.093, 'CA' :  7.444,
                 'CG' : 15.702, 'CE' : 10.286, 'N'  : 15.480,
                 'O'  :  2.250, 'SD' : 19.502},
         'ASP': {'C'  :  7.110, 'CB' : 22.634, 'CA' : 14.257,
                 'CG' :  3.331, 'O'  : 14.934, 'N'  :  5.535,
                 'OD1': 14.057, 'OD2':  1.949},
         'LEU': {'C'  : 16.627, 'CB' :  5.280, 'CA' :  8.389,
                 'CG' : 14.299, 'O'  :  5.688, 'N'  : 30.423,
                 'CD1': 29.766, 'CD2':  1.777},
         'ARG': {'C'  : 18.112, 'CB' :  5.241, 'CA' : 16.730,
                 'CG' : 16.616, 'NE' : 13.828, 'O'  : 15.959,
                 'CD' :  1.912, 'CZ' : 21.727, 'NH1': 21.265,
                 'NH2':  7.178, 'N'  :  1.769},
         'TRP': {'C'  : 22.097, 'CZ2':  9.303, 'CB' :  4.131,
                 'CA' : 14.342, 'CG' :  5.499, 'O'  : 16.075,
                 'N'  : 15.575, 'CH2':  1.621, 'CE3': 22.125,
                 'CE2':  8.268, 'CD2':  2.468, 'CZ3': 13.935,
                 'NE1':  7.413, 'CD1': 15.173},
         'GLU': {'C'  : 22.925, 'CB' :  2.084, 'CA' : 15.389,
                 'CG' : 20.656, 'O'  : 21.058, 'CD' :  7.101,
                 'OE2':  2.090, 'N'  : 14.298, 'OE1':  6.603},
         'TYR': {'C'  :  5.766, 'CD2': 12.756, 'OH' : 13.124,
                 'CB' : 15.263, 'CA' :  1.525, 'CG' : 18.792,
                 'O'  :  6.204, 'N'  : 18.807, 'CZ' : 35.593,
                 'CD1':  7.103, 'CE1': 14.230, 'CE2': 14.968}}

ranAS = {'CYS': {'C'  :  0.833, 'CB' : 34.358, 'CA' : 7.768,
                 'O'  : 29.027, 'N'  :  5.159, 'SG' : 76.311},
         'GLN': {'C'  :  0.694, 'CB' : 28.846, 'CA' :  3.896,
                 'CG' : 32.484, 'O'  : 28.716, 'CD' :  2.467,
                 'N'  :  7.736, 'NE2': 54.470, 'OE1': 32.637},
         'ILE': {'C'  :  0.336, 'CB' : 10.203, 'CA' :  2.117,
                 'O'  : 26.670, 'N'  :  4.723, 'CD1': 63.709,
                 'CG1': 32.543, 'CG2': 50.456},
         'SER': {'C'  :  1.003, 'OG' : 47.009, 'CB' : 44.666,
                 'CA' :  8.688, 'O'  : 29.511, 'N'  :  6.051},
         'VAL': {'C'  :  0.491, 'CB' :  9.140, 'CA' :  5.240,
                 'O'  : 23.612, 'N'  :  2.422, 'CG1': 65.475,
                 'CG2': 59.022},
         'LYS': {'C'  :  0.713, 'CB' : 27.111, 'CA' :  4.234,
                 'CG' : 24.319, 'CE' : 35.072, 'CD' : 23.182,
                 'NZ' : 76.704, 'O'  : 27.814, 'N'  :  7.147},
         'ASN': {'C'  :  0.671, 'CB' : 33.549, 'CA' :  6.662,
                 'CG' :  2.403, 'O'  : 27.881, 'N'  :  7.635,
                 'OD1': 36.200, 'ND2': 50.093},
         'PRO': {'C'  :  0.696, 'CB' : 34.966, 'CA' :  4.790,
                 'CG' : 46.846, 'O'  : 34.264, 'CD' : 32.028,
                 'N'  :  0.314},
         'GLY': {'CA' : 43.741, 'C'  :  2.298, 'O'  : 31.107,
                 'N'  : 12.702},
         'THR': {'C'  :  0.709, 'CB' : 15.180, 'CA' : 7.358,
                 'OG1': 38.230, 'O'  : 27.475, 'N'  : 4.937,
                 'CG2': 62.021},
         'PHE': {'C'  :  0.876, 'CE1': 33.974, 'CB' : 30.431,
                 'CA' :  4.790, 'CG' :  1.388, 'O'  : 28.205,
                 'N'  :  6.352, 'CZ' : 34.444, 'CD1': 19.752,
                 'CD2': 18.580, 'CE2': 32.812},
         'ALA': {'CB' : 70.070, 'CA' : 11.256, 'C'  :  0.717,
                 'O'  : 27.891, 'N'  :  7.782},
         'HIS': {'C'  :  0.923, 'CE1': 32.391, 'CB' : 31.658,
                 'CA' :  6.665, 'CG' :  1.280, 'O'  : 27.750,
                 'N'  :  7.365, 'CD2': 27.492, 'ND1': 15.030,
                 'NE2': 36.247},
         'MET': {'C'  :  0.590, 'CB' : 28.322, 'CA' :  3.961,
                 'CG' : 29.567, 'CE' : 74.891, 'N'  :  7.345,
                 'O'  : 27.658, 'SD' : 29.550},
         'ASP': {'C'  :  0.972, 'CB' : 37.245, 'CA' : 10.246,
                 'CG' :  3.665, 'O'  : 25.235, 'N'  :  3.251,
                 'OD1': 37.525, 'OD2': 33.852},
         'LEU': {'C'  :  0.793, 'CB' : 21.950, 'CA' :  2.852,
                 'CG' : 14.129, 'O'  : 27.978, 'N'  :  5.231,
                 'CD1': 61.814, 'CD2': 59.042},
         'ARG': {'C'  :  0.829, 'CB' : 28.313, 'CA' :  2.930,
                 'CG' : 27.099, 'NE' : 23.452, 'O'  : 28.419,
                 'CD' : 23.936, 'CZ' :  1.719, 'NH1': 50.063,
                 'NH2': 46.094, 'N'  :  7.109},
         'TRP': {'C'  :  0.856, 'CZ2': 31.924, 'CB' : 28.556,
                 'CA' :  3.339, 'CG' :  1.337, 'O'  : 28.269,
                 'N'  :  5.626, 'CH2': 34.844, 'CE3': 20.285,
                 'CE2':  4.521, 'CD2':  3.335, 'CZ3': 34.289,
                 'NE1': 32.834, 'CD1': 23.722},
         'GLU': {'C'  :  0.987, 'CB' : 27.019, 'CA' :  5.321,
                 'CG' : 41.443, 'O'  : 29.082, 'CD' :  3.502,
                 'OE2': 35.857, 'N'  :  4.529, 'OE1': 32.074},
         'TYR': {'C'  :  0.777, 'CD2': 20.732, 'OH' : 56.712,
                 'CB' : 29.172, 'CA' :  4.380, 'CG' :  1.517,
                 'O'  : 28.824, 'N'  :  5.930, 'CZ' :  6.352,
                 'CD1': 19.637, 'CE1': 30.502, 'CE2': 30.358}}


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
        self.ranMS = ranMS
        self.ranAS = ranAS

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


    def __relExposure( self, key='AS' ):
        """
        Calculate how exposed an atom is relative to the same
        atom in a GLY-XXX-GLY tripeptide, an approximation of
        the unfolded state.
        key - MS or AS
        -> relAS - list of relative solvent accessible surfaces
           relMS - list of relative molecular surface exposure
        """
        if not key=='MS' and not key=='AS':
            raise SurfaceRacer_Error,\
                  'Incorrect key for relative exposiure: %s '%key
             
        rel = []
        i=0
        for a in self.model.atoms:
            atom = a['name']
            ## a little crude, we don't care if it is terminal
            ## => unrealistic exposures for termini
            if atom == 'OXT':
                atom='O'
            resi = a['residue_name']
            if self.result[key][i] != 0:
                rel += [ (self.result[key][i] /
                          self.ranAS[resi][atom]) *100 ]
            else:
                rel += [0.0]
            i += 1
        
        self.result['rel'+key ] = Numeric.array(rel)
        
                    
    
    def isFailed( self ):
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_result( self.output )
        if self.probe == 1.4:
            self.__relExposure('MS')
            self.__relExposure('AS')
        else:
            print "No relative accessabilities calculated when using other prob radius then 1.4 A"

#######
## test
if __name__ == '__main__':

    from Biskit import PDBModel
    import Biskit.tools as T
    import glob
    import Biskit.mathUtils as MA

    print "Loading PDB..."

    f = glob.glob( T.testRoot()+'/lig_pcr_00/pcr_00/*_1_*pdb' )[1]
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

    print 'Relative MS of atoms 10 to 20 atoms:', r['relMS'][10:20]

    print 'Relative AS of atoms 10 to 20 atoms:', r['relAS'][10:20]
