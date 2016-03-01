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
## last $Date$
## $Revision$
## last $Author:  
"""
Display a Ramachandran plot for a list of PDBModels with
the same atom contents.
"""

from Biskit import ColorSpectrum as CS
from Biskit.MatrixPlot import Legend
from Biskit.PDBDope import PDBDope 
from Biskit import EHandler

import Biskit.tools as T

import numpy as N

try:
    import biggles
except:
    biggles = 0


class Ramachandran:

    def __init__( self, models, name=None, profileName='relAS', verbose=1 ):
        """
        @param models: List of models display a Ramachandran plot for
        @type  models: [ PDBModel ] OR PDBModel
        @param name: model name, will show up in plot
        @type  name: str
        @param profileName: name of profile to use for coloring
                            (default: 'relAS')
        @type  profileName: str
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        """
        if not biggles:
            raise ImportError, 'biggles module could not be imported.'

        if type(models) != type([]):
            models = [ models ]

        self.psi = []
        self.phi = []

        self.gly = []
        self.pro = []

        self.prof=[]
        self.profileName =  profileName

        self.name=name

        self.verbose = verbose

        # calculate angles, profiles ...
        self.calc( models )
        
        self.prof = N.ravel(self.prof)
        self.gly = N.ravel(self.gly)
        self.pro = N.ravel(self.pro)


    def calc( self, models ):
        """
        Calculate angles, profiles and other things needed.
        
        @param models: List of models 
        @type  models: [ PDBModel ]
        """              
        res_count = 0
        for m in models:
            
            ## add profile if not there
            if self.profileName:
                self.prof += [ self.calcProfiles( m ) ]
                     
            ## calclate phi and psi angles for model
            self.phi_and_psi( m )

            ## get list with GLY and PRO residue indices
            gly_atomInd = m.indices(lambda a: a['residue_name']=='GLY')
            gly_resInd  = N.array( m.atom2resIndices( gly_atomInd ) )
            pro_atomInd = m.indices(lambda a: a['residue_name']=='PRO')
            pro_resInd  = N.array( m.atom2resIndices( pro_atomInd ) )    
            self.gly.append( gly_resInd + res_count )
            self.pro.append( pro_resInd + res_count )
            res_count += m.lenResidues()
            

    def calcProfiles( self, m ):
        """
        Calculate needed profiles.

        @param m: PDBModel to calculate data for
        @type  m: PDBModel
        """
        if self.verbose: print "Initiating PDBDope..."
        d = PDBDope( m )
                
        if not self.profileName in m.atoms.keys():
            
            if self.profileName in ['MS', 'AS', 'curvature', 'relAS', 'relMS']:
                if self.verbose: print "Adding SurfaceRacer profile...",
                d.addSurfaceRacer()
                            
            if self.profileName in ['relASA']:
                if self.verbose: print "Adding WhatIf ASA...",
                d.addASA()
                        
            if self.profileName in ['density']:
                if self.verbose: print "Adding surface density...",
                d.addDensity()
                                     
        if not self.profileName in m.residues.keys():
                    
            if self.profileName in ['cons_abs', 'cons_max', 'cons_ent']:
                if self.verbose: print "Adding conservation data...",
                d.addConservation()
                                               
            if self.profileName in ['ASA_total', 'ASA_sc', 'ASA_bb']:
                if self.verbose: print "Adding WhatIf ASA...",
                d.addASA()
      
        if self.verbose: print 'Done.'

        ## convert atom profiles to average residue profile
        if self.profileName in m.atoms.keys():
            prof = []
            aProfile = m.profile( self.profileName )
            resIdx =  m.resIndex().tolist()
            resIdx += [ m.lenAtoms()]
            for i in range(len(resIdx)-1):
                prof += [ N.average( N.take(aProfile, range(resIdx[i],
                                                            resIdx[i+1]), 0 ), 0 )]
        else:
            prof = m.profile( self.profileName )

        return prof

    
    def phi_and_psi( self, model ):
        """
        Calculate phi and psi torsion angles for all
        residues in model::
        
          phi - rotation about the N-CA bond
              - last position in a chain = None
          psi - rotation about CA-C
              - first position in a chain = None          

        @param model: PDBModel
        @type  model: PDBModel 
        """
        for c in range( model.lenChains(breaks=1) ):
            cModel = model.takeChains( [c], breaks=1 )

            xyz = cModel.xyz

            xyz_CA =  N.compress( cModel.maskCA(), xyz, 0 )
            xyz_N  =  N.compress( cModel.mask( ['N'] ), xyz, 0 )
            xyz_C  =  N.compress( cModel.mask( ['C'] ), xyz, 0 )

            ## phi: c1 - N
            ##      c2 - CA
            ##      c3 - C
            ##      c4 - N of next residue
            for i in range( len(xyz_N)-1 ):
                self.phi += [self.dihedral( xyz_N[i], xyz_CA[i],
                                            xyz_C[i], xyz_N[i+1] )]
            self.phi += [None]

            ## psi: c1 - C of previous residue  
            ##      c2 - N
            ##      c3 - CA
            ##      c4 - C
            self.psi += [None]
            for i in range( 1, len(xyz_N) ):
                self.psi += [self.dihedral( xyz_C[i-1], xyz_N[i],
                                            xyz_CA[i], xyz_C[i] )]


    def dihedral( self, coor1, coor2, coor3, coor4 ):
        """
        Calculates the torsion angle of a set of four atom coordinates.
        The dihedral angle returned is the angle between the projection
        of i1-i2 and the projection of i4-i3 onto a plane normal to i2-i3.

        @param coor1: coordinates
        @type  coor1: [float]
        @param coor2: coordinates
        @type  coor2: [float]
        @param coor3: coordinates
        @type  coor3: [float]
        @param coor4: coordinates
        @type  coor4: [float]        
        """
        vec21 = coor2 - coor1
        vec32 = coor3 - coor2
        L = N.cross( vec21, vec32 )
        L_norm = N.sqrt(sum(L**2))

        vec43 = coor4 - coor3
        vec23 = coor2 - coor3
        R = N.cross( vec43, vec23 )
        R_norm = N.sqrt(sum(R**2))

        S     = N.cross( L, R )
        angle = sum( L*R ) / ( L_norm * R_norm )

        ## sometimes the value turns out to be ever so little greater than 
        ## one, to prevent N.arccos errors for this, set angle = 1.0
        if angle >  1.0: angle = 1.0
            
        if angle < -1.0: angle = -1.0

        angle = N.arccos(angle) *180/N.pi
        if sum(S*vec32) < 0.0:
            angle = -angle

        return angle

    
    def ramachandran( self ):
        """
        Create all the ramachandran plot points.

        @return: list of biggles.Point objects (all the points of the
                 plot)and a biggles.Inset object (property scale).
        @rtype: [ biggles.Point ], biggles.Inset
        """
        p = []

        ## calculate colors and create a legend if a property is given
        if self.profileName:
            palette = CS('plasma', 0, 100)
            col     = palette.color_array( self.prof )

            legend = Legend(  palette.legend() )
            inset  = biggles.Inset((1.1, 0.60), (1.2, .97), legend)

        else:
            col = ['black']*len(self.phi)
            inset = None

        ## add data points to plot
        for i in range(len(self.phi)):
            ## don't add termini - has missing angles
            if self.phi[i] and self.psi[i]:
                if i in self.gly:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="star", size=1, color=col[i] )]
                elif i in self.pro:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="filled square", size=1,
                                         color=col[i] )]
                else:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="filled circle", size=1,
                                         color=col[i] )]
        return p, inset


    def ramachandran_background( self ):
        """
        Creates a background (favoured regions) for a ramachandran plot.

        @return: list of biggles.Point objects
        @rtype: [ biggles.Point ]
        """
        bg = []
        mat = biggles.read_matrix( T.dataRoot() +
                                   '/biggles/ramachandran_bg.dat')
        x, y = N.shape(mat)
        for i in range(x):
            for j in range(y):
                if mat[i,j] < 200:
                    a = (360./y)*j    - 180
                    b = (360./x)*(x-i)- 180
                    bg += [ biggles.Point( a, b, type="dot" )]
        return bg


    def show( self, fileName=None ):
        """
        Show ramachandran plot.
        """
        plot = biggles.FramedPlot()
        plot.xrange = (-180., 180.)
        plot.yrange = (-180., 180.)
        plot.xlabel = "$\Phi$"
        plot.ylabel = "$\Psi$"      
        
        if self.name:
            plot.title = self.name

        ## add allowed regions
        bg_plot = self.ramachandran_background( )
        for p in bg_plot:
            plot.add( p )

        ## add ramachandran phi, psi valies
        points, inset = self.ramachandran(  )
        for p in points:
            plot.add(p)
        if inset:
            plot.add( inset )
      
        plot.add( biggles.PlotLabel( 1.14, 0.55, self.profileName, size=2) )
        plot.add( biggles.PlotLabel( 1.1, 0.45, "GLY star", size=2) )
        plot.add( biggles.PlotLabel( 1.12, 0.40, "PRO square", size=2) )
      
        plot.show()

        if fileName:
            plot.write_eps( fileName )


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    def test_Ramachandran(self):
        """Ramachandran test"""
        self.traj = T.load( T.testRoot()+'/lig_pcr_00/traj.dat' )

        self.traj.ref.atoms.set('mass', self.traj.ref.masses() ) 

        self.mdl = [ self.traj[0], self.traj[11] ]
        self.mdl = [ md.compress( md.maskProtein() ) for md in self.mdl ]

        self.rama = Ramachandran( self.mdl , name='test', profileName='mass',
                                  verbose=self.local)

        self.psi = N.array( self.rama.psi )

        if self.local:
            self.rama.show()
            
        r = N.sum( N.compress( N.logical_not(N.equal(self.psi, None)),
                               self.psi ) )
        self.assertAlmostEqual( r, -11717.909796797909, 2 )

 
if __name__ == '__main__':

    BT.localTest()
    
    

    
