##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2009 Peter Schmidtke, Raik Gruenberg
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


import re,numpy as npy
import sys,os
from math import *

from biskit.core.localpath import LocalPath
import biskit.tools as tools


class EZDParser:
    """ 
    Parse Electron Density Maps (EZD format).

    Electron Density Maps can be obtained from the EDS in Upsala:
    http://eds.bmc.uu.se/eds/
    """

    #: keys of all atom profiles that are read directly from the PDB file
#    EZD_KEYS = ['comment', 'type', 'origin', 'extent',
#                'scale', 'alpha', 'beta', 'gamma',
#                'a', 'b', 'c']

    def __init__( self, source=None ):
        """
        @param source: path to EZD file
        @type  source: str OR file handle
        """

        ## Raik says: this is always overridden by the "self.source=source" 
        if type( source ) is str and os.path.isfile( source ):
            self.source = LocalPath( source )

        self.source = source

        ## my suggestion (LocalPath is a bit overkill here):
        self.source = tools.absfile( source )

        #try :
        f=tools.gzopen(self.source)
        self.__jumpComment(f)
        cell=self.__readCell(f)
        self.a=cell[0]
        self.b=cell[1]
        self.c=cell[2]
        self.alpha=cell[3]
        self.beta=cell[4]
        self.gamma=cell[5]

        self.origin=self.__readOrigin(f)
        self.extent=self.__readExtent(f)
        self.gridsize=self.__readGridSize(f)
        self.scale=self.__readScale(f)


        self.EDGrid=self.__readGrid(f)
        self.__setCartesianCoords()
        self.sd=npy.std(self.EDGrid)
        self.mean=npy.mean(self.EDGrid)
        f.close()
        #except Exception:
            #   print "Tes un gros naze...voila"

    def __jumpComment(self,f):
        for i in range(2) : f.readline()         #move the file pointer to the cell entry

    def __readCell(self,f):
        line=f.readline().split()       #change this in re....
        if(line[0]=='CELL'):
            return npy.array([float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6])])
        return None

    def __readOrigin(self,f):
        line=f.readline().split()
        if(line[0]=='ORIGIN'):          #change this to re
            return(npy.array([float(line[1]),float(line[2]),float(line[3])]))
        return None

    def __readExtent(self,f):
        line=f.readline().split()  #change to re
        if(line[0]=='EXTENT'):
            return npy.array([int(line[1]),int(line[2]),int(line[3])])
        return None

    def __readGridSize(self,f):
        line=f.readline().split()
        if(line[0]=='GRID'):        #change to re
            return npy.array([int(line[1]),int(line[2]),int(line[3])])
        return None

    def __readScale(self,f):
        line=f.readline().split()
        if(line[0]=='SCALE'):   #change to re
            tmp=f.readline()
            return float(line[1])
        return None

    def __readGrid(self,f):
        grid=npy.zeros(self.extent)
        #initiate xyz counter for reading the grid data
        z=0
        y=0
        x=0

        r=re.compile('END')
        #for count in range(n_entries/7) :i
        while 1:
            rl=f.readline()
            if len(r.findall(rl))>0 or len(rl)==0 :
                break
            c=rl.split(" ")
    #       if(len(c)!=7) : 
    #               print "last line...nearly done...oooohh"


            for i in range(len(c)):
                grid[x][y][z]=float(c[i])
                x+=1
                if x >= self.extent[0]:
                    x=0
                    y+=1
                    if y >=self.extent[1]:
                        y=0
                        z+=1
        return grid

    def getEDGrid(self):
        return self.EDGrid

    def __setCartesianCoords(self):
        abc=npy.array([self.a,self.b,self.c])
        abg=npy.array([self.alpha,self.beta,self.gamma])*(pi/180.0)
        delta=abc/self.gridsize
        self.resolution=delta
        #abc=abc/gridsize
        self.gridspacing=delta        
        #next according to C source code of VMD for origin calculation
        xaxis=npy.zeros(3)
        yaxis=npy.zeros(3)
        zaxis=npy.zeros(3)

        xaxis[0]=delta[0]

        yaxis[0]=cos(abg[2])*delta[1]
        yaxis[1]=sin(abg[2])*delta[1]

        z1=cos(abg[1])
        z2=sin(cos(abg[0]-cos(abg[1])*cos(abg[2])))/sin(abg[2])
        z3=sqrt(1.0-z1*z1-z2*z2)

        zaxis[0]=z1*delta[2]
        zaxis[1]=z2*delta[2]
        zaxis[2]=z3*delta[2]
        originc=npy.zeros(3)
        originc[0]=xaxis[0]*self.origin[0]+yaxis[0]*self.origin[1]+zaxis[0]*self.origin[2]
        originc[1]=yaxis[1]*self.origin[1]+zaxis[1]*self.origin[2]
        originc[2]=zaxis[2]*self.origin[2]

        self.__cartesianOrigin=originc.copy()

#         xaxis=xaxis*(self.extent[0]-1)
#         yaxis=yaxis*(self.extent[1]-1)
#         zaxis=zaxis*(self.extent[2]-1)
#         
        ext=self.extent.copy()
        ext.resize(4)
        ext[3]=4
        self.__cartesianPositionGrid=npy.zeros(ext)

        #next according to C source code on http://www.ccl.net/chemistry/resources/messages/2000/11/22.004-dir/index.html

        #optimize this with an index array!!!!!!
        #o=open("out.pdb","w")

        for x_i in range(self.extent[0]):
            for y_i in range(self.extent[1]):
                for z_i in range(self.extent[2]):
                    fx=(float(x_i))
                    fy=(float(y_i))
                    fz=(float(z_i))
                    cx=originc[0]+(xaxis[0]*fx+yaxis[0]*fy+zaxis[0]*fz)
                    cy=originc[1]+(yaxis[1]*fy+zaxis[1]*fz)
                    cz=originc[2]+(zaxis[2]*fz)
                    self.__cartesianPositionGrid[x_i][y_i][z_i]=npy.array([fx,fy,fz,float(self.EDGrid[x_i][y_i][z_i])])
        #o.close()

    def getCartesianOrigin(self):
        return self.__cartesianOrigin

    def getCartesianPositionGrid(self):
        return self.__cartesianPositionGrid

    def getIdxFromCart(self,coords):
        """Get numpy array indices for a cartesian position in the ED grid"""
        return npy.round((coords-self.__cartesianOrigin)/self.resolution)


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.LONG ]

    def prepare( self ):
        self.f = os.path.join( tools.testRoot(), 'eedensity/3tgi.ezd.gz' )
        self.p = None

    def test_all( self ):
        """EZDParser test"""

        if self.local:
            print("parsing ...", end=' ')
        self.p = self.p or EZDParser( self.f )

        if self.local:
            print("Done")

        self.origin = self.p.getCartesianOrigin()
        target = [ -38.92757372, -115.8054999 ,  -36.27633333]

        ## compare result with expected result but allow for numeric deviations
        self.assertTrue( npy.sum(self.origin - target) < 1e-8 )

        self.grid = self.p.getCartesianPositionGrid()
        self.assertEqual( npy.shape( self.grid ), (129, 105, 104, 4) )

        self.i = self.p.getIdxFromCart( [10, 10, 10] )
        self.assertTrue( npy.all( self.i == [  82.,  212.,   80.] ) )



if __name__ == '__main__':

    BT.localTest()
