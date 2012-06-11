##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2012 Peter Schmidtke, Daniel Alvarez Garcia, Raik Gruenberg
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
## last $Author: Peter Schmidtke & Daniel Alvarez $
### last $Date: 11/12/2011$
## $Revision: $



"""This module includes currently all CUDA related calculation

Using cuda with Biskit and BiskitFF is fully optional but if you have a CUDA
enabled GPU and cuda installed on your machine we recommend you to try it as
you can achieve very nice performance increases especially on distance matrix
calculations, energy calculations with grids and so on

"""


import pycuda.driver as drv
import pycuda.tools
import pycuda.autoinit
import numpy as npy
import numpy.linalg as la
from pycuda.compiler import SourceModule
import sys


#import pycuda.driver as drv

# Initialize CUDA
drv.init()
import pycuda.tools
version=drv.get_version()
num_dev=drv.Device.count()
print >>sys.stderr, "Identified CUDA devices %d"%(num_dev)
selection=0
free=npy.zeros([num_dev],dtype="int32")
for dnumber in range(num_dev):
        dev = drv.Device(dnumber)
        context=dev.make_context()
        device = context.get_device()
        mem=drv.mem_get_info()
        free[dnumber]=mem[0]/(1024.**2)
        context.detach()
#print mem
w=int(npy.argmax(free))
print >>sys.stderr, "Using GPU device : ",w
dev=drv.Device(w)

context=dev.make_context()
device=context.get_device()
import atexit
atexit.register(context.pop)

bsizes1D=npy.arange(1,512)

def __getBlockNGridSizes(n1,n2,force=None):
    """A small helper function to determine possible grid and block sizes for CUDA based on the size of the input vectors in 2D"""
    bsizes=npy.arange(1,23) #list of possible block sizes (x,y) without having the risk to go beyond 512 threads per bloc, not that 23 is not included, so max 484 threads are executed
    if not force is None and n1%force==0 and n2%force==0:
        bsizes=npy.array([force])
    bx=bsizes[max(npy.argwhere(n1%bsizes==0))] # take the maximum blocksize that is possible
    by=bsizes[max(npy.argwhere(n2%bsizes==0))]
    
    gx=n1/bx                        #determine inherent grid size
    gy=n2/by

    return((int(bx),int(by),int(gx),int(gy)))
    

def __getBlockNGridSizes1D(n1):
    """Helper function for getting good blockSizes (1D) and grid sizes (1D)"""
    bx=bsizes1D[max(npy.argwhere(n1%bsizes1D==0))]
    gx=n1/bx
    
    return((int(bx),1, int(gx), 1))
    
    
def getCoulomb1D(xyz1,xyz2,c1,c2):
    """Calculate the pairwise Coulomb potential between xyz1 and all positions in xyz2
    @param : a 3d numpy array (1 element)
    @param : a 3d numpy array (n elements)"""
    
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    c1=c1.astype(npy.dtype("float32"))      #cast the array to float32
    c2=c2.astype(npy.dtype("float32"))      #cast the array to float32
    x1=xyz1[0]
    y1=xyz1[1]
    z1=xyz1[2]
    #print x1
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])  
    
    (bX,gX)=__getBlockNGridSizes1D(len(xyz2))
    
    #input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    #output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    #itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    #mem=drv.mem_get_info()  #get total and free memory information
    #free=mem[0]             #get only the free memory on the card
   
    #if itemsize+len(xyz1)*len(xyz2)*4<free :
        #Here find the C source code for the CUDA distance calculation
    coulombGrid = SourceModule("""
    __global__ void coulombRow(float *dest, float x1,float *x2, float y1,float *y2, float z1,float *z2, float c1, float *c2) //__global__ indicates that it is a kernel
    {
    const int x = blockIdx.x*blockDim.x+threadIdx.x;
    float d=sqrtf((x1-x2[x])*(x1-x2[x])+(y1-y2[x])*(y1-y2[x])+(z1-z2[x])*(z1-z2[x])); //Amber Coulomb potential
    dest[x]=(c1*c2[x])/d; //Amber Coulomb potential -> see http://ambermd.org/Questions/units.html
    }
    """)
    
    coulombGrid = coulombGrid.get_function("coulombRow")
    
    #create an empty array storing the resulting distance matrix
    coulombMat = npy.zeros(len(x2),dtype="float32")
    
    #the cuda call
    coulombGrid(drv.Out(coulombMat), x1,drv.In(x2), y1,drv.In(y2),z1,drv.In(z2),c1,drv.In(c2),block=(bX,1,1), grid=(gX,1))
    
    #reshape the resulting array to a distance matrix
    
    return coulombMat


def getvdWMat(xyz1,xyz2,nbpidx,atypes1,atypes2,a_coefs,b_coefs,hbond_a_coefs, hbond_b_coefs,ntypes,excluded_idx_start_list,excluded_atoms,sym=1,correct=0):
    """Calculate the pairwise Coulomb potential
    @param : a 3d numpy array
    @param : a 3d numpy array """
    
    excluded_idx_start_list=excluded_idx_start_list.astype(npy.dtype("uint32"))
    excluded_atoms=excluded_atoms.astype(npy.dtype("int32"))
    nbpidx=nbpidx.astype(npy.dtype("int16"))
    #print nbpidx
    atypes1=atypes1.astype(npy.dtype("uint16"))
    atypes2=atypes2.astype(npy.dtype("uint16"))
    #excluded_atoms-=excluded_atoms
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    a_coefs=a_coefs.astype(npy.dtype("float32"))      #cast the array to float32
    b_coefs=b_coefs.astype(npy.dtype("float32"))      #cast the array to float32
    hbond_a_coefs=hbond_a_coefs.astype(npy.dtype("float32"))      #cast the array to float32
    hbond_b_coefs=hbond_b_coefs.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    #print __getBlockNGridSizes(len(xyz1),len(xyz2))
    
    input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
     
    #if itemsize+len(xyz1)*len(xyz2)*4<free :
    vdWGrid = SourceModule("""
    __global__ void vdWGridSym(float *dest, float *x1,float *x2, float *y1,float *y2, float *z1,float *z2, short *ico,unsigned short *atypes1, unsigned short *atypes2,float *acoefs, float *bcoefs, float *hbond_a_coef, float *hbond_b_coef, unsigned int *ex_start, int *ex_list, unsigned short ntypes) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        float acoef,bcoef,d;
        int i;
        dest[x+yt]=0.0;
        if(y>x){
            for(i=ex_start[x];i<ex_start[x+1];i++){
                    if(ex_list[i]-1==y) return;
            }
            const int icotmp=ico[ntypes*(atypes1[x]-1)+atypes2[y]-1];
            
            if(icotmp <0) {
                acoef=hbond_a_coef[icotmp*(-1)-1];
                bcoef=hbond_b_coef[icotmp*(-1)-1];
                d=sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y]));
                dest[x+yt]=acoef/powf(d,12)-bcoef/powf(d,10); //Lennard Jone potential -> see http://ambermd.org/Questions/vdwequation.pdf
            }
            else {
                acoef=acoefs[icotmp-1];
                bcoef=bcoefs[icotmp-1];
                d=sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y]));
                dest[x+yt]=acoef/powf(d,12)-bcoef/powf(d,6); //Lennard Jone potential -> see http://ambermd.org/Questions/vdwequation.pdf
            }
            
        }
    }
    __global__ void vdWGridAsym(float *dest, float *x1,float *x2, float *y1,float *y2, float *z1,float *z2, short *ico,unsigned short *atypes1, unsigned short *atypes2,float *acoefs, float *bcoefs, float *hbond_a_coef, float *hbond_b_coef, unsigned int *ex_start, int *ex_list, unsigned short ntypes) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        const int totX= gridDim.x*blockDim.x;   //get the dimension of the whole grid on the x axis
        float acoef,bcoef,d;
        int i;
        dest[x+yt]=0.0;
        int end=x+1;
        if(x+1>=totX) end=totX-1;
        for(i=ex_start[x];i<ex_start[end];i++){
                if(ex_list[i]-1==y) return;
        }
        const int icotmp=ico[ntypes*(atypes1[x]-1)+atypes2[y]-1];
        if(icotmp <0) {
                acoef=hbond_a_coef[icotmp*(-1)-1];
                bcoef=hbond_b_coef[icotmp*(-1)-1];
                d=sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y]));
                dest[x+yt]=acoef/powf(d,12)-bcoef/powf(d,10); //Lennard Jone potential -> see http://ambermd.org/Questions/vdwequation.pdf
        }
        else {
                acoef=acoefs[icotmp-1];
                bcoef=bcoefs[icotmp-1];
                d=sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y]));
                dest[x+yt]=acoef/powf(d,12)-bcoef/powf(d,6); //Lennard Jone potential -> see http://ambermd.org/Questions/vdwequation.pdf
        }
        
        
    }
    """)
    
    if sym==1 :
        vdWGrid = vdWGrid.get_function("vdWGridSym")
    else :
        vdWGrid = vdWGrid.get_function("vdWGridAsym")
    
    #create an empty array storing the resulting distance matrix
    vdWMat = npy.zeros(len(x1)*len(x2),dtype="float32")
    
    
    #the cuda call
    vdWGrid(drv.Out(vdWMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),drv.In(nbpidx), drv.In(atypes1),drv.In(atypes2),drv.In(a_coefs),drv.In(b_coefs),drv.In(hbond_a_coefs),drv.In(hbond_b_coefs),drv.In(excluded_idx_start_list),drv.In(excluded_atoms),npy.uint16(ntypes),block=(bX,bY,1), grid=(gX,gY))
    
    #reshape the resulting array to a distance matrix
    #vdWMat=vdWMat.reshape(len(x2),len(x1))
    return vdWMat



def getCoulombMat(xyz1,xyz2,c1,c2,excluded_idx_start_list,excluded_atoms,sym=1,correct=0,ystart=0):
    """Calculate the pairwise Coulomb potential
    @param : a 3d numpy array
    @param : a 3d numpy array """
    
    excluded_idx_start_list=excluded_idx_start_list.astype(npy.dtype("uint32"))
    excluded_atoms=excluded_atoms.astype(npy.dtype("int32"))
    #excluded_atoms-=excluded_atoms
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    c1=c1.astype(npy.dtype("float32"))      #cast the array to float32
    c2=c2.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    #print __getBlockNGridSizes(len(xyz1),len(xyz2))
    
    input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
   
    #if itemsize+len(xyz1)*len(xyz2)*4<free :
    coulombGrid = SourceModule("""
    __global__ void coulombGridSym(float *dest, float *x1,float *x2, float *y1,float *y2, float *z1,float *z2, float *c1, float *c2,unsigned int *ex_start, int *ex_list, int ystart) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        const int yglobal=y+ystart;
        int i;
        dest[x+yt]=0.0;
        if(y>x){
            for(i=ex_start[x];i<ex_start[x+1];i++){
                    if(ex_list[i]-1==yglobal) return;
            }
            dest[x+yt]=(c1[x]*c2[y])/sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y])); //Amber Coulomb potential -> see http://ambermd.org/Questions/units.html
        }
    }
    __global__ void coulombGridAsym(float *dest, float *x1,float *x2, float *y1,float *y2, float *z1,float *z2, float *c1, float *c2,unsigned int *ex_start, int *ex_list,int ystart) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        const int yglobal=y+ystart;
        const int totX= gridDim.x*blockDim.x;   //get the dimension of the whole grid on the x axis
        int i;
        dest[x+yt]=0.0;
        int end=x+1;
        if(x+1>=totX) {
            end=totX-1;
            for(i=ex_start[x];i<=ex_start[end];i++) if(ex_list[i]-1==yglobal){
//                dest[x+yt]=ex_list[i];
                return;
            }
        }
        else for(i=ex_start[x];i<ex_start[end];i++) if(ex_list[i]-1==yglobal) {
//            dest[x+yt]=ex_list[i];
            return;
        }
        //dest[x+yt]=ex_list[i];
        dest[x+yt]=(c1[x]*c2[y])/sqrtf((x1[x]-x2[y])*(x1[x]-x2[y])+(y1[x]-y2[y])*(y1[x]-y2[y])+(z1[x]-z2[y])*(z1[x]-z2[y])); //Amber Coulomb potential -> see http://ambermd.org/Questions/units.html
    }
    """)
    
    if sym==1 :
        coulombGrid = coulombGrid.get_function("coulombGridSym")
    else :
        coulombGrid = coulombGrid.get_function("coulombGridAsym")
    
    #create an empty array storing the resulting distance matrix
    coulombMat = npy.zeros(len(x1)*len(x2),dtype="float32")
    
    #the cuda call
    coulombGrid(drv.Out(coulombMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),drv.In(c1),drv.In(c2),drv.In(excluded_idx_start_list),drv.In(excluded_atoms),npy.uint32(ystart),block=(bX,bY,1), grid=(gX,gY))
    
    #reshape the resulting array to a distance matrix
    coulombMat=coulombMat.reshape(len(x2),len(x1))
    return coulombMat


def getGridEnergies(protein, ACoefs, BCoefs, cGrid, vGrid, xResGrid,pocket):
    """2d grid computation of coulomb and vdw potential"""
    
    if cGrid.data.shape!= vGrid.data.shape or npy.any(cGrid.origin!=vGrid.origin):
        print "The vdw grid and coulomb grid don't have the same shape and/or origin. Breaking up cuda calculation"
        return
    
    xyz=protein.xyz
    origin=npy.array(cGrid.origin,"float32")
    delta=cGrid.delta[0]

    xyz = xyz.astype(npy.float32)
    x=npy.copy(xyz[:,0])
    y=npy.copy(xyz[:,1])
    z=npy.copy(xyz[:,2])
    pxyz = pocket.astype(npy.float32)
    px=npy.copy(pxyz[:,0])
    py=npy.copy(pxyz[:,1])
    pz=npy.copy(pxyz[:,2])
    #print x.shape
    #x.shape=(len(x),1)
    #y.shape=(len(y),1)
    #z.shape=(len(z),1)
    #print z.shape
    charges_gpu=protein.AmberCharges.astype(npy.float32)
    acoefs_gpu=ACoefs.astype(npy.float32)
    bcoefs_gpu=BCoefs.astype(npy.float32)
    bZ=1
    (bX,bY,gX,gY)=__getBlockNGridSizes(cGrid.data.shape[1],cGrid.data.shape[2])
    #print __getBlockNGridSizes(cGrid.data.shape[1],cGrid.data.shape[2])
    #(bX,bY,bZ,gX,gY)=(8,8,2,8,8)  #this is fixed for 64 cubed grids here!!!
    
    rescGrid=npy.zeros([cGrid.data.shape[1],cGrid.data.shape[2]],"float32")
    resvGrid=npy.zeros([vGrid.data.shape[1],vGrid.data.shape[2]],"float32")
    
    sqGridEnergies = SourceModule("""
    __global__ void getsqGridEnergies(float* cOut, float *vOut, float *xcoords,float *ycoords,float *zcoords, float *px,float *py,float *pz,float *charges, float *acoefs, float *bcoefs, unsigned int natoms, unsigned int np,float *origin, float delta,unsigned short gridX){
        const int gx = blockIdx.x*blockDim.x+threadIdx.x;
        const int gy = blockIdx.y*blockDim.y+threadIdx.y;
        const int gyt = gy*gridDim.x*blockDim.x;
        int i;
        float d=10.0;
        float posX=origin[0]+(gridX)*delta;
        float posY=origin[1]+gy*delta;
        float posZ=origin[2]+gx*delta;
        cOut[gx+gyt]=0.0;
        vOut[gx+gyt]=0.0;
        float vdwTmp=0.0;
        float minD=100.0;
        for(i=0;i<np;i++){
            d=sqrtf(powf(posX-px[i],2)+powf(posY-py[i],2)+powf(posZ-pz[i],2));
            if(d<minD) minD=d;
        }
        if(minD>3.0) {
                cOut[gx+gyt]=0.0;
                vOut[gx+gyt]=100.0;
                return;
        }
        for(i=0;i<natoms && d>1.0;i++){
            d=sqrtf(powf(posX-xcoords[i],2)+powf(posY-ycoords[i],2)+powf(posZ-zcoords[i],2));
                cOut[gx+gyt]+=charges[i]/d;
                vdwTmp+=acoefs[i]/powf(d,12)-bcoefs[i]/powf(d,10);
        }
        if(d<1.0) vOut[gx+gyt]=100;
        else vOut[gx+gyt]=vdwTmp;
    }
    """)
    function = sqGridEnergies.get_function("getsqGridEnergies")
    function(drv.Out(rescGrid),drv.Out(resvGrid), drv.In(x),drv.In(y), drv.In(z),drv.In(px),drv.In(py),drv.In(pz),drv.In(charges_gpu),drv.In(acoefs_gpu),drv.In(bcoefs_gpu),npy.uint(len(xyz)),npy.uint(len(pocket)),drv.In(origin),npy.float32(delta),npy.uint16(xResGrid),block=(bX,bY,bZ), grid=(gX,gY))
    return rescGrid,resvGrid
    #(bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    
    
def getGridEnergies2(fullSystem, cGrid, vGrid, xResGrid,pocket,probe_atype,subSystemMask=None):
    """2d grid computation of coulomb and vdw potential version 2 for big systems with solvent"""
    
    if cGrid.data.shape!= vGrid.data.shape or npy.any(cGrid.origin!=vGrid.origin):
        print "The vdw grid and coulomb grid don't have the same shape and/or origin. Breaking up cuda calculation"
        return
    if subSystemMask is None:
        subSystemMask=npy.ones(length(fullSystem),dtype="bool")
    
    
    xyz=fullSystem.xyz.copy()
    nproteinAtoms=npy.sum(subSystemMask)
    nsolventAtoms=len(xyz)-nproteinAtoms
    #solventXyz=fullSystem.xyz[]
    origin=npy.array(cGrid.origin,"float32")
    delta=cGrid.delta[0]

    xyz = xyz.astype(npy.float32)
    x=npy.copy(xyz[:,0])
    y=npy.copy(xyz[:,1])
    z=npy.copy(xyz[:,2])
    pxyz = pocket.astype(npy.float32)
    px=npy.copy(pxyz[:,0])
    py=npy.copy(pxyz[:,1])
    pz=npy.copy(pxyz[:,2])
    #print x.shape
    #x.shape=(len(x),1)
    #y.shape=(len(y),1)
    #z.shape=(len(z),1)
    #print z.shape
    charges_gpu=fullSystem.AmberCharges.astype(npy.float32)
    acoefs_gpu=fullSystem.AmberTopology_a_coefs.copy()
    acoefs_gpu=acoefs_gpu.astype(npy.float32)
    bcoefs_gpu=fullSystem.AmberTopology_b_coefs.copy()
    bcoefs_gpu=bcoefs_gpu.astype(npy.float32)
    ico_gpu=fullSystem.AmberTopology_ico.copy()
    atypes_gpu=fullSystem.AmberTopology_A_types.copy()
    ntypes=fullSystem.AmberTopology_n_a_types
    #self.AmberTopology_N_of_excluded_atoms (n)
    #self.AmberTopology_Excluded_idx_start_list (n)
    #self.AmberTopology_Excluded_atoms
    
    #self.AmberTopology_A_types (atype indices) (n)
    #self.AmberTopology_ico (ntypesxntypes)
    #self.AmberTopology_a_coefs
    #self.AmberTopology_b_coefs
    
    bZ=1
    (bX,bY,gX,gY)=__getBlockNGridSizes(cGrid.data.shape[1],cGrid.data.shape[2])
    #print __getBlockNGridSizes(cGrid.data.shape[1],cGrid.data.shape[2])
    #(bX,bY,bZ,gX,gY)=(8,8,2,8,8)  #this is fixed for 64 cubed grids here!!!
    
    rescGrid=npy.zeros([cGrid.data.shape[1],cGrid.data.shape[2]],"float32")
    resvGrid=npy.zeros([vGrid.data.shape[1],vGrid.data.shape[2]],"float32")
    
    sqGridEnergies = SourceModule("""
    __global__ void getsqGridEnergies(float* cOut, float *vOut, float *xcoords,float *ycoords,float *zcoords, float *px,float *py,float *pz,float *charges,  short *atom_types,short *ico,float *acoefs,float *bcoefs,float *origin,unsigned short gridX,float delta,unsigned int np,unsigned int nProteinAtoms,unsigned int nSolventAtoms,short probe_atom_type_index,unsigned short ntypes){
        const int gx = blockIdx.x*blockDim.x+threadIdx.x;
        const int gy = blockIdx.y*blockDim.y+threadIdx.y;
        const int gyt = gy*gridDim.x*blockDim.x;
        
        /** init the results arrays with 0.0*/
        cOut[gx+gyt]=0.0;
        vOut[gx+gyt]=0.0;
        
        float posX=origin[0]+(gridX)*delta;
        float posY=origin[1]+gy*delta;
        float posZ=origin[2]+gx*delta;
        int i;
        int j;
        float d=0.0;
        float minD=100.0;
        float acoef=0.0;
        float bcoef=0.0;
        float eps=78.4;// #dielectric constant
        float A=-8.5525;
        float l=0.003627;
        float k=7.7839;
        float B=eps-A;
        cOut[gx+gyt]=0.0;
        /** get the minimum distance of the current grid point from the predefined pocket */
        for(i=0;i<nProteinAtoms;i++){
            d=sqrtf(powf(posX-xcoords[i],2)+powf(posY-ycoords[i],2)+powf(posZ-zcoords[i],2));
            if(minD>d)minD=d;
        }
        /*for(i=0;i<np;i++){
            d=sqrtf(powf(posX-px[i],2)+powf(posY-py[i],2)+powf(posZ-pz[i],2));
            if(minD>d)minD=d;
        }*/

        /**if the grid point is too far away from the pocket, lets put a penalty on the vdw term*/
        if(minD>3.0) {
                cOut[gx+gyt]=0.0;
                vOut[gx+gyt]=1000.0;
                return;
        }
        else {
        /**calculate the electric field on a griven grid point */
            for(i=0;i<nProteinAtoms;i++){
                d=sqrtf(powf(posX-xcoords[i],2)+powf(posY-ycoords[i],2)+powf(posZ-zcoords[i],2));
                if(d>2.2){
                    const int icotmp=ico[ntypes*(atom_types[i]-1)+probe_atom_type_index-1];
                    if(icotmp<0) {  /*normally we should read the 10-12 potential from asol and bsol here, but I never saw a topology where this was used and where the coeff was non zero...so claimed useless for now, but should be implemented*/
                        acoef=0.0;
                        bcoef=0.0;
                    }
                    acoef=acoefs[icotmp-1];
                    bcoef=bcoefs[icotmp-1];
    //                cE=

                    cOut[gx+gyt]+=-(charges[i]*((B+A)*expf(2*l*d*B)+(k*l*d*(B*B)+k*B+2.*k*A)*expf(l*d*B)+(k*k)*A)/(((d*d)*(B*B)+2.*(d*d)*A*B+(d*d)*(A*A))*expf(2.0*l*d*B)+(2.*k*(d*d)*A*B+2*k*(d*d)*(A*A))*expf(l*d*B)+(k*k)*(d*d)*(A*A)));
                    vOut[gx+gyt]+=acoef/powf(d,12)-bcoef/powf(d,10);
                }
                else {
                    cOut[gx+gyt]+=0.0;
                    vOut[gx+gyt]+=100.;
                }
            }
        }
    }
    
    
    """)
    #print origin
    function = sqGridEnergies.get_function("getsqGridEnergies")
    function(drv.Out(rescGrid),drv.Out(resvGrid), drv.In(x),drv.In(y), drv.In(z),drv.In(px),drv.In(py),drv.In(pz),drv.In(charges_gpu),drv.In(atypes_gpu),drv.In(ico_gpu),drv.In(acoefs_gpu),drv.In(bcoefs_gpu),drv.In(origin),npy.uint16(xResGrid),npy.float32(delta),npy.uint32(len(pocket)),npy.uint32(nproteinAtoms),npy.uint32(nsolventAtoms),npy.int16(probe_atype),npy.uint16(ntypes),block=(bX,bY,bZ), grid=(gX,gY))
    return rescGrid,resvGrid
    #print rescGrid,resvGrid

#    cE=-(charges[i]*((B+A)*npy.exp(2.*l*d*B)+(k*l*d*(B**2)+k*B+2.*k*A)*npy.exp(l*d*B)+(k**2)*A)/(((d**2)*(B**2)+2.*(d**2)*A*B+(d**2)*(A**2))*npy.exp(2.0*l*d*B)+(2.*k*(d**2)*A*B+2*k*(d**2)*(A**2))*npy.exp(l*d*B)+(k**2)*(d**2)*(A**2))
    #(bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))


def mapSolventProtInteractionEn2Grid(sGrid,fullSystem,protein):
    """"""
    import math,copy
    nProteinAtoms=len(protein.xyz)
    nSolventAtoms=len(fullSystem.xyz)-nProteinAtoms

    print "preparing for %d x %d calculations"%(nProteinAtoms,nSolventAtoms)
    ressGrid=npy.zeros([sGrid.data.shape[0]*sGrid.data.shape[1]*sGrid.data.shape[2]],"float32")
#    ressGrid=npy.ones_like(sGrid.data)-1
    nGridDim=sGrid.data.shape[1]
    origin=npy.array(sGrid.origin,"float32")
    delta=sGrid.delta[0]
    
    cudaCode="""
    __device__ inline float myAtomicAdd(float* address, float value)
    {
        float old = value;  
        while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
    };
    
    __global__ void mapSolventProteinIntEnergy(float *sOut, float *origin, float delta,unsigned short nResGridDim,float *proteinXyz, float *solventXyz,unsigned int nProteinAtoms, unsigned int nSolventAtoms, float *solventCharges,float *proteinCharges, unsigned short *atypesSolvent, unsigned short *atypesProtein, float *acoefs, float *bcoefs, short *ico,unsigned short ntypes){
        const int gx = blockIdx.x*blockDim.x+threadIdx.x;
        const int gy = blockIdx.y*blockDim.y+threadIdx.y;
        const int gyt = gy*gridDim.x*blockDim.x;
        const int curIdx=gx+gyt;
        
        
        
        float protX=proteinXyz[gy*3];
        float protY=proteinXyz[gy*3+1];
        float protZ=proteinXyz[gy*3+2];
        
        float solvX=solventXyz[gx*3];
        float solvY=solventXyz[gx*3+1];
        float solvZ=solventXyz[gx*3+2];
        
        float solventCharge=solventCharges[gx];
        float proteinCharge=proteinCharges[gy];
        
        float originX=origin[0];
        float originY=origin[1];
        float originZ=origin[2];        
        
        float acoef=0.0;
        float bcoef=0.0;
        
        float d=sqrtf(powf(protX-solvX,2)+powf(protY-solvY,2)+powf(protZ-solvZ,2));
        
        int outX=rintf((solvX-originX)/delta);
        int outY=rintf((solvY-originY)/delta);
        int outZ=rintf((solvZ-originZ)/delta);
        
        const int icotmp=ico[ntypes*(atypesProtein[gy]-1)+atypesSolvent[gx]-1];
        if(icotmp<0) {  /*normally we should read the 10-12 potential from asol and bsol here, but I never saw a topology where this was used and where the coeff was non zero...so claimed useless for now, but should be implemented*/
            acoef=0.0;  
            bcoef=0.0;
        }
        acoef=acoefs[icotmp-1];
        bcoef=bcoefs[icotmp-1];

        if(outX>=0 && outX<nResGridDim && outY>=0 && outY<nResGridDim && outZ>=0 && outZ<nResGridDim && d>2.0){
            myAtomicAdd(&sOut[outX*nResGridDim*nResGridDim+outY*nResGridDim+outZ],solventCharge*proteinCharge/d+(acoef/powf(d,12)-bcoef/powf(d,6)));

            __syncthreads();
            //sOut[outX*nGridDim+outY*nGridDim+outZ]=-1.0;
            
        }
        else return;
    }
    """
    module=SourceModule(cudaCode)
    function = module.get_function("mapSolventProteinIntEnergy")
    
    acoefs_gpu=fullSystem.AmberTopology_a_coefs.copy()
    acoefs_gpu=acoefs_gpu.astype(npy.float32)
    bcoefs_gpu=fullSystem.AmberTopology_b_coefs.copy()
    bcoefs_gpu=bcoefs_gpu.astype(npy.float32)
    ico_gpu=fullSystem.AmberTopology_ico.copy()
    ntypes=fullSystem.AmberTopology_n_a_types
    
    
    # split nProt and nSolvent in bunches of 16x16 2d blocks
    bZ=1
    base=16
    maxSolventMultiplicator=200
    maxProteinMultiplicator=maxSolventMultiplicator
    nStepsProt=math.floor(nProteinAtoms/float(base))
    nStepsSolvent=math.floor(nSolventAtoms/float(base))
    
    
    solventStartIndex=nProteinAtoms
    curSolventStartIndex=solventStartIndex
    if(base*maxSolventMultiplicator<nSolventAtoms):
        curSolventEndIndex=curSolventStartIndex+base*maxSolventMultiplicator
    else :
        curSolventEndIndex=solventStartIndex+nSolventAtoms
    solventEndIndex=nProteinAtoms+nSolventAtoms
    
    while(curSolventEndIndex<=solventEndIndex):
        solventXyz=fullSystem.xyz[curSolventStartIndex:curSolventEndIndex].copy()
        solventCharges=fullSystem.AmberCharges[curSolventStartIndex:curSolventEndIndex].copy()
        solventAtypes=fullSystem.AmberTopology_A_types[curSolventStartIndex:curSolventEndIndex].copy()      
            
        solventXyz=solventXyz.astype("float32")
        solventCharges=solventCharges.astype("float32")
        solventAtypes=solventAtypes.astype("uint16")
        
        proteinStartIndex=0
        curProteinStartIndex=proteinStartIndex
        if(base*maxProteinMultiplicator<nProteinAtoms):
            curProteinEndIndex=curProteinStartIndex+base*maxProteinMultiplicator
        else :
            curProteinEndIndex=nProteinAtoms
        proteinEndIndex=nProteinAtoms
        
        while(curProteinEndIndex<=proteinEndIndex):
            proteinXyz=fullSystem.xyz[curProteinStartIndex:curProteinEndIndex].copy()
            proteinCharges=copy.deepcopy(fullSystem.AmberCharges[curProteinStartIndex:curProteinEndIndex])
            proteinAtypes=fullSystem.AmberTopology_A_types[curProteinStartIndex:curProteinEndIndex].copy()
            proteinXyz=proteinXyz.astype("float32")
            proteinCharges=proteinCharges.astype("float32")
            proteinAtypes=proteinAtypes.astype("uint16")
            
            #print  __getBlockNGridSizes(len(solventXyz),len(proteinXyz),force=base)
            (bX,bY,gX,gY) = __getBlockNGridSizes(len(solventXyz),len(proteinXyz),force=base)
            print nGridDim
            function(drv.InOut(ressGrid),drv.In(origin),npy.float32(delta),npy.uint16(nGridDim),drv.In(proteinXyz),drv.In(solventXyz),npy.uint32(len(proteinXyz)),npy.uint32(len(solventXyz)),drv.In(solventCharges),drv.In(proteinCharges),drv.In(solventAtypes),drv.In(proteinAtypes),drv.In(acoefs_gpu),drv.In(bcoefs_gpu),drv.In(ico_gpu),npy.uint16(ntypes),block=(bX,bY,bZ), grid=(gX,gY))
            
            curProteinStartIndex=curProteinEndIndex
            if(curProteinEndIndex==nProteinAtoms):
                curProteinEndIndex+=1
            else :
                if(base*maxProteinMultiplicator+curProteinStartIndex<nProteinAtoms):
                    curProteinEndIndex=curProteinStartIndex+base*maxProteinMultiplicator
                else :
                    curProteinEndIndex=nProteinAtoms

        curSolventStartIndex=curSolventEndIndex        
        if curSolventEndIndex==nProteinAtoms+nSolventAtoms:
            curSolventEndIndex+=1
        else :
            if(base*maxSolventMultiplicator+curSolventStartIndex<nSolventAtoms):
                curSolventEndIndex=curSolventStartIndex+base*maxSolventMultiplicator
            else :
                curSolventEndIndex=nProteinAtoms+nSolventAtoms
    
#    ressGrid.shape=(nGridDim,nGridDim,nGridDim)
    return ressGrid.reshape((nGridDim,nGridDim,nGridDim))



def getDistanceMatrixSq3D(xyz1,xyz2,maxCalcSize=3000) :
    """Get a pairwise distance matrix using GPU based calculations
   @param : a 3d numpy array
   @param : a 3d numpy array """
    distMat = npy.zeros([len(xyz2),len(xyz1)],dtype="float32")
    
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32

    
    #input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    #output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    #itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
    #mem=drv.mem_get_info()  #get total and free memory information
    #free=mem[0]             #get only the free memory on the card
   
    #if itemsize+len(xyz1)*len(xyz2)*4<free :
        #Here find the C source code for the CUDA distance calculation
    distMatGrid = SourceModule("""
    __global__ void distMatGrid(float *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2) //__global__ indicates that it is a kernel
    {
    const int x = blockIdx.x*blockDim.x+threadIdx.x;
    const int y = blockIdx.y*blockDim.y+threadIdx.y;
    //const int xt = x*gridDim.y*blockDim.y;
    const int yt = y*gridDim.x*blockDim.x;
    dest[x+yt]=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]); //euclidean squared distance between two points in 3D space
    //dest[y+xt]=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]); //euclidean squared distance between two points in 3D space
    
    }
    """)
    
    distMatGrid = distMatGrid.get_function("distMatGrid")
    
    #create an empty array storing the resulting distance matrix
    
    
    #the cuda call
    #print "total : ",len(xyz1),len(xyz2)
    curXStartIdx=0
    curXEndIdx=len(xyz1)
    
    if(curXEndIdx>maxCalcSize):
        curXEndIdx=maxCalcSize
    cnt=0
    while curXEndIdx<=len(xyz1):
        
        curYStartIdx=0
        curYEndIdx=len(xyz2)    
        if(curYEndIdx>maxCalcSize):
            curYEndIdx=maxCalcSize
        
        while curYEndIdx<=len(xyz2):
            x1=npy.copy(xyz1[curXStartIdx:curXEndIdx,0])
            y1=npy.copy(xyz1[curXStartIdx:curXEndIdx,1])
            z1=npy.copy(xyz1[curXStartIdx:curXEndIdx,2])
            
            x2=npy.copy(xyz2[curYStartIdx:curYEndIdx,0])
            y2=npy.copy(xyz2[curYStartIdx:curYEndIdx,1])
            z2=npy.copy(xyz2[curYStartIdx:curYEndIdx,2])
            
            tmpMat=npy.zeros([len(x2),len(x1)],dtype="float32")
            (bX,bY,gX,gY)=__getBlockNGridSizes(len(x1),len(x2))
            
            #print (bX,bY,gX,gY), len(x1),len(x2)
            cnt+=1
            distMatGrid(drv.Out(tmpMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),block=(bX,bY,1), grid=(gX,gY))
            
            distMat[curYStartIdx:curYEndIdx,curXStartIdx:curXEndIdx]=tmpMat.copy()
            curYStartIdx=curYEndIdx
            if (curYEndIdx==len(xyz2)):
                curYEndIdx+=1
            else :
                if(len(xyz2)-curYEndIdx)<maxCalcSize:
                    curYEndIdx=len(xyz2)
                else :
                    curYEndIdx+=maxCalcSize
            
            
        curXStartIdx=curXEndIdx
        if(curXEndIdx==len(xyz1)):
            curXEndIdx+=1
        else :
            if (len(xyz1)-curXEndIdx)<maxCalcSize:
                curXEndIdx=len(xyz1)
            else :
                curXEndIdx+=maxCalcSize
        
      
    #reshape the resulting array to a distance matrix
    #distMat=distMat.reshape(len(x2),len(x1))
    return distMat
    
    #else:
        #print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        #return None
    
    

def getDistanceMatrixSq3DDEPRECATED(xyz1,xyz2) :
    """Get a pairwise distance matrix using GPU based calculations
   @param : a 3d numpy array
   @param : a 3d numpy array """
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    
    print "dmat between %d x %d "%(len(xyz1),len(xyz2))
    
    input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
   
    #if itemsize+len(xyz1)*len(xyz2)*4<free :
        #Here find the C source code for the CUDA distance calculation
    distMatGrid = SourceModule("""
    __global__ void distMatGrid(float *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2) //__global__ indicates that it is a kernel
    {
    const int x = blockIdx.x*blockDim.x+threadIdx.x;
    const int y = blockIdx.y*blockDim.y+threadIdx.y;
    const int yt = y*gridDim.x*blockDim.x;
    dest[x+yt]=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]); //euclidean squared distance between two points in 3D space
    }
    """)
    
    distMatGrid = distMatGrid.get_function("distMatGrid")
    
    #create an empty array storing the resulting distance matrix
    distMat = npy.zeros(len(x1)*len(x2),dtype="float32")
    
    #the cuda call
    distMatGrid(drv.Out(distMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),block=(bX,bY,1), grid=(gX,gY))
    
    #reshape the resulting array to a distance matrix
    distMat=distMat.reshape(len(x2),len(x1))
    return distMat
    
    #else:
        #print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        #return None


def getAsaOverlap((spntsPocket,elAsaPocket),(spntsFoundPocket,elAsaFoundPocket),maxCalcSize=10000):
    """given a set of surface points (3d coords and associated portion of surface), this function compares this set (coordinates) with a second set of coordinates and sends back the intersection between both in terms of ASA"""
    xyz1=spntsPocket.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=spntsFoundPocket.astype(npy.dtype("float32"))      #cast the array to float32
    
    asa1=elAsaPocket.astype(npy.dtype("float32"))
    #asa2=elAsaFoundPocket.astype(npy.dtype("float32"))
    
    
    source=SourceModule("""
    __global__ void comparePoints(float *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2, float *asa1) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        const float xsq=(a1[x]-a2[y])*(a1[x]-a2[y]);
        const float ysq=(b1[x]-b2[y])*(b1[x]-b2[y]);
        const float zsq=(c1[x]-c2[y])*(c1[x]-c2[y]);
        const float thresh=1e-10;
        if(xsq <thresh && ysq<thresh && zsq<thresh) dest[x+yt]=asa1[x];
    }
    """)
    
    comparePoints = source.get_function("comparePoints")
    
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])

    totalRes=0.0
    curStartIdx=0
    curEndIdx=len(xyz1)
    if curEndIdx>maxCalcSize:
        curEndIdx=maxCalcSize
    #here we consider only the pocket as potentially big (can also happen with the ligand defined pocket .,...but I have no time.... and no neurons left
    while curEndIdx<=len(xyz2):
        x2=npy.copy(xyz2[curStartIdx:curEndIdx,0])
        y2=npy.copy(xyz2[curStartIdx:curEndIdx,1])
        z2=npy.copy(xyz2[curStartIdx:curEndIdx,2])
        (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(x2))
        res=npy.zeros([len(xyz1),len(x2)],dtype="float32")
        comparePoints(drv.InOut(res),drv.In(x1),drv.In(x2),drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),drv.In(asa1),block=(bX,bY,1), grid=(gX,gY))
    
        totalRes+=npy.sum(res)
        curStartIdx=curEndIdx
        
        if(curEndIdx==len(xyz2)):
            curEndIdx+=1
        else :
            if (len(xyz2)-curEndIdx)<maxCalcSize:
                curEndIdx=len(xyz2)
            else :
                curEndIdx+=maxCalcSize

    return totalRes


def getGoodGridPointsAlphaSpheres(receptor,foundPocket,pocketGrid,minPocket,gridSpacing,origin,result):
    """given a 3d input grid with ones inside (pocketGrid), this function checks of each grid point as within an alpha sphere (foundPocket) and within at minimum the equilibrium vdw distance from nearby receptor atoms """
    
    refDist=5.5
    #calculate the distance matrix between the pocket and the receptor
    dm=getDistanceMatrixSq3D(receptor.xyz,foundPocket.xyz)
    pocket_ref=receptor.take(npy.where(npy.sum(dm<=refDist**2,axis=0))[0])
    eqDist=pocket_ref["vdwEqDistance"]
    
    xyz1=pocket_ref.xyz.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    xyz2=foundPocket.xyz.astype(npy.dtype("float32"))
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    radii=foundPocket["radius"]
    minPocket=minPocket.astype(npy.dtype("float32"))
    origin=origin.astype(npy.dtype("float32"))
    source=SourceModule("""__global__ void checkGridPoints(short *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2, float *radii, float* eqDist, float gridSpacing, float *origin, short npocket, short nreceptor,short runningX,short pGridYDim, short pGridZDim) 
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt= x*pGridZDim;
        const int destIdx = runningX*pGridYDim*pGridZDim+yt+y;
        /*const int y = blockIdx.x*blockDim.x+threadIdx.x;
        const int x = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt= x*gridDim.y*blockDim.y;
        const int destIdx = runningX*pGridYDim*pGridZDim+yt+y;*/
        const float curx=origin[0]+runningX*gridSpacing;    //x level is given through arguments
        const float cury=origin[1]+x*gridSpacing;       //x (CUDA) = y (resultGrid)
        const float curz=origin[2]+y*gridSpacing;       //y(CUDA) = z (resultGrid
        int i;
        short count=0;
        
        for(i=0;i<npocket;i++){
            if(((a2[i]-curx)*(a2[i]-curx)+(b2[i]-cury)*(b2[i]-cury)+(c2[i]-curz)*(c2[i]-curz))<(radii[i]*radii[i])) {
                count++;
                break;
            }
        }
        if(count){
            for(i=0;i<nreceptor;i++){
                if(((a1[i]-curx)*(a1[i]-curx)+(b1[i]-cury)*(b1[i]-cury)+(c1[i]-curz)*(c1[i]-curz))<(eqDist[i]*eqDist[i])) {
                    return;
                }
            }
            dest[destIdx]=1;
        }
        // else do nothing        
        
    }
    """)
    comparePoints = source.get_function("checkGridPoints")

    (bX,bY,gX,gY)=__getBlockNGridSizes(pocketGrid.shape[1],pocketGrid.shape[2])
    xDim=pocketGrid.data.shape[0]
    yDim=pocketGrid.data.shape[1]
    zDim=pocketGrid.data.shape[2]
    pocketGrid.data=pocketGrid.data.reshape([pocketGrid.data.size])
    for x in xrange(xDim):
        
    #tmpGrid=npy.zeros([pocketGrid.data.shape[1],pocketGrid.data.shape[2]],dtype="float32")
        comparePoints(drv.InOut(pocketGrid.data),drv.In(x1),drv.In(x2),drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),drv.In(radii),drv.In(eqDist),npy.float32(gridSpacing),drv.In(origin), npy.int16(len(foundPocket)),npy.int16(len(pocket_ref)),npy.int16(x),npy.int16(yDim),npy.int16(zDim),block=(bX,bY,1), grid=(gX,gY))
    
    
    
    pocketGrid.data=pocketGrid.data.reshape([pocketGrid.shape[0],pocketGrid.shape[1],pocketGrid.shape[2]])

    #pocketGrid.writeDX("pocket.dx")
    #print tmpGrid,len(foundPocket)
    #print tmpGrid
    #pocketGrid[x,:,:]=npy.copy(tmpGrid)
    
    #now transform this to cartesian coordinates
    ok=npy.where(pocketGrid.data>0)
    result.xyz=npy.zeros([len(ok[0]),3])
    #print len(ok[0]),ok[0]
    for i in xrange(len(ok[0])):
        #print origin[1]+ok[1][i]*gridSpacing
        result.xyz[i,0]=origin[0]+ok[0][i]*gridSpacing
        result.xyz[i,1]=origin[1]+ok[1][i]*gridSpacing
        result.xyz[i,2]=origin[2]+ok[2][i]*gridSpacing

    #print result.xyz[:,1]

def getSurfacePts(xyz1, molxyz, radii, excluded_idx, lig_ref=None, pocket=None) :
    maxCalcSize=50000
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    molxyz=molxyz.astype(npy.dtype("float32"))      #cast the array to float32
    radii=radii.astype(npy.dtype("float32"))
    suppMask=None    
        
    if not lig_ref==None and not pocket==None :
        lig_ref = lig_ref.astype(npy.dtype("float32"))
        pocket  = pocket.astype(npy.dtype("float32"))
    excluded_idx=excluded_idx.astype(npy.dtype("uint16"))


    #input_size=xyz1.size*xyz1.dtype.itemsize+molxyz.size*molxyz.dtype.itemsize  #get the total memory necessary for the input vector
    #output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    #itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
    
    #mem=drv.mem_get_info()  #get total and free memory information
    #free=mem[0]             #get only the free memory on the card
   
    #if itemsize<free :
    #Here find the C source code for the CUDA distance calculation
    isASAPoint = SourceModule("""
    __global__ void is_ASA_point(short *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2, float *radii, unsigned short nb,unsigned short *ex_idx) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        unsigned short y;
        float correct=-0.0001;
        dest[x]=0;
        for(y=0;y<nb;y++){
            if(x!=ex_idx[x]){
                float d=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]);
                if(d<radii[y]*radii[y]+correct) {
                    dest[x]=0;
                    return;
                }
            }
        }
        dest[x]=1;
    }
    __global__ void is_ASA_point_with_ref(short *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2,
        float *radii, unsigned short nb,unsigned short *ex_idx, float *pocket,float *ligand, unsigned short len_ligand, float *p) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        unsigned short y;
        float correct=-0.0001;
        dest[x]=0;
        float d=0.0;
        float x1,y1,z1,x2,y2,z2;
        float dmin=15000.0;
        unsigned int ligidx=0;
        unsigned short pidx=0;
        pidx=ex_idx[x];
        x2=pocket[pidx*3];
        y2=pocket[pidx*3+1];
        z2=pocket[pidx*3+2];
        for(y=0;y<len_ligand;y++){
            x1=ligand[y*3];
            y1=ligand[y*3+1];
            z1=ligand[y*3+2];
            d=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
            if(d<dmin) {
                dmin=d;
                ligidx=y;
            }
        }
        for(y=0;y<nb;y++){
            if(x!=ex_idx[x]){
                d=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]);
                    if(d<=radii[y]*radii[y]+correct) {
                        dest[x]=0;
                        return;
                }
            }
        }
        
        x1=ligand[ligidx*3];
        y1=ligand[ligidx*3+1];
        z1=ligand[ligidx*3+2];
        d=(x1-a1[x])*(x1-a1[x])+(y1-b1[x])*(y1-b1[x])+(z1-c1[x])*(z1-c1[x]);
        if(dmin+20>=d) dest[x]=1;
        else dest[x]=0;
    }
    """)
    #create an empty array storing the resulting mask
    mask = npy.zeros(len(xyz1),dtype="int16")
    
    if not lig_ref==None and not pocket==None :
        isASAPoint = isASAPoint.get_function("is_ASA_point_with_ref")
    else:
        isASAPoint = isASAPoint.get_function("is_ASA_point")
    curStartIdx=0
    curEndIdx=len(xyz1)
    if curEndIdx>maxCalcSize:
        curEndIdx=maxCalcSize
    while curEndIdx<=len(xyz1):
        x1=npy.copy(xyz1[curStartIdx:curEndIdx,0])
        y1=npy.copy(xyz1[curStartIdx:curEndIdx,1])
        z1=npy.copy(xyz1[curStartIdx:curEndIdx,2])
        
        x2=npy.copy(molxyz[:,0])
        y2=npy.copy(molxyz[:,1])
        z2=npy.copy(molxyz[:,2])
        tmpexidx=npy.copy(excluded_idx[curStartIdx:curEndIdx])
        (bX,bY,gX,gY)=__getBlockNGridSizes1D(x1.shape[0])
        #print curEndIdx,len(xyz1),tmpexidx.shape,x1.shape,molxyz.shape
        
        tmpmask=npy.zeros(len(x1),dtype="int16")
        if not lig_ref==None and not pocket==None :
            isASAPoint(drv.Out(tmpmask), drv.In(x1),drv.In(x2),    drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2), drv.In(radii), npy.uint16(molxyz.shape[0]), drv.In(tmpexidx),drv.In(pocket),drv.In(lig_ref),
                    npy.uint16(lig_ref.shape[0]),  block=(bX,bY,1), grid=(gX,gY))
        else:
            isASAPoint(drv.Out(tmpmask), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),                         drv.In(radii), npy.uint16(molxyz.shape[0]), drv.In(tmpexidx),block=(bX,bY,1), grid=(gX,gY))
            
        mask[curStartIdx:curEndIdx]=tmpmask.copy()
        curStartIdx=curEndIdx
        
        if(curEndIdx==len(xyz1)):
            curEndIdx+=1
        else :
            if (len(xyz1)-curEndIdx)<maxCalcSize:
                curEndIdx=len(xyz1)
            else :
                curEndIdx+=maxCalcSize
    
    
    ##thresh=100000
    ##while bX*gX>thresh:
        ##if not suppMask==None:
            ##npy.hstack((suppMask,getSurfacePts(xyz1[:thresh],molxyz,radii, excluded_idx,lig_ref,pocket)))
        ##else:
            ##suppMask=getSurfacePts(xyz1[:thresh],molxyz,radii, excluded_idx,lig_ref,pocket)
        ##xyz1=npy.copy(xyz1[thresh:])
        ##(bX,bY,gX,gY)=__getBlockNGridSizes1D(xyz1.shape[0])

    
    
    
    
    #if not suppMask==None:
        #mask=npy.hstack((suppMask,mask))
    return mask
    
    #else:
        #print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        #return None

def getSurfacePtsDEPRECATED(xyz1, molxyz, radii, excluded_idx, lig_ref=None, pocket=None) :
    """Get a pairwise distance matrix using GPU based calculations
   @param : a 3d numpy array
   @param : a 3d numpy array
   @param : a 1d numpy array 
   @param : a 1d numpy array """
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    molxyz=molxyz.astype(npy.dtype("float32"))      #cast the array to float32
    radii=radii.astype(npy.dtype("float32"))
    maskSmall=None
    if not lig_ref==None and not pocket==None :
        lig_ref = lig_ref.astype(npy.dtype("float32"))
        pocket  = pocket.astype(npy.dtype("float32"))
    excluded_idx=excluded_idx.astype(npy.dtype("uint16"))
    (bX,bY,gX,gY)=__getBlockNGridSizes1D(xyz1.shape[0])
    if bX==1 and gX>65536:
        shift=7
        (bX,bY,gX,gY)=__getBlockNGridSizes1D(xyz1.shape[0]-shift)
        maskSmall=getSurfacePts(xyz1[-shift:, :], molxyz, radii[-shift:], excluded_idx[-shift:], lig_ref, pocket)
        xyz1=xyz1[:-shift, :]
        radii=radii[-shift:]
        excluded_idx=excluded_idx[-shift:]
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(molxyz[:,0])
    y2=npy.copy(molxyz[:,1])
    z2=npy.copy(molxyz[:,2])

    input_size=xyz1.size*xyz1.dtype.itemsize+molxyz.size*molxyz.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
    
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
   
    #if itemsize<free :
    #Here find the C source code for the CUDA distance calculation
    isASAPoint = SourceModule("""
    __global__ void is_ASA_point(short *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2, float *radii, unsigned short nb,unsigned short *ex_idx) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        unsigned short y;
        float correct=-0.0001;
        dest[x]=0;
        for(y=0;y<nb;y++){
            if(x!=ex_idx[x]){
                float d=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]);
                if(d<radii[y]*radii[y]+correct) {
                    dest[x]=1;
                    return;
                }
            }
        }
        dest[x]=0;
    }
    __global__ void is_ASA_point_with_ref(short *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2,
        float *radii, unsigned short nb,unsigned short *ex_idx, float *pocket,float *ligand, unsigned short len_ligand, float *p) //__global__ indicates that it is a kernel
    {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        unsigned short y;
        float correct=-0.0001;
        dest[x]=0;
        float d=0.0;
        float x1,y1,z1,x2,y2,z2;
        float dmin=15000.0;
        unsigned int ligidx=0;
        unsigned short pidx=0;
        pidx=ex_idx[x];
        x2=pocket[pidx*3];
        y2=pocket[pidx*3+1];
        z2=pocket[pidx*3+2];
        for(y=0;y<len_ligand;y++){
            x1=ligand[y*3];
            y1=ligand[y*3+1];
            z1=ligand[y*3+2];
            d=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
            if(d<dmin) {
                dmin=d;
                ligidx=y;
            }
        }
        for(y=0;y<nb;y++){
            if(x!=ex_idx[x]){
                d=(a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y]);
                    if(d<radii[y]*radii[y]+correct) {
                        dest[x]=1;
                        return;
                }
            }
        }
        x1=ligand[ligidx*3];
        y1=ligand[ligidx*3+1];
        z1=ligand[ligidx*3+2];
        d=(x1-a1[x])*(x1-a1[x])+(y1-b1[x])*(y1-b1[x])+(z1-c1[x])*(z1-c1[x]);
        if(d<dmin+4.0) dest[x]=0;
        else dest[x]=1;
    }
    """)
    #create an empty array storing the resulting mask
    mask = npy.zeros(len(x1),dtype="int16")
    
    if not lig_ref==None and not pocket==None :
        isASAPoint = isASAPoint.get_function("is_ASA_point_with_ref")
        isASAPoint(drv.Out(mask), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),
                    drv.In(radii), npy.uint16(molxyz.shape[0]), drv.In(excluded_idx),drv.In(pocket),drv.In(lig_ref),
                    npy.uint16(lig_ref.shape[0]),  block=(bX,bY,1), grid=(gX,gY))
    else:
        isASAPoint = isASAPoint.get_function("is_ASA_point")
        isASAPoint(drv.Out(mask), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2), 
                    drv.In(radii), npy.uint16(molxyz.shape[0]), drv.In(excluded_idx),block=(bX,bY,1), grid=(gX,gY))
    
    #reshape the resulting array to a distance matrix
    if not maskSmall==None:
        mask=npy.hstack((mask, maskSmall))
    return mask
    
    #else:
        #print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        #return None


def getDistanceMatrix3D(xyz1,xyz2) :
    """Get a pairwise distance matrix using GPU based calculations
   @param : a 3d numpy array
   @param : a 3d numpy array """
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    
    input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
   
    if itemsize+len(xyz1)*len(xyz2)*4<free :
        #Here find the C source code for the CUDA distance calculation
        distMatGrid = SourceModule("""
        __global__ void distMatGrid(float *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2) //__global__ indicates that it is a kernel
        {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        dest[x+yt]=sqrtf((a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y])); //euclidean squared distance between two points in 3D space
        }
        """)
        
        distMatGrid = distMatGrid.get_function("distMatGrid")
        
        #create an empty array storing the resulting distance matrix
        distMat = npy.zeros(len(x1)*len(x2),dtype="float32")
        
        #the cuda call
        distMatGrid(drv.Out(distMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),block=(bX,bY,1), grid=(gX,gY))
        
        #reshape the resulting array to a distance matrix
        distMat=distMat.reshape(len(x2),len(x1))
        return distMat
        
    else:
        print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        return None
    
    
#def getDistanceMatrix3D(xyz1,xyz2,stride=1000):
    #"""helper function for calculation of distance matrix"""
    
    #dm=npy.zeros([len(xyz2)][len(xyz2)])
    
    
    #return dm    

def getDistanceMatrix3DKernel(xyz1,xyz2) :
    """Get a pairwise distance matrix using GPU based calculations
   @param : a 3d numpy array
   @param : a 3d numpy array """
    xyz1=xyz1.astype(npy.dtype("float32"))      #cast the array to float32
    xyz2=xyz2.astype(npy.dtype("float32"))      #cast the array to float32
    x1=npy.copy(xyz1[:,0])
    y1=npy.copy(xyz1[:,1])
    z1=npy.copy(xyz1[:,2])
    
    x2=npy.copy(xyz2[:,0])
    y2=npy.copy(xyz2[:,1])
    z2=npy.copy(xyz2[:,2])
    
    (bX,bY,gX,gY)=__getBlockNGridSizes(len(xyz1),len(xyz2))
    
    input_size=xyz1.size*xyz1.dtype.itemsize+xyz2.size*xyz2.dtype.itemsize  #get the total memory necessary for the input vector
    output_size=(x1.size*x2.size)*x1.dtype.itemsize                         #get the total memory necessary for the output vector
    itemsize=input_size+output_size                                         #get the total amount of memory necessary for the operation
   
    mem=drv.mem_get_info()  #get total and free memory information
    free=mem[0]             #get only the free memory on the card
   
    if itemsize+len(xyz1)*len(xyz2)*4<free :
        #Here find the C source code for the CUDA distance calculation
        distMatGrid = SourceModule("""
        __global__ void distMatGrid(float *dest, float *a1,float *a2, float *b1,float *b2, float *c1,float *c2) //__global__ indicates that it is a kernel
        {
        const int x = blockIdx.x*blockDim.x+threadIdx.x;
        const int y = blockIdx.y*blockDim.y+threadIdx.y;
        const int yt = y*gridDim.x*blockDim.x;
        dest[x+yt]=sqrtf((a1[x]-a2[y])*(a1[x]-a2[y])+(b1[x]-b2[y])*(b1[x]-b2[y])+(c1[x]-c2[y])*(c1[x]-c2[y])); //euclidean squared distance between two points in 3D space
        }
        """)
        
        distMatGrid = distMatGrid.get_function("distMatGrid")
        
        #create an empty array storing the resulting distance matrix
        distMat = npy.zeros(len(x1)*len(x2),dtype="float32")
        
        #the cuda call
        distMatGrid(drv.Out(distMat), drv.In(x1),drv.In(x2), drv.In(y1),drv.In(y2),drv.In(z1),drv.In(z2),block=(bX,bY,1), grid=(gX,gY))
        
        #reshape the resulting array to a distance matrix
        distMat=distMat.reshape(len(x2),len(x1))
        return distMat
        
    else:
        print("Calculation can not be performed because you want to allocate "+str(itemsize/1024./1024.)+" Mb but your Graphics card has only "+str(free/1024./1024.)+" Mb of free memory left")
        return None
