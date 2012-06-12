##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2009 Peter Schmidtke, Daniel Alvarez Garcia, Raik Gruenberg
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


import re,numpy as npy
import sys,os
import math
import fnmatch
from LocalPath import LocalPath
import tools
import Biskit as B
import subprocess

#if we define the environment NOCUDA then by default no cuda and pycua is used
# otherwise we check if we can import the cudm module.
if not os.environ.has_key("NOCUDA"): 
    try :
        import cudm as cu
        _cudaAvail = True
    except ImportError :
        print "cudm was not imported"
        _cudaAvail = False
        pass
else :
    _cudaAvail=False

import time
import Biskit.tools as tools
import geometry as geom


project_path=tools.projectRoot()

lj_coefficients=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"LJPairwiseABCoefficients.pck") #loads LJCoefs dictionary
#lj_coefficients=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"LJPairwiseABCoefficients.pck") #loads LJCoefs dictionary
vdw_parameters=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"vdwParameters.pck") #loads LJCoefs dictionary
gaff_vdw_parameters=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"gaff_vdwParameters.pck") #loads LJCoefs dictionary
vdw_parameters_deNovo=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"vdwParameters_deNovo.pck") #loads LJCoefs dictionary
#lj_atom_types=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"LJAtomTypes.pck")

def mapLJParams2model(model):
    for el in model.AmberTopology_atom_names:
      print el
    model.AmberTopology_vdwRadius=npy.array([float(vdw_parameters[el][0]) for el in model.AmberTopology_atom_names])
    model.AmberTopology_vdwwellDepth=npy.array([float(vdw_parameters[el][1]) for el in model.AmberTopology_atom_names])

def mapGaffLJParams2Mol2model(model):
    model.AmberTopology_vdwRadius=npy.array([float(gaff_vdw_parameters[el][0]) for el in model.atomTypesAmberGaff])
    model.AmberTopology_vdwwellDepth=npy.array([float(gaff_vdw_parameters[el][1]) for el in model.atomTypesAmberGaff])

def getLJcoefs(names1,name2,coef):
    tmp=[str(n+":"+name2) for n in names1]
    
    return npy.array([lj_coefficients[coef][param] for param in tmp],dtype="float64")

def getSubInteractionEnergy(struct, idxrng1,idxrng2,dict=1,eps=78.4):
    """Calculates the interaction energy between sub parts of one single topology
    @ params : 
    @ PDBModel : the pdb structure opened with Biskit and a valid Amber Topology loaded
    @ index range 1 : a list of two elements defining the range of the first residue (or part) in terms of atom indices
    @ index range 2 : a list of two elements defining the range of the second residue (or part) in terms of atom indices
    @ dict (optional) : output energies as dictionary (default) = 1, else if you want a list output, provide 0
    @ eps (optional) : float for dielectric constant, default set to 78.4
    """
    idxs1=npy.arange(idxrng1[0],idxrng1[1]+1)
    idxs2=npy.arange(idxrng2[0],idxrng2[1]+1)
    
    coulomb=npy.zeros([len(idxs1), len(idxs2)])
    vdw=npy.zeros([len(idxs1), len(idxs2)])
    coulomb_dddf=npy.zeros([len(idxs1), len(idxs2)])
    
    #coefs derived from JACS,127,15248-15256, 2005 Agarwal et al
    #eps=78.4 #dielectric constant
    A=-8.5525
    l=0.003627
    k=7.7839
    B=eps-A
    
    tct=0.0
    tc1=0.0
    tc2=0.0
    tvt=0.0
    tv1=0.0
    tv2=0.0
    
    for i, pos2 in enumerate(idxs2):
        
        d=npy.sqrt(((struct.xyz[idxs1]-struct.xyz[pos2])**2).sum(axis=1))
        #print npy.min(d)
        coulomb[:,i]=(struct.AmberCharges[idxs1]*struct.AmberCharges[pos2])/d
        coulomb_dddf[:,i]=(struct.AmberCharges[idxs1]*struct.AmberCharges[pos2])/(d*(A+B/(1.0+k*npy.exp(-l*B*d))))
        #print getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"A12")


        id1 = struct['amber_atype_index'][idxs1] - 1
        id2 = struct['amber_atype_index'][pos2]

        icoIDs=(struct.AmberTopology_n_a_types * id1) + id2 - 1
        coefIDs = struct.AmberTopology_ico[icoIDs] - 1
        Acoef = struct.AmberTopology_a_coefs[coefIDs]
        Bcoef = struct.AmberTopology_b_coefs[coefIDs]

        #ico=struct.AmberTopology_ico[struct.AmberTopology_n_a_types*(struct.AmberTopology_A_types[idxs1-1])+struct.AmberTopology_A_types[pos2]-1]-1
#        acoefs=struct.AmberTopology_a_coefs[ico]
#        bcoefs=struct.AmberTopology_a_coefs[ico]
        
        vdw[:,i]=Acoef/(d**12)-Bcoef/(d**6)

    #GPU PART TODO
#    if len(struct1>2300):
#        xmsize=2300
#    else :
#        xmsize=len(struct1)
#    if len(struct2>)
#    n=len(struct)
#    ncalc=int(math.floor(float(n/msize)))
#    nrest=n-ncalc*msize
#    
#    ttotal=0.0
#    sum=0.0
#    sym=0
#    #for sq matrix
#    for x in range(ncalc):
#        xrnge=range(x*msize,(x+1)*msize)
#        for y in range(ncalc):
#            yrnge=range(y*msize,(y+1)*msize)
#            #print xrnge
#            #print yrnge
#            if x<=y:
#                sym=0
#                if y==x :
#                    sym=1    
#                tmp=struct.AmbertTopology_Excluded_atoms - yrnge[0]
#                t0=time.time()
#                coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],tmp,sym=sym)
#                t1=time.time()
#                ttotal+=t1-t0
#                sum+= npy.sum(coulombMat)
#                #print sum
#            #now get the rest of the tiles
#        if xrnge[-1]<=yrnge[-1]:
#            yrnge =range(yrnge[-1]+1,n)
#            sym=0
#            if xrnge[-1]==yrnge[-1]:
#                sym=1
#            tmp=struct.AmbertTopology_Excluded_atoms - yrnge[0]
#            t0=time.time()
#            coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],tmp,sym=sym)
#            
#            t1=time.time()
#            ttotal+=t1-t0
#            sum+= npy.sum(coulombMat)
#            #print sum
#    
#    print ttotal
    vdwTot=npy.sum(vdw)
    coulombTot=npy.sum(coulomb)
    coulombDDDF=npy.sum(coulomb_dddf)
    if dict==1 :
        return {"EEL":coulombTot,"VDW":vdwTot, "EEL-DDDF":coulombDDDF, "TotalClassic":coulombTot+vdwTot, "TotalwDDDF":coulombDDDF+vdwTot}
    else :
        return [coulombTot,vdwTot,coulombDDDF,coulombTot+vdwTot,coulombDDDF+vdwTot]

def getInteractionEnergyCpu(struct1, struct2, soft = 0, coulombSwitch=1, dict=1, atomsStruct=False,maxCalcSize=3000):
    """
    Get the interaction (pure non bonded) energy between two PDBModels struct1 and struct2
    Calculations are performed on the CPU, consider the GPU functions for big systems, and CPU for small systems
    @arg1: PDBModel: Amber Topology preloaded
    @arg2: PDBModel: Amber Topology preloaded
    @arg3: int (0 or 1); If 1 use VdW soft potentials instead of normal ones.
    @arg4: int(0 or 1): If 1 use Distance Dependant Electrostatic potential calculation instead of classic one.
    @arg5: int(0 or 1): If 1, return a dictionary, else return a list.
    @arg6: False or int(1 or 2): If False, return totals, if 1 return atomic energies for struct1, if 2 return atomic energies for struct2.
    """
    if len(struct2)>len(struct1):
        if atomsStruct:
            if atomsStruct==2: atomsStruct=1
            elif atomsStruct==1:atomsStruct=2
        struct3=struct2
        struct2=struct1
        struct1=struct3
    coulomb=npy.zeros([len(struct1), len(struct2)])
    vdw=npy.zeros([len(struct1), len(struct2)])
    #coulomb_dddf=npy.zeros([len(struct1), len(struct2)])
    #vdwSoft=npy.zeros([len(struct1), len(struct2)])
    
    #coefs derived from JACS,127,15248-15256, 2005 Agarwal et al
    eps=78.4 #dielectric constant
    A=-8.5525
    l=0.003627
    k=7.7839
    B=eps-A
    
    tct=0.0
    tc1=0.0
    tc2=0.0
    tvt=0.0
    tv1=0.0
    tv2=0.0
    if _cudaAvail is True:
        dm=cu.getDistanceMatrix3D(struct1.xyz,struct2.xyz)
    
    for i, pos2 in enumerate(struct2):
        
        if _cudaAvail is True:
            d=dm[i,:]
        else:
            d=npy.sqrt(((struct1.xyz-struct2.xyz[i])**2).sum(axis=1))
        #print npy.min(d)
        if coulombSwitch == 0:  #Normal Coulumb potential
            coulomb[:,i]=(struct1.AmberCharges*struct2.AmberCharges[i])/d
        else:                   #Distance dependant Coulumb potential
            coulomb[:,i]=(struct1.AmberCharges*struct2.AmberCharges[i])/(d*(A+B/(1.0+k*npy.exp(-l*B*d))))
        #print getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"A12")
        if soft == 0:   #Normal VdW potentials
            vdw[:,i]=getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"A12")/(d**12)-getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"B6")/(d**6)

            #cannot use this here as it's based on two distinct topologies
            #ico=struct.ico[struct.AmberTopology_n_a_types*(struct.amber_atype_index+struct.amber_atype_index[pos2])]
            #acoefs=struct.AmberTopology_a_coefs[ico]
            #bcoefs=struct.AmberTopology_a_coefs[ico]


            #ljPot=acoefs/(d**12)-bcoefs/(d**6)
            #np.testing.assert_array_equal(ljPot,vdw[:,i])


        else:   #Soft VdW potentials
            vdw[:,i]=getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"A12")/(d**12)-getLJcoefs(struct1.AmberTopology_atom_names,struct2.AmberTopology_atom_names[i],"B10")/(d**10)
    #GPU PART TODO
#    if len(struct1>2300):
#        xmsize=2300
#    else :
#        xmsize=len(struct1)
#    if len(struct2>)
#    n=len(struct)
#    ncalc=int(math.floor(float(n/msize)))
#    nrest=n-ncalc*msize
#    
#    ttotal=0.0
#    sum=0.0
#    sym=0
#    #for sq matrix
#    for x in range(ncalc):
#        xrnge=range(x*msize,(x+1)*msize)
#        for y in range(ncalc):
#            yrnge=range(y*msize,(y+1)*msize)
#            #print xrnge
#            #print yrnge
#            if x<=y:
#                sym=0
#                if y==x :
#                    sym=1    
#                tmp=struct.AmbertTopology_Excluded_atoms - yrnge[0]
#                t0=time.time()
#                coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],tmp,sym=sym)
#                t1=time.time()
#                ttotal+=t1-t0
#                sum+= npy.sum(coulombMat)
#                #print sum
#            #now get the rest of the tiles
#        if xrnge[-1]<=yrnge[-1]:
#            yrnge =range(yrnge[-1]+1,n)
#            sym=0
#            if xrnge[-1]==yrnge[-1]:
#                sym=1
#            tmp=struct.AmbertTopology_Excluded_atoms - yrnge[0]
#            t0=time.time()
#            coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],tmp,sym=sym)
#            
#            t1=time.time()
#            ttotal+=t1-t0
#            sum+= npy.sum(coulombMat)
#            #print sum
#    
#    print ttotal
    if not atomsStruct:
        totVdw=npy.sum(vdw)
        totCoulomb=npy.sum(coulomb)
    else:
        if atomsStruct == 1: ax=1
        else: ax=0
        totVdw = vdw.sum(axis=ax)
        totCoulomb = coulomb.sum(axis=ax)
    #totEELDDDF=npy.sum(coulomb_dddf)
    #totVdwSoft=npy.sum(vdwSoft)
    if dict==1 :
        #return {"EEL":totCoulomb,"VDW":totVdw, "EEL-DDDF":totEELDDDF,"TotalClassic":totCoulomb+totVdw,"totalWEELDDDF":totEELDDDF+totVdw,"totalWEELDDDFSoft":totVdwSoft+totEELDDDF}
        return {"EEL":totCoulomb,"VDW":totVdw, "Total":totCoulomb+totVdw}
    else :
        #return [totCoulomb,totVdw,totEELDDDF,totCoulomb+totVdw,totEELDDDF+totVdw,totEELDDDF+totVdwSoft]
        return [totCoulomb, totVdw, totCoulomb+totVdw]


def getEnergyOfSystem(struct):
    import cudm as cu
    print """Warning, this function yields potentially wronge energies for now, so don't use it!!!!"""
    msize=5000
    n=len(struct)
    if n<msize :
      msize=n
    print n, "number of atoms"
    ncalc=int(math.floor(float(n/msize)))
    nrest=n-ncalc*msize
    ttotal=0.0
    sum=0.0
    sym=0
    #for sq matrix
    for x in range(ncalc):
        xrnge=range(x*msize,(x+1)*msize)
        for y in range(ncalc):
            yrnge=range(y*msize,(y+1)*msize)
            #print xrnge
            #print yrnge
            if x<=y:
                sym=0
                if y==x :
                    sym=1    
                #tmp=struct.AmberTopology_Excluded_atoms - yrnge[0]
                t0=time.time()
                coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],struct.AmberTopology_Excluded_atoms,sym=sym,ystart=yrnge[0])
                t1=time.time()
                ttotal+=t1-t0
                sum+= npy.sum(coulombMat)
                print coulombMat
                #print sum
            #now get the rest of the tiles
        if xrnge[-1]<=yrnge[-1]:
            yrnge =range(yrnge[-1]+1,n)
            #print yrnge
            #print xrnge
            if yrnge != []:
                sym=0
                if xrnge[-1]==yrnge[-1]:
                    sym=1
                #tmp=struct.AmberTopology_Excluded_atoms - yrnge[0]
                t0=time.time()
                #print yrnge[0]
                print struct.AmberTopology_Excluded_atoms 
                #print struct.AmberTopology_Excluded_idx_start_list[xrnge]
                coulombMat=cu.getCoulombMat(struct.xyz[xrnge],struct.xyz[yrnge],struct.AmberCharges[xrnge],struct.AmberCharges[yrnge],struct.AmberTopology_Excluded_idx_start_list[xrnge],struct.AmberTopology_Excluded_atoms,sym=sym,ystart=yrnge[0])
                print coulombMat
                t1=time.time()
                ttotal+=t1-t0
                sum+= npy.sum(coulombMat)
            #print sum
    
#    print ttotal
    return {"EEL":sum}


def getInteractionForcesCpu(struct1, struct2, soft=0, coulombSwitch=1, dict=1):
    """Get the interaction (pure non bonded) forces between two PDBModels struct1 and struct2
    Calculations are performed on the CPU, consider the GPU functions for big systems, and CPU for small systems
    @ argument 1 : PDBModel : 1st structure (needs preloaded Amber Topology)
    @ argument 2 : PDBModel : 2nd structure (needs preloaded Amber Topology)
    @ argument 3: soft: int( 0 or 1): if soft is 1, returns soft potentials calculations only. If 0 as default returns normal potential calculations results.
    @ argument 4: coulombSwitch: int(0 or 1). If 1 use distance dependant Coulomb potentials. 
    @ optional argument 4 : dict : set this to 1 if you want a dictionary as output, to 0 if you prefer a list, both are ordered the same way.
    """
    if len(struct2)>len(struct1):
        struct3=struct2
        struct2=struct1
        struct1=struct3

    vdwForce=npy.zeros([len(struct1), len(struct2)],dtype="float64")
    coulombForce=npy.zeros([len(struct1), len(struct2)],dtype="float64")

    #coefs derived from JACS,127,15248-15256, 2005 Agarwal et al
    eps=78.4 #dielectric constant
    A=-8.5525
    l=0.003627
    k=7.7839
    B=eps-A
    
    tct=0.0
    tc1=0.0
    tc2=0.0
    tvt=0.0
    tv1=0.0
    tv2=0.0
    vectors=[]
    
    t0=time.time()
    if _cudaAvail is True:
        dm=cu.getDistanceMatrix3D(struct1.xyz,struct2.xyz)    
        
    for i, pos2 in enumerate(struct2):
        if _cudaAvail is True:
            d1 = dm[i,:]
        else:
            d1 = npy.sqrt(((struct1.xyz-struct2.xyz[i])**2).sum(axis=1))

        dmask = d1<5.0
        d = d1[dmask]
        A12=getLJcoefs(struct1.AmberTopology_atom_names[dmask],struct2.AmberTopology_atom_names[i],"A12")
                
        #print npy.min(d)
        #print "type",A12.dtype, B10, B6
        
        if soft == 0:
            B6=getLJcoefs(struct1.AmberTopology_atom_names[dmask],struct2.AmberTopology_atom_names[i],"B6")
            vdwTmp=-(12.0*A12)/d**13 + (6.0*B6)/d**7
        else:
            B10=getLJcoefs(struct1.AmberTopology_atom_names[dmask],struct2.AmberTopology_atom_names[i],"B10")
            vdwTmp=-(12.0*A12)/d**13 + (10.0*B10)/d**11
            
        vdwForce[dmask,i] = vdwTmp
        
        if coulombSwitch == 0:
            elEnergyTmp=(struct1.AmberCharges[dmask]*struct2.AmberCharges[i])/(d**2)
        else:
            #wrongcoulombForce_dddf[:,i]=(struct1.AmberCharges*struct2.AmberCharges[i])/(((B+A)*npy.exp(2.*l*d*B)+(k*l*d*(B**2)+k*B+2.*k*A)*npy.exp(l*d*B)+(k**2)*A)/(npy.exp(2.*l*d*B)+2.*k*npy.exp(l*d*B)+k**2))
            elEnergyTmp=-((struct1.AmberCharges[dmask]*struct2.AmberCharges[i])*((B+A)*npy.exp(2.*l*d*B)+(k*l*d*(B**2)+k*B+2.*k*A)*npy.exp(l*d*B)+(k**2)*A)/(((d**2)*(B**2)+2.*(d**2)*A*B+(d**2)*(A**2))*npy.exp(2.0*l*d*B)+(2.*k*(d**2)*A*B+2*k*(d**2)*(A**2))*npy.exp(l*d*B)+(k**2)*(d**2)*(A**2)))
        
        coulombForce[dmask,i] = elEnergyTmp
        

        vdwTmp.shape=(len(vdwTmp),1)
        elEnergyTmp.shape=(len(elEnergyTmp),1)
        d.shape=(len(d),1)
        totalForces = vdwTmp+elEnergyTmp

        vectors.append(npy.sum(((struct1.xyz[dmask]-struct2.xyz[i])/d)*totalForces, axis=0))
        #print npy.sum(totalForces)
        #vectSoft.append(npy.sum(((struct1.xyz[dmask]-struct2.xyz[i])/d)*totalSoftForces,axis=0))
    
    totVdwForce=npy.sum(vdwForce)
    totCoulombForce=npy.sum(coulombForce)

    t1=time.time()
    #print "needed %.3f s on the CPU"%(t1-t0)
    
    if str(totCoulombForce)=="nan":
      totCoulombForce=0.0
    
    if dict==1 :
        return {"EEL":totCoulombForce, "VDW":totVdwForce, "Total":totCoulombForce+totVdwForce, "Vectors":vectors}
    else :
        return [totCoulombForce, totVdwForce, totCoulombForce+totVdwForce, vectors]


#
#def geom.superpose3D(ref, target, weights=None,refmask=None,targetmask=None,returnRotMat=False):
#    """geom.superpose3D performs 3d superposition using a weighted Kabsch algorithm : http://dx.doi.org/10.1107%2FS0567739476001873 & doi: 10.1529/biophysj.105.066654
#    definition : geom.superpose3D(ref, target, weights,refmask,targetmask)
#
#    @parameter 1 :  ref - xyz coordinates of the reference structure (the ligand for instance)
#    @type 1 :       float64 numpy array (nx3)
#    ---
#    @parameter 2 :  target - theoretical target positions to which we should move (does not need to be physically relevant.
#    @type 2 :       float 64 numpy array (nx3)
#    ---
#    @parameter 3:   weights - numpy array of atom weights (usuallly between 0 and 1)
#    @type 3 :       float 64 numpy array (n)
#
#    @parameter 4:   mask - a numpy boolean mask for designating atoms to include
#    Note ref and target positions must have the same dimensions -> n*3 numpy arrays where n is the number of points (or atoms)
#    Returns a set of new coordinates, aligned to the target state as well as the rmsd
#    """
#
#    if weights == None :
#        weights=1.0
#    if refmask == None :
#        refmask=npy.ones(len(ref),"bool")
#    if targetmask == None :
#        targetmask=npy.ones(len(target),"bool")
#
#    #first get the centroid of both states
#    ref_centroid = npy.mean(ref[refmask]*weights,axis=0)
#    #print ref_centroid
#    refCenteredCoords=ref-ref_centroid
#    #print refCenteredCoords
#
#    target_centroid=npy.mean(target[targetmask]*weights,axis=0)
#    targetCenteredCoords=target-target_centroid
#    #print targetCenteredCoords
#
#    #the following steps come from : http://www.pymolwiki.org/index.php/OptAlign#The_Code and http://en.wikipedia.org/wiki/Kabsch_algorithm
#
#    # Initial residual, see Kabsch.
#    E0 = npy.sum( npy.sum(refCenteredCoords[refmask] * refCenteredCoords[refmask]*weights,axis=0),axis=0) + npy.sum( npy.sum(targetCenteredCoords[targetmask] * targetCenteredCoords[targetmask]*weights,axis=0),axis=0)
#    reftmp=npy.copy(refCenteredCoords[refmask])
#    targettmp=npy.copy(targetCenteredCoords[targetmask])
#
#    #print refCenteredCoords[refmask]
#    #single value decomposition of the dotProduct of both position vectors
#    try:
#        dotProd = npy.dot( npy.transpose(reftmp), targettmp* weights)
#
#        V, S, Wt = npy.linalg.svd(dotProd )
#
#    except Exception:
#        try:
#            dotProd = npy.dot( npy.transpose(reftmp), targettmp)
#            V, S, Wt = npy.linalg.svd(dotProd )
#        except Exception:
#            print >> sys.stderr,"Couldn't perform the Single Value Decomposition, skipping alignment"
#
#        return ref, 0
#
#    # we already have our solution, in the results from SVD.
#    # we just need to check for reflections and then produce
#    # the rotation.  V and Wt are orthonormal, so their det's
#    # are +/-1.
#    reflect = float(str(float(npy.linalg.det(V) * npy.linalg.det(Wt))))
#
#    if reflect == -1.0:
#        S[-1] = -S[-1]
#        V[:,-1] = -V[:,-1]
#
#    rmsd = E0 - (2.0 * sum(S))
#    rmsd = npy.sqrt(abs(rmsd / len(ref[refmask])))   #get the rmsd
#
#    #U is simply V*Wt
#    U = npy.dot(V, Wt)  #get the rotation matrix
#
#    # rotate and translate the molecule
#    new_coords = npy.dot((refCenteredCoords), U)+ target_centroid  #translate & rotate
#    #new_coords=(refCenteredCoords + target_centroid)
#    #print U
#    if returnRotMat :
#        return new_coords,rmsd, U
#    return new_coords,rmsd
#

def minimize(protein, orig_ligand, soft=0, stepSize=0.1, energyCut=0., maxSteps=1000, verbose=0):
    """
    Steepest descent minimization algorithm.
    definition : minimize(protein, ligand, soft=0, stepSize=0.1, energyCut=0., maxSteps=1000) 
    
    @parameter 1 :  protein - protein (or fixed structure)
    @type 1 :       PDBModel
    ---
    @parameter 2 :  ligand - PDBModel of the ligand (or mobile structure (to minimize))
    @type 2 :       PDBModel
    ---
    @parameter 3 :  (default 0) : soft - Flag for using soft potentials, if 1 yes, if 0 no. Will use distance dependent Electrostatics Potential by default.
    @type 3 :       int
    ---
    @parameter 4 :  (default 0.1) : stepSize - the displacement in angstroms following the force vectors
    @type 4 :       float
    ---
    @parameter 5 :  (default 0.0) : energyCut - stop the minimization when energyCut is reached
    @type 5 :       float
     ---
    @parameter 6 :  (default 1000) : maxSteps - stop the minimization when reached maximum steps even not converged
    @type 6 :       int
     ---
    @parameter 7 :  (default 0) : verbose : flag to know if output of intermediate energies (1) or not (0)
    @type 7 :       int
    """
    ligand=orig_ligand.clone()
    ligand.AmberCharges = orig_ligand.AmberCharges
    ligand.AmberTopology_atom_names = orig_ligand.AmberTopology_atom_names
    xyzBackUp = ligand.xyz
    i = 0
    energy = 1000.
    E1 = energy
    convergenceTrack = []
    convergenceId = 0
    if verbose:
        print "Starting minimisation"
    while energy > energyCut and i<maxSteps:
        
        x = getInteractionForcesCpu(protein, ligand, soft=soft, dict=1)
        forceVectors = x['Vectors']
        
        target_position = ligand.xyz + stepSize * ( npy.array(forceVectors)  / npy.sqrt((npy.array(forceVectors)**2).sum(axis=1)).reshape(len(forceVectors), 1 ) )
                
        #Calculate weights upon the forces. The bigger the force over one atom, the bigger the weight
        totals = npy.sqrt((npy.array(forceVectors)**2).sum(axis=1))
        weights = (totals / totals.max()).reshape(len(totals), 1)
        (newCoords,rmsd) = geom.superpose3D(xyzBackUp,target_position, weights)
        
        #Update coordinates
        xyzBackUp = newCoords
        ligand.xyz = newCoords
        
        #Check energies
        #ligand.writePdb(tools.testRoot()+os.sep+"biskitFF"+os.sep+"kk_normal"+str(i)+".pdb")
        energy = getInteractionEnergyCpu(protein, ligand, soft= soft, dict=1)['Total']
        if verbose:
            print("""Total Energy : %.2f kcal/mol
RMSD to previous pose : %.2f Angstroms
Current Stepsize for molecular displacement : %.2f"""%(energy,rmsd,stepSize))
            print("-"*10)
        E2 = energy
        
        #StepSize modulation
        if E2 < E1: stepSize += 0.20*stepSize
        elif E2 > E1: stepSize -= 0.5*stepSize
        
        #Convergence tracking
        convergenceTrack.append(E2)
        convergenceId += 1
        
        #Check if it is converged calculating standard deviation
        if convergenceId == 20:
            convergenceTrack = npy.array(convergenceTrack)
            deviation = npy.std(convergenceTrack)
            mean=npy.mean(convergenceTrack)
            if abs(deviation/mean) < 0.05 and mean > 10.: 
                break
            elif abs(deviation/mean) < 0.001:
                break
            convergenceTrack = []
            convergenceId = 0
               
        E1 = E2    
        i+=1

    ligand.xyz = newCoords
    if verbose:
        print "Minimisation finished"
    #ligand.writePdb(tools.testRoot()+os.sep+"biskitFF"+os.sep+"kk_final.pdb")
    return (ligand, energy)


def getEquilibriumDistance(PDBModel, atom1ID, atom2ID, verbose=False):
    """
    This function returns the equilibrium distance between 
    two atoms of a PDBModel using the topology information.
    Equilibrium distance can be obtained form Lennard Jones
    coefficients given a pairwise interaction with the expression:
        r = (2 * (A/B)) ** (1/6)
    It involves a 6th root.
    
    parm1   PDBMOdel    System with Amber Topology loaded
    type                Biskit.PDBModel
    
    parm2   atom1ID     index of the first atom || name of the atomtype
    type                int || string
    
    parm3   atom2ID     index of the second atom || name of the atomtype
    type                int || string
    """
    if type(atom1ID) is int:
        id1 = PDBModel[atom1ID]['amber_atype_index'] - 1
    else:
        tmp = npy.where(PDBModel['amber_atype'] == atom1ID)[0][0]
        id1 = PDBModel[tmp]['amber_atype_index'] - 1
    
    if type(atom2ID) is int:
        id2 = PDBModel[atom2ID]['amber_atype_index']
    else:
        tmp = npy.where(PDBModel['amber_atype'] == atom2ID)[0][0]  
        id2 = PDBModel[tmp]['amber_atype_index']
    
    icoID = (PDBModel.AmberTopology_n_a_types * id1) + id2 - 1
    coefID = PDBModel.AmberTopology_ico[icoID] - 1
    A = PDBModel.AmberTopology_a_coefs[coefID]
    B = PDBModel.AmberTopology_b_coefs[coefID]
    
    if verbose:
        return {"A":A, "B":B, "distance":(2*(A/B))**(1/6.)}
    else:
        return (2*(A/B))**(1/6.)

def equilibriumMatrix(PDBModel, buff=0., keepPolars=True):
    """
    This function returns a ditionary with all the possible 
    equilibrium distances between the different atom types 
    in the PDBModel.
    buff: Increase equilibrium distance by this buffer
    keepPolars = prevent the increasing of buff in polar-polar 
    interactions.
    """
    eqDict = {}
    difTypes = npy.unique(PDBModel['amber_atype'])
    polarList = []
    for an in ['O','N']:
        [polarList.append(atn) for atn in fnmatch.filter(difTypes, '%s*'%(an))]
    for t in difTypes:
        for m in difTypes:
            eqDict[t+m] = getEquilibriumDistance(PDBModel, t, m) + buff
            if (keepPolars is True and (t in polarList or m in polarList)): eqDict[t+m] -= buff
    
    return eqDict

#############
##  TESTING        
##############
import Biskit.test as BT
#
class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.LONG ]
#
    def prepare( self ):
#        self.f = os.path.join( tools.testRoot(), 'eedensity/3tgi.ezd.gz' )
        print("Loading protein from pdb file :")
        self.protein=B.PDBModel(tools.testRoot()+os.sep+"biskitFF"+os.sep+"protein.pdb")
        print("Loading ligand from pdb file :")
        self.ligand=B.PDBModel(tools.testRoot()+os.sep+"biskitFF"+os.sep+"rotation5.pdb")

        print("Reading Amber Topology for protein :")
        self.protein.loadAmberTopology(tools.testRoot()+os.sep+"biskitFF"+os.sep+"protein.top")
        print("Reading Amber Topology for ligand :")
        self.ligand.loadAmberTopology(tools.testRoot()+os.sep+"biskitFF"+os.sep+"ligand.top")

        

#        print energyNormal

#        minimLigandSoft.writePdb(tools.testRoot()+os.sep+"biskitFF"+os.sep+'minimizedSoft.pdb')
#        minimLigandNormal.writePdb(tools.testRoot()+os.sep+"biskitFF"+os.sep+'minimizedNormal.pdb')

#        self.p = None
#
    def test_LJ_parameters(self):
        """Testing reading pre-calculated Lennard Jones parameters for atom pairs from pickled files"""
        lj_coefficients=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"LJPairwiseABCoefficients.pck") #loads LJCoefs dictionary

        self.assert_(lj_coefficients["B6"]["CZ:OS"]-519.16333070072085<1e-8)
        self.assert_(lj_coefficients["B6"]["IM:S"]-2522.5811073172445<1e-8)

    def test_VDW_parameters(self):
        """Testing reading pickled van der Waals parameters for each Amber Atom Type"""
        vdw_parameters=tools.load(project_path+os.sep+"Biskit"+os.sep+"data"+os.sep+"amber"+os.sep+"leap"+os.sep+"vdwParameters.pck") #loads LJCoefs dictionary
        self.assert_(npy.all( vdw_parameters['HZ']==[ 1.459,  0.015,  0.   ,  1.   ]))

    def test_minimizationSoft( self ):
        """BiskitFF minimsation testing routines on soft potentials"""
        print "Testing Soft Potential Minimization routines:"
        print "="*20
        self.minimLigandSoft, energySoft = minimize(self.protein, self.ligand, soft=1,verbose=1)
        print "Soft potential minimisation --> energy : ",energySoft
        self.assert_((energySoft+1.38942145266)<1e-3)

    def test_minimizationNormal(self):
        print "Testing Normal Minimization routines:"
        #best practices: first minimise using a soft potential, then take the minimised pose and use normal potentials
        print "="*20
        minimLigandSoft, energySoft = minimize(self.protein, self.ligand, soft=1,verbose=1)
        self.assert_((energySoft+1.38942145266)<1e-3)
        minimLigandNormal, energyNormal = minimize(self.protein, minimLigandSoft, soft=0, stepSize=0.1, energyCut=-20.5, maxSteps=500,verbose=1)
        print "Normal potential minimisation --> energy : ",energyNormal
        self.assert_((energyNormal+12.0369455356)<1e-3)

    def test_getSubInteractionPotential(self):
        print "Testing the generic interaction potential routine"
        print "="*20
        #get the interaction energy between residue 1 and 30 of the protein
        #CAUTION: we currently don't retain excluded atom lists normally for the interaction energy calculation!!
        idxs1=npy.where(self.protein.maskFrom("residue_number",1))[0]
        idxs2=npy.where(self.protein.maskFrom("residue_number",30))[0]
        
        print getSubInteractionEnergy(self.protein,[idxs1[0],idxs1[-1]],[idxs2[0],idxs2[-1]])
        
#

#        if self.local:
#            print "parsing ...",
#        self.p = self.p or EZDParser( self.f )
#
#        if self.local:
#            print "Done"
#        
#        self.origin = self.p.getCartesianOrigin()
#        target = [ -38.92757372, -115.8054999 ,  -36.27633333]
#
#        ## compare result with expected result but allow for numeric deviations
#        self.assert_( npy.sum(self.origin - target) < 1e-8 )
#
#        self.grid = self.p.getCartesianPositionGrid()
#        self.assertEquals( npy.shape( self.grid ), (129, 105, 104, 4) )
#
#        self.i = self.p.getIdxFromCart( [10, 10, 10] )
#        self.assert_( npy.all( self.i == [  82.,  212.,   80.] ) )
#        


if __name__ == '__main__':

    BT.localTest()
