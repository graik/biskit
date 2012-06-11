#import pickle as pck
import numpy as npy
import Biskit as B
import string,sys
"""script to create a pickled dictionary with A & B coefficients for calculation of the Lennard Jones Potential for all atom types versus all atom types"""

f=open("parm99.dat","r")
parms={}
f.readline()
for line in f.readlines():
  content=line.split()
  parms[content[0]]=npy.array([float(el) for el in content[1:5]],dtype="float64")
f.close()

f=open("gaff.dat","r")
gaffparms={}
f.readline()
for line in f.readlines():
  content=line.split()
  gaffparms[content[0]]=npy.array([float(el) for el in content[1:3]],dtype="float64")
f.close()

#make the following correspondances :
"""N =   NA  N2  N*  NC  NB  NT  NY
C* = CA  CB  CC  CD  CK  CM  CN  CQ  CR  CV  CW  CY  CZ
"""
for el in ["NA",  "N2",  "N*",  "NC",  "NB",  "NT",  "NY"] :
    parms[el]=parms["N"]

for el in ["CA",  "CB", "CC",  "CD",  "CK" , "CM",  "CN",  "CQ",  "CR",  "CV",  "CW",  "CY",  "CZ"]:
    parms[el]=parms["C*"]
    
LJCoefs={"A12":{}, "B6":{},"B10":{},"A13":{},"B7":{},"B11":{}}
gaffLJCoefs={"A12":{}, "B6":{},"B10":{},"A13":{},"B7":{},"B11":{}}
#calculate A and B coefficients for all pairwise atom type interactions
for atomtype1 in parms.keys():
    for atomtype2 in parms.keys():
#        print atomtype1
#            print atomtype2
        LJCoefs["A12"][atomtype1+":"+atomtype2]=npy.sqrt(parms[atomtype1][1]*parms[atomtype2][1])*(parms[atomtype1][0]+parms[atomtype2][0])**12
        LJCoefs["B10"][atomtype1+":"+atomtype2]=2.0*npy.sqrt(parms[atomtype1][1]*parms[atomtype2][1])*(parms[atomtype1][0]+parms[atomtype2][0])**10
        LJCoefs["B6"][atomtype1+":"+atomtype2]=2.0*npy.sqrt(parms[atomtype1][1]*parms[atomtype2][1])*(parms[atomtype1][0]+parms[atomtype2][0])**6
        
B.tools.dump(parms,"vdwParameters.pck")
B.tools.dump(LJCoefs,"LJPairwiseABCoefficients.pck")
B.tools.dump(parms.keys(),"LJAtomTypes.pck")


for atomtype1 in gaffparms.keys():
    for atomtype2 in gaffparms.keys():
        gaffLJCoefs["A12"][atomtype1+":"+atomtype2]=npy.sqrt(gaffparms[atomtype1][1]*gaffparms[atomtype2][1])*(gaffparms[atomtype1][0]+gaffparms[atomtype2][0])**12
        gaffLJCoefs["B10"][atomtype1+":"+atomtype2]=2.0*npy.sqrt(gaffparms[atomtype1][1]*gaffparms[atomtype2][1])*(gaffparms[atomtype1][0]+gaffparms[atomtype2][0])**10
        gaffLJCoefs["B6"][atomtype1+":"+atomtype2]=2.0*npy.sqrt(gaffparms[atomtype1][1]*gaffparms[atomtype2][1])*(gaffparms[atomtype1][0]+gaffparms[atomtype2][0])**6
        
B.tools.dump(gaffparms,"gaff_vdwParameters.pck")
B.tools.dump(gaffLJCoefs,"gaff_LJPairwiseABCoefficients.pck")
B.tools.dump(gaffparms.keys(),"gaff_LJAtomTypes.pck")

