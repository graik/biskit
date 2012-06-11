##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2009-10-21 14:51:08 +0200 (Wed, 21 Oct 2009) $
## $Revision: 936 $

"""
Store and manipulate coordinates and atom information.
"""

import tools as T
import molUtils
import mathUtils
import match2seq
import rmsFit
from LocalPath import LocalPath
from Errors import BiskitError
from Biskit import EHandler
from ProfileCollection import ProfileCollection, ProfileError
from PDBParserFactory import PDBParserFactory
from PDBParseFile import PDBParseFile
import Biskit as B

import numpy.oldnumeric as N
import numpy as npy
import numpy.oldnumeric.mlab as MLab
import os, sys
import copy
import time
import string
import types
import Scientific.IO.PDB as IO

class MOL2Collection:
    """Stores mutltiple MOL2Models"""
    def __init__(self):
        self.mols=[]
    
    def __init__(self, source):
        self.mols=[]
        if type( source ) is str and os.path.isfile( source ) :
            self.mols=self.readMol2Collection(source)
            self.source=source
        else : 
            print "file : "+source+" not found or not a valid input file!"

    def __repr__(self):
        return "Mol2 Collection containing %d molecules"%(len(self.mols))
        
    def readMol2Collection(self,source):
        mols=[]
        f=open(source,"r")
        while not f.closed:
            mols.append(MOL2Model(f))
        mols=mols[:-1]
        return mols
        
class MOL2Model:
    """
    Store and manipulate coordinates and atom infos stemming from a PDB file.
    Coordinates are stored in the numpy array 'xyz'; the additional atom infos
    from the PDB (name, residue_name, and many more) are efficiently stored in
    a L{PDBProfiles} instance 'atoms' which can be used to also associate
    arbitrary other data to the atoms. Moreover, a similar collection
    'residues' can hold data associated to residues (but is initially empty).
    A normal dictionary 'info' accepts any information about the whole model.

    For detailed documentation,
    see http://biskit.pasteur.fr/doc/handling_structures/PDBModel

    @todo:
       * outsource validSource into PDBParserFactory
       * prevent repeated loading of test PDB for each test
    """

    #: keys of all atom profiles that are read directly from the PDB file
    PDB_KEYS = ['name', 'residue_number', 'insertion_code', 'alternate',
                'name_original', 'chain_id', 'occupancy', 'element',
                'segment_id', 'charge', 'residue_name', 'after_ter',
                'serial_number', 'type', 'temperature_factor']

    def __init__( self, fhandle):
        if type(fhandle) is str :
            fhandle=open(fhandle,"r")
            
        if type(fhandle) is file :
            f=fhandle
            
            tmp=f.readline()
            while tmp.strip()!="@<TRIPOS>MOLECULE" and tmp:
                tmp=f.readline()
                
            f.readline()
            try:
                (self.nAtoms,self.nBonds)=[int(el) for el in f.readline().split()[:2]]
            except ValueError:
                f.close()
                return
            self.source=fhandle.name
            """@<TRIPOS>ATOM
               1 C       43.9700    33.1090    31.4590 C.ar  1 <1>  -0.1440 
            """
            tmp=f.readline()
            while tmp.strip()!="@<TRIPOS>ATOM" and tmp:
                tmp=f.readline()
            self.xyz=npy.zeros([self.nAtoms,3],dtype="float64")
            self.atomNames=npy.zeros([self.nAtoms],dtype="S2")
            self.atomTypesTripos=npy.zeros([self.nAtoms],dtype="S5")
            self.charges=npy.zeros([self.nAtoms],dtype="float32")
            
            for nat in xrange(self.nAtoms):
                line=f.readline().strip()
                line=line.split()
                #print line[1]
                try :
                    self.atomNames[nat]=npy.str(line[1])
                    self.xyz[nat,:]=npy.array([line[2:5]],dtype="float64")
                    self.atomTypesTripos[nat]=line[5]
                    self.charges[nat]=float(line[-1])
                except Exception,error:
                    print "error during parsing of mol2 line : "
                    print line
                    print error
                    break
            while tmp.strip()!="@<TRIPOS>BOND" and tmp:
                tmp=f.readline()
            self.bonds=npy.zeros([self.nBonds,2],dtype="int")
            self.bondTypes=npy.zeros([self.nBonds],dtype="S3")
            
            for nb in xrange(self.nBonds):
                line=f.readline().strip()
                line=line.split()
                try:
                    self.bonds[nb,:]=line[1:3]
                    self.bondTypes[nb]=line[3]
                except Exception,error:
                    print "error during parsing of mol2 line : "
                    print line
                    print error
                    break
            
            while tmp.strip()!="@<TRIPOS>SUBSTRUCTURE" and tmp:
                tmp=f.readline()
            tmp=f.readline()
            try :
                self.molName=tmp.split()[1]
            except Exception,err:
                self.molName="UNKOWN"
            
            self.atomTypesAmberGaff=npy.zeros([self.nAtoms],dtype="S2")
            
            for i,atype in enumerate(self.atomTypesTripos):
                if atype[:2]=="O.":
                    #print i, atype
                    asel1= self.bonds[(self.bonds[:,0]==(i+1)),1]
                    asel2= self.bonds[(self.bonds[:,1]==(i+1)),0]
                    
                    for a in asel1:
                        if self.atomTypesTripos[a-1]=="H":
                            self.atomTypesAmberGaff[i]="OH"
                            break
                        elif self.atomTypesTripos[a-1][:2]=="S.":
                            self.atomTypesAmberGaff[i]="OS"
                            break
                    if len(self.atomTypesAmberGaff[i])<1:
                        for a in asel2:
                            if self.atomTypesTripos[a-1]=="H":
                                self.atomTypesAmberGaff[i]="OH"
                                break
                            elif self.atomTypesTripos[a-1][:2]=="S.":
                                self.atomTypesAmberGaff[i]="OS"
                                break
                    if len(self.atomTypesAmberGaff[i])<1:
                        self.atomTypesAmberGaff[i]="O"
                    #for triposType in self.atomTypesTripos[]:
                        #print triposType
                        
                elif  atype[:1]=="H":
                    asel1= self.bonds[(self.bonds[:,0]==(i+1)),1]
                    asel2= self.bonds[(self.bonds[:,1]==(i+1)),0]
                    
                    for a in asel1:
                        if self.atomTypesTripos[a-1][:2]=="N.":
                            self.atomTypesAmberGaff[i]="HN"
                            break
                        elif self.atomTypesTripos[a-1][:2]=="S.":
                            self.atomTypesAmberGaff[i]="HS"
                            break
                        elif self.atomTypesTripos[a-1][0]=="P":
                            self.atomTypesAmberGaff[i]="HP"
                            break
                        elif self.atomTypesTripos[a-1][:2]=="O.":
                            self.atomTypesAmberGaff[i]="HO"
                            break
                    if len(self.atomTypesAmberGaff[i])<1:
                        for a in asel2:
                            if self.atomTypesTripos[a-1][:2]=="N.":
                                self.atomTypesAmberGaff[i]="HN"
                                break
                            elif self.atomTypesTripos[a-1][:2]=="S.":
                                self.atomTypesAmberGaff[i]="HS"
                                break
                            elif self.atomTypesTripos[a-1][0]=="P":
                                self.atomTypesAmberGaff[i]="HP"
                                break
                            elif self.atomTypesTripos[a-1][:2]=="O.":
                                self.atomTypesAmberGaff[i]="HO"
                                break
                    if len(self.atomTypesAmberGaff[i])<1:
                        self.atomTypesAmberGaff[i]="H"    
                       
                elif atype[:2]=="C.":
                    self.atomTypesAmberGaff[i]="C"
                elif atype[:2]=="N.":
                    self.atomTypesAmberGaff[i]="N"
                elif atype[:2]=="S.":
                    self.atomTypesAmberGaff[i]="S"
                elif atype[0]=="P":
                    self.atomTypesAmberGaff[i]="P"
                elif atype[0]=="F":
                    self.atomTypesAmberGaff[i]="F"
                elif atype[:2]=="Cl":
                    self.atomTypesAmberGaff[i]="CL"
                elif atype[:2]=="Br":
                    self.atomTypesAmberGaff[i]="BR"
                elif atype[0]=="I":
                    self.atomTypesAmberGaff[i]="I"
        else :
            if fhandle==None:
                self.atomTypesAmberGaff=None
                self.atomTypesTripos=None
                self.bonds=None
                self.bondTypes=None
                self.xyz=None
                self.atomNames=None
        
        self.setRings()
                
    def copy(self):
        new=self.__init__()
        #TODO continue here
        return new

    def __repr__(self):
        return "Mol2Model with %d atoms and %d bonds"%(len(self.xyz),len(self.bonds))
    
    def setRings(self):
        """set rings found in the molecule to the self.rings table. rings will contain per table element a numpy array with indices to the atoms arrays (0 first atom). So to get the positions of the first ring : 
        You can simply type self.xyz[self.rings[0],:]
      """
        self.rings=[]
        u_atoms =npy.unique(self.bonds[:, 0])
        for atomId in u_atoms:
            idxs=npy.where(self.bonds[:, 1]==atomId)[0]
            if len(idxs)>1:
                #print atomId,self.bonds[idxs,0]
                self.rings.append(npy.arange(self.bonds[idxs,0][0]-1,atomId))
                
    
    def writeMol2(self,fname):
        f=open(fname,"w")
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("Molecule\n")
        f.write("%d %d %d %d %d\n"%(self.nAtoms, self.nBonds,1,0,0))
        f.write("SMALL\nUSER_CHARGES\n\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i,atom in enumerate(self.atomNames):
            f.write("%3d %-2s      %4.4f    %4.4f    %4.4f %-5s 1 <1>  %+3.4f\n"%(i+1,atom,self.xyz[i,0],self.xyz[i,1],self.xyz[i,2],self.atomTypesTripos[i],self.charges[i]))
        
        f.write("@<TRIPOS>BOND\n")
        for i, bond in enumerate(self.bonds):
            f.write("%3d %3d %3d  %s\n"%(i+1,bond[0],bond[1],self.bondTypes[i]))
        
        f.close()