from quaternion import rotquat, rotmat, qmult
from emath import normalized, norm,dot, vectorangle, vectorcos
import numpy as N


def planarize(model = None, orig_plane =[[1,0,0],[0,0,1]],orig_perp = False, target_plane=[[1,0,0],[0,0,-1]],target_perp = True,):
    
    ## We want perpendicular plane-defining vectors
    
    if orig_perp == False:
        aux = N.cross(orig_plane[0],orig_plane[1])
        orig_plane[1] = N.cross(orig_plane[0],aux)

    if target_perp == False:
        aux = N.cross(target_plane[0],target_plane[1])
        target_plane[1] = N.cross(target_plane[0],aux)
    
    
    v1 = N.array([N.array(normalized(orig_plane[0])),N.array([0,0,1]),N.array([0,0,1])])
    v2 = N.array([N.array(normalized(target_plane[0])),N.array([0,0,1]),N.array([0,0,1])])
    q  = rotquat(v1,v2)
    R = N.transpose(N.matrix(rotmat(q[0])))
    
    orig_plane[1] = orig_plane[1] * R
    
    v1 = N.array([N.array(normalized(orig_plane[1])),N.array([0,0,1]),N.array([0,0,1])])
    v2 = N.array([N.array(normalized(target_plane[1])),N.array([0,0,1]),N.array([0,0,1])])
    q  = rotquat(v1,v2)

    R = R * N.transpose(N.matrix(rotmat(q[0])))
    
    for i in range(len(model.xyz)):
        model.xyz[i] = model.xyz[i] * R
   
    
    return model , R
    
def calcPlane (residue):
    
    sel_atoms = N.where(residue.maskFrom( 'name', ['N','CA','C'] ))
    Npos = sel_atoms[0][0]
    Capos= sel_atoms[0][1]
    Cpos = sel_atoms[0][2]
    
    Ni  = residue.xyz[Npos]
    Ca = residue.xyz[Capos]
    C  = residue.xyz[Cpos]
    
    return  C-Ni,Ca-Ni

def doAAReorientation (residue = None, vectors = []):
    """
    Reorients a residue along the C-Ca axis and C-N axis.
    It's used for having an standard base orientation for all the residues.
    
    @param residue: Residue to be reoriented.
    @type residue: PDBModel
    @param iNC: Atom indexes of the N, CA,C atoms.
    @type iNC: list of int
    
    @return: The residue reoriented.
    @rtype: PDBModel
    """
    sel_atoms = N.where(residue.maskFrom( 'name', ['N','CA','C'] ))
    Npos = sel_atoms[0][0]
    Capos= sel_atoms[0][1]
    Cpos = sel_atoms[0][2]
    
    Ni  = residue.xyz[Npos]
    Ca = residue.xyz[Capos]
    C  = residue.xyz[Cpos]
    
    ## Center
    residue.xyz = residue.xyz - Ni
    
    ## Rotation
    
    nc = C-Ni
    nca = Ca-Ni
    
    residue, R = planarize(residue,N.array([nc,nca]))
    
    resvectors = []
    for v in vectors:
        resvectors.append(v*R)

    return residue, R, resvectors

def flip (model = None, orig_plane =[[1,0,0],[0,0,1]],vectors = []):
    target_plane = [[orig_plane[0][0],orig_plane[0][1],orig_plane[0][2]],[-orig_plane[1][0],-orig_plane[1][1],-orig_plane[1][2]-0.5]]
    
    model, R = planarize(model,N.array(orig_plane),False,N.array(target_plane),False)
    
    
    resvectors = []
    for v in vectors:
        resvectors.append(v*R)

    return model, resvectors
    
#~ def mirrorH
#~ def mirrorF
def flipAA(residue,vectors = []):
    sel_atoms = N.where(residue.maskFrom( 'name', ['N','CA','C'] ))
    Npos = sel_atoms[0][0]
    Capos= sel_atoms[0][1]
    Cpos = sel_atoms[0][2]
    
    Ni  = residue.xyz[Npos]
    Ca = residue.xyz[Capos]
    C  = residue.xyz[Cpos]
    
    residue.xyz = residue.xyz - Ni
    
    nc = C-Ni
    nca = Ca-Ni
    
    return  flip(residue,[nc,nca],vectors)



def stickAAs( a,b , flip ):
    iC= where(a.maskFrom( 'name', ['C'] ))
    #~ print a
    #~ print b
    
    if flip :
        
        b = self.rotateAA(b,[1.,0.,0.],[0.,0.,1.],True,'lol.pdb')
        b = self.rotateAA(b,[0.,0.,1.],[-1.,0.,0.],True,'lol.pdb')
        
        T = array([0,0.56,1.27])
    
    else:
        T = array([0,-0.56,-1.27])
        
    C = a.xyz[iC[0][-1]]
    
    for i in range(len(b.xyz)):
        b.xyz[i] = b.xyz[i] + T + C
    b.update()
    
    
    m = a.concat(b)
    #~ print m
    m.renumberResidues()
    
    m['serial_number'] = range(1,m.lenAtoms()+1)
    m['chain_id'] = m.lenAtoms() * ['A']
    
    return m
    
def concatenateAAChains (a,b) :
    # temptative function for chain concatenation
    # chains may not have TER o OXT
    
    #reorientation of the other molecule in the or. of the other (Nn->Cn to b's N1->C1)
    iCa= where(a.maskFrom( 'name', ['C'] ))
    iNa= where(a.maskFrom( 'name', ['N'] ))
    
    iCb= where(b.maskFrom( 'name', ['C'] ))
    iNb= where(b.maskFrom( 'name', ['N'] ))
    
    #Residues are supposed to have the correct backbone at least
    aCn = a.xyz[iCa[0][-1]] # Cn
    aNn = a.xyz[iNa[0][-1]] # Nn
    
    aCnNn = emath.normalized(aCn - aNn)
    
    bC1 = b.xyz[iCb[0][0]] # C1
    bN1 = b.xyz[iNb[0][0]] # N1
    
    bC1N1 = emath.normalized(bC1 - bN1)
    
    #center
    a.xyz = a.xyz - aNn
    b.xyz = b.xyz - bN1
    #~ print "centered-",bN1,bC1,iCb[0][0],iCa[0][0]
    #Reorient b 
    b = self.rotateAA(b,bC1N1,aCnNn)
    
    c =self.stickAAs(a,b,False)
    #~ print a 
    #~ print b
    #~ print c
    return c
    


##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit.PDBModel import PDBModel
import os

class Test(BT.BiskitTest):
    """ Test cases for Polyfret"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    def test_vectors(self):
        """test some vector operations"""
        import random
        a = N.array([random.randrange(5),random.randrange(5),random.randrange(5) ])
        b = N.array([random.randrange(5),random.randrange(5),random.randrange(5) ])
        c = N.array([random.randrange(5),random.randrange(5),random.randrange(5) ])
           
        ab = a-b
        cb = b-c 
        
        aux = N.cross(ab,cb)
        cc2 = N.cross(ab,aux)
        
        self.assertEqual(N.dot(ab,cc2),0)
        
    def test_Planarize(self):
        """Planarize test case"""
        p = PDBModel(T.testRoot()+"/polysys/allala.pdb")
        residues = p.resModels()
        f = open(T.testRoot()+"/polysys/tots.pdb","w")
        
        for i in range(len(residues)):
            print i
            residues[i].xyz = residues[i].xyz - residues[i].xyz[0]
            nc = residues[i].xyz[2]-residues[i].xyz[0]
            nca = residues[i].xyz[1]-residues[i].xyz[0]
            planarize(residues[i],N.array([nc,nca]))[0].writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
        
        for i in range(len(residues)):
            print i
            doAAReorientation(residues[i])[0].writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
            
        for i in range(len(residues)):
            print i
            res = residues[i].clone()
            flipAA(residues[i])[0].concat(flipAA(flipAA(res)[0])[0]).writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
            
        f.close()
        
        
    def test_CalcPlane(self):
        """CalcPlane test case"""
        
        print "TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"

if __name__ == '__main__':
    BT.localTest()    
        