from quaternion import rotquat, rotmat, qmult,quaternion
import numpy as N
from vectors import normalized, norm,dot, angle, vcos, cross,sub,matrix2list,vcos
from tools import pad , cutdepth
import Biskit.molUtils as MU 

def orientPlanes(modela = None, orig_plane =[[1,0,0],[0,0,1]],orig_perp = False, target_plane=[[1,0,0],[0,0,-1]],target_perp = True):
    
    ## We want perpendicular plane-defining vectors
    model = modela.clone()
    
    if orig_perp == False:
        aux = cross(orig_plane[0],orig_plane[1])
        orig_plane[1] = cross(aux,orig_plane[0])
    
    if target_perp == False:
        aux = cross(target_plane[0],target_plane[1])
        target_plane[1] = N.cross(aux,target_plane[0])
       
    v1 = N.array([normalized(orig_plane[0]),[0,0,1],[0,0,1]])
    v2 = N.array([normalized(target_plane[0]),[0,0,1],[0,0,1]])
    q  = rotquat(v1,v2)
    R = N.transpose(N.matrix(rotmat(q[0])))
    
        
    orig_plane[1] = matrix2list(orig_plane[1] * R)
    
    mydot = vcos(orig_plane[1],target_plane[1])

    if 1 - abs(mydot) < 0.01:
        
        ## Entonces son o iguales o opuestos y el mod de cuaterniones falla
        if mydot >=0.99 and mydot <=1.01:
            ## Si son iguales el angulo es 0 y cos0 = 1
            ## No hay que modificar la matriz
            
            pass
        elif mydot <=-0.99 and mydot >=-1.01:
            ## Si son opuestos el angulo es -1
            ## Rotamos 180 alrededor del eje target 1 (voltearlo)
            
            R2 = N.transpose(N.matrix(rotmat(quaternion(target_plane[0],3.14159))))
            
            R = R * R2
    else:
        v1 = N.array([normalized(orig_plane[1]),[0,0,1],[0,0,1]])
        v2 = N.array([normalized(target_plane[1]),[0,0,1],[0,0,1]])
        q  = rotquat(v1,v2)

        R = R * N.transpose(N.matrix(rotmat(q[0])))
    
    
    
    for i in range(len(model.xyz)):
        model.xyz[i] = model.xyz[i] * R

    return model , R

    
def orientVectors(model,vorig,vdest):
    """
    Orient first vector so that has same direction as second vector.
    model - PDBModel
    vorig - (x,y,z), origin vector
    vdest - (x,y,z), second vector
    
    """
    mydot = vcos(vorig,vdest)
    
    if 1 - abs(mydot) < 0.01: 
        
        ## Entonces son o iguales o opuestos y el mod de cuaterniones falla
        if mydot >=0.99 and mydot <=1.01:
            ## Si son iguales el angulo es 0 y cos0 = 1
            ## No hay que modificar la matriz
            
           return N.matrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
            
        elif mydot <=-0.99 and mydot >=-1.01:
            ## Si son opuestos el angulo es -1
            ## Rotamos 180 alrededor de un eje arbitrario
            
            R = N.transpose(N.matrix(rotmat(quaternion(axis,3.14159))))
            
    else:  
        v1 = N.array([normalized(vorig),[0,0,1],[0,0,1]])
        v2 = N.array([normalized(vdest),[0,0,1],[0,0,1]])
        q  = rotquat(v1,v2)
        R = N.transpose(N.matrix(rotmat(q[0])))    
    if model != None:
        for i in range(len(model.xyz)):
            model.xyz[i] = model.xyz[i] * R
    
    return model,R
    

def scorePlane( planea,planeb):
    return  angle(cross(planea[0],planea[1]),cross(planeb[0],planeb[1]))
                   
def getN (residue):
    
    return residue.xyz[N.where(residue.maskFrom( 'name', ['N'] ))[0][0]]
    
def calcPlane (residue,Ca=None,C = None, do_norm = True):
    
    if C == None and Ca == None:
        sel_atoms = N.where(residue.maskFrom( 'name', ['N','CA','C'] ))
        Npos = sel_atoms[0][0]
        Capos= sel_atoms[0][1]
        Cpos = sel_atoms[0][2]
        Ni  = residue.xyz[Npos]
        Ca = residue.xyz[Capos]
        C  = residue.xyz[Cpos]
    else:
        Ni = residue
    
    if not do_norm:   
        return  [C-Ni,Ca-Ni]
    else:
        v1 = cutdepth(C-Ni,0.001)
        v2 = cutdepth(Ca-Ni,0.001)
        aux = cutdepth(cross(v1,v2),0.001)
        v2 = cutdepth(cross(aux,v1),0.001)
        return [cutdepth(normalized(v1),0.001),cutdepth(normalized(v2),0.001)]



def genPoints(model,chain = 0,atomname ='N'):
    model = model.compress( model.maskProtein())
    
    try:
        chain = model.takeChains([chain])
    except:
        return []
    
    Nfirst = chain.compress(chain.maskFrom( 'name', 'N' )).xyz[0]
    Cafirst = chain.compress(chain.maskFrom( 'name', 'CA' )).xyz[0]
    Cfirst = chain.compress(chain.maskFrom( 'name', 'C' )).xyz[0]
    
    chain = chain.compress(chain.maskFrom( 'name', atomname ))
    
    ## Center
    #~ print "**"
    aux = [chain.xyz[0][0],chain.xyz[0][1],chain.xyz[0][2]]
    for j in range(len(chain.xyz)):
        chain.xyz[j] = chain.xyz[j] - aux
    
    ## Rotate
    chain, R_chain = planarize(chain,calcPlane(Nfirst,Cafirst,Cfirst))
   
    points = chain.xyz
        
    #~ print points
    
    return  points


def expandedPoints(numAA,dist=2.65):
    ## It creates numAA+1 points to fit a structure
    
    point = [0.,0.,0.]
    tr = [dist,0.,0.]
    points=[N.array(point)]
    for i in range(1,numAA+1):
        points.append(N.array([0.,tr[0]*i,0.]))
        
    return points
    
def points2Pdb ( points,chain='A'):
    res = ""
    i=1
    for p in points:
        
        res += "ATOM    %s  CA  GLY %s %s     %s %s %s  1.00 60.55           C\n"%(pad(i,3),chain,pad(i,3),pad(p[0]),pad(p[1]),pad(p[2]))
                #%-d      33.852   5.514  53.987"
        i=i+1
    return res 


def doAAReorientation (residue = None):
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
    
    residue, R = planarize(residue,calcPlane(residue))
    
    

    return residue, R

def flip (model = None, orig_plane =[[1,0,0],[0,0,1]],vectors = []):
    target_plane = [[orig_plane[0][0],orig_plane[0][1],orig_plane[0][2]],[-orig_plane[1][0],-orig_plane[1][1],-orig_plane[1][2]]]
    
    model, R = planarize(model,orig_plane,False,target_plane,False)
    
    
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

def aalist2intlist(list):
    aas = MU.allAA()
    cor ={}
    for i in range(len(aas)):
        cor[aas[i]]=i
    del aas
    
    l = []
    for i in list:
        if i == 'WAT':
            l.append(-1)
        else:
            l.append(cor[MU.singleAA([i])[0]])
    del cor
    return l

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
            residues[i].xyz = residues[i].xyz - residues[i].xyz[0]
            nc = residues[i].xyz[2]-residues[i].xyz[0]
            nca = residues[i].xyz[1]-residues[i].xyz[0]
            planarize(residues[i],[nc,nca])[0].writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
        
        for i in range(len(residues)):
            doAAReorientation(residues[i])[0].writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1+11))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
        
        f.close()
        
        f = open(T.testRoot()+"/polysys/tots2.pdb","w")
        
        for i in range(len(residues)):
            res = residues[i].clone()
            flipAA(residues[i])[0].concat(flipAA(flipAA(res)[0])[0]).writePdb(T.testRoot()+"/polysys/"+str(i)+".pdb")
            f2 =  open(T.testRoot()+"/polysys/"+str(i)+".pdb","r")
            f.write("MODEL    %4d\n"%(i+1+22))
            f.writelines(f2.readlines()[:-1])
            f.write("ENDMDL\n")
            f2.close()
            
        f.close()
        
    def test_points2pdb(self):
        """points2pdb testcase"""
        
        points = [[0,0,0],[4,0,0],[8,0,0],[12,0,0],[16,0,0],[20,0,0],[20,4,0],[20,8,0],[20,12,0],[20,16,0],[20,20,0]]
        print points2Pdb(points)
        self.assertEqual(points2Pdb(points)[0:79],"ATOM      1  CA  GLY A   1           0       0       0  1.00 60.55           C\n")
    
    def test_CalcPlane(self):
        """CalcPlane test case"""
        res = PDBModel(T.testRoot()+"/polysys/res.pdb")
        
        plane = calcPlane(res)
        self.assertEqual(plane,[[0.99999999977321929, -0, 0.0], [0.0, 0.0, -1.0]] )
        
        
    def test_planarize2(self):
        """Second planarization test"""
        
        ## Try:
        ## planarize(planarize(X,a->b),b->a) == x
        
        res = PDBModel(T.testRoot()+"/polysys/res.pdb")
        res.xyz[0] = [0.,0.,0.]
        res.xyz[1] = [0.,-3.,4.]
        res.xyz[2] = [0.,1.,0.]
        plane = calcPlane(res)
        lastres= planarize(res,plane)[0]
        lastplane = calcPlane(lastres)
        finalres =  planarize(lastres,lastplane,False,plane,False)[0]
        self.assertAlmostEqual(res.xyz[0][0],finalres.xyz[0][0],3)
        self.assertAlmostEqual(res.xyz[0][1],finalres.xyz[0][1],3)
        self.assertAlmostEqual(res.xyz[0][2],finalres.xyz[0][2],3)
        self.assertAlmostEqual(res.xyz[1][0],finalres.xyz[1][0],3)
        self.assertAlmostEqual(res.xyz[1][1],finalres.xyz[1][1],3)
        self.assertAlmostEqual(res.xyz[1][2],finalres.xyz[1][2],3)
        self.assertAlmostEqual(res.xyz[2][0],finalres.xyz[2][0],3)
        self.assertAlmostEqual(res.xyz[2][1],finalres.xyz[2][1],3)
        self.assertAlmostEqual(res.xyz[2][2],finalres.xyz[2][2],3)
        
        
        
        res = PDBModel(T.testRoot()+"/polysys/res.pdb")
        res.xyz[0] = [0.,0.,0.]
        res.xyz[1] = [0.,0.,1.]
        res.xyz[2] = [0.,1.,0.]
        plane = calcPlane(res)
        lastres= planarize(res,plane)[0]
        lastplane = calcPlane(lastres)
        finalres =  planarize(lastres,lastplane,False,plane,False)[0]
        self.assertAlmostEqual(res.xyz[0][0],finalres.xyz[0][0],3)
        self.assertAlmostEqual(res.xyz[0][1],finalres.xyz[0][1],3)
        self.assertAlmostEqual(res.xyz[0][2],finalres.xyz[0][2],3)
        self.assertAlmostEqual(res.xyz[1][0],finalres.xyz[1][0],3)
        self.assertAlmostEqual(res.xyz[1][1],finalres.xyz[1][1],3)
        self.assertAlmostEqual(res.xyz[1][2],finalres.xyz[1][2],3)
        self.assertAlmostEqual(res.xyz[2][0],finalres.xyz[2][0],3)
        self.assertAlmostEqual(res.xyz[2][1],finalres.xyz[2][1],3)
        self.assertAlmostEqual(res.xyz[2][2],finalres.xyz[2][2],3)
        
        
        res = PDBModel(T.testRoot()+"/polysys/res.pdb")
        res.xyz[0] = [0.,0.,0.]
        res.xyz[1] = [1.,0.,0.]
        res.xyz[2] = [0.,0.,1.]
        plane = calcPlane(res)
        lastres= planarize(res,plane)[0]
        lastplane = calcPlane(lastres)
        finalres =  planarize(lastres,lastplane,False,plane,False)[0]
        self.assertAlmostEqual(res.xyz[0][0],finalres.xyz[0][0],3)
        self.assertAlmostEqual(res.xyz[0][1],finalres.xyz[0][1],3)
        self.assertAlmostEqual(res.xyz[0][2],finalres.xyz[0][2],3)
        self.assertAlmostEqual(res.xyz[1][0],finalres.xyz[1][0],3)
        self.assertAlmostEqual(res.xyz[1][1],finalres.xyz[1][1],3)
        self.assertAlmostEqual(res.xyz[1][2],finalres.xyz[1][2],3)
        self.assertAlmostEqual(res.xyz[2][0],finalres.xyz[2][0],3)
        self.assertAlmostEqual(res.xyz[2][1],finalres.xyz[2][1],3)
        self.assertAlmostEqual(res.xyz[2][2],finalres.xyz[2][2],3)
        
    def test_aalist(self):
        """aalist2intlist test"""
           
        print aalist2intlist(['WAT', 'LYS', 'LYS', 'LYS', 'WAT', 'WAT', 'ASP', 'WAT'])
if __name__ == '__main__':
    BT.localTest()    
        
