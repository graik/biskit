from math import cos, sin, tan, sqrt, atan2,acos, pi
import numpy as N



## EASY vector handling
def sub(a,b):
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]]

def add(a,b):
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]]
    
def norm ( v ):
    return sqrt(v[0]*v[0] + v[1]*v[1]+ v[2]*v[2])

def dot (v1,v2):
    return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

def minus(v1):
    return [-v1[0],-v1[1],-v1[2]]

def cross (a , b):
    return [(a[1]*b[2]) - (a[2]*b[1]), (a[2]*b[0]) - (a[0]*b[2]), (a[0]*b[1]) - (a[1]*b[0])]

def normalized(a):
    n = norm(a)
    return [a[0]/n,a[1]/n,a[2]/n]
    
def angle( v1, v2 ):
    
    a = (norm(v1)*norm(v2))
    
    b = dot (v1,v2)
    
    return acos(max(-1.0,min(b/a,1.0)))

def vcos( v1, v2 ):
    
    a = (norm(v1)*norm(v2))
    
    b = dot (v1,v2)
    
    return max(-1.0,min(b/a,1.0))


def rotationMatrix (alpha=0.,beta=0.,gamma=0.):
    #alpha -> z gamma -> x beta ->y
    cos_alpha = cos(alpha); sin_alpha = sin(alpha)
    cos_beta  = cos(beta);  sin_beta  = sin(beta)
    cos_gamma = cos(gamma); sin_gamma = sin(gamma)
    R = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    
    R[0][0] = cos_alpha*cos_beta
    R[0][1] = cos_alpha*sin_beta*sin_gamma - sin_alpha*cos_gamma
    R[0][2] = cos_alpha*sin_beta*cos_gamma + sin_gamma * sin_alpha

    R[1][0] =  sin_alpha * cos_beta
    R[1][1] =  sin_alpha * sin_beta * sin_gamma  + cos_gamma * cos_alpha
    R[1][2] =  cos_gamma * sin_beta * sin_alpha - cos_alpha * sin_gamma

    R[2][0] = -sin_beta 
    R[2][1] = cos_beta * sin_gamma
    R[2][2] = cos_beta * cos_gamma
    
    return R
    
def matrix2list(m):
    return N.array([m[0,0],m[0,1],m[0,2]])    
    
def matrixMult(a,b):
    res = [0.,0.,0.]
    
    res[0] = a[0]*b[0][0]+a[1]*b[1][0]+a[2]*b[2][0]
    res[1] = a[0]*b[0][1]+a[1]*b[1][1]+a[2]*b[2][1]
    res[2] = a[0]*b[0][2]+a[1]*b[1][2]+a[2]*b[2][2]
    
    return res

def sphericalAngles( coord= [1.,0.,0.]):
    
    aX= atan2(abs(coord[2] ),abs(coord[1]))
    if coord[2] < 0. :
        aX =  - aX
    
    norm = sqrt(coord[1]*coord[1] + coord[0]*coord[0])
    
    if norm < 0.00001 :
        aZ = 0.
        
    else:
        norm2 =  coord[1] /norm
        if  norm2>1:
            norm2 = 1
        if  norm2<-1:
            norm2 = -1
        
        aZ = acos(norm2)
        if coord[0] <0:
            aZ = 2*pi - aZ
            
    #~ print "angles",aX*180./pi,aZ*180./pi
    return (aX ,aZ)
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T


class Test(BT.BiskitTest):
    """ Test cases for fold data creation"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    
    def test_basics(self):
        """Basic ops test case"""
        a = [1.,0.,0.]
        b = [1.,0.,1.]
        c = [0.,1.,0.]
        
        self.assertEqual(sub(a,b),[0.,0.,-1])
        self.assertEqual(add(a,b),[2.,0.,1])
        self.assertEqual(cross(a,c),[0.0, 0.0, 1.0])
        self.assertEqual(dot(a,c),0)
        self.assertEqual(normalized(b),[0.70710678118654746, 0.0, 0.70710678118654746])
        print sphericalAngles([-1,-1,-1])
        print angle([3.3629, -0.92, -0.26], [2.4797924, 0.86277348, 1.0560791] )
    
    def test_3D(self):
        """3D ops test case"""
        a = [1.,0.,0.]
        b = [0.,0.,1.]
        c = [0.,1.,0.]
        
        R = rotationMatrix(pi/2,0.,0.)
        d = matrixMult(a,R)
        print d
        
        R = rotationMatrix(0.,pi/2,0.)
        d = matrixMult(a,R)
        print d
        
        R = rotationMatrix(0.,0.,pi/2)
        d = matrixMult(a,R)
        print d
    

    
if __name__ == '__main__':
    BT.localTest()    
 
