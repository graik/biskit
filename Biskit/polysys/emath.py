
from math import cos, sin, tan, sqrt, atan2,acos, pi
from numpy import matrix
import numpy as N

def norm ( v ):
	return N.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def dot (v1,v2):
    return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

def vectorangle( v1, v2 ):
	
	a = (norm(v1)*norm(v2))
	
	b = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
	
	return acos(max(-1.0,min(b/a,1.0)))

def vectorcos( v1, v2 ):
	
	a = (norm(v1)*norm(v2))
	
	b = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
	
	return max(-1.0,min(b/a,1.0))
def rotation (alpha=0.,beta=0.,gamma=0.):
	#alpha -> z gamma -> x beta ->y
	cos_alpha = cos(alpha); sin_alpha = sin(alpha)
	cos_beta  = cos(beta);  sin_beta  = sin(beta)
	cos_gamma = cos(gamma); sin_gamma = sin(gamma)
	R = matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	
	R[0,0] = cos_alpha*cos_beta
	R[0,1] = cos_alpha*sin_beta*sin_gamma - sin_alpha*cos_gamma
	R[0,2] = cos_alpha*sin_beta*cos_gamma + sin_gamma * sin_alpha

	R[1,0] = sin_alpha * cos_beta
	R[1,1] =  sin_alpha * sin_beta * sin_gamma  + cos_gamma * cos_alpha
	R[1,2] =  cos_gamma * sin_beta * sin_alpha - cos_alpha * sin_gamma

	R[2,0] = -sin_beta 
	R[2,1] = cos_beta * sin_gamma
	R[2,2] = cos_beta * cos_gamma
	
	return R
	
def normalized(a):
	norm = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
	return [a[0]/norm,a[1]/norm,a[2]/norm]

def vectpermatrix(a,b):
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



""" 
cbrt(x) = x^{1/3}, if x >= 0 = -|x|^{1/3}, if x < 0 
""" 

def cbrt(x): 
	from math import pow 
	if x >= 0: 
		return pow(x, 1.0/3.0) 
	else:
		return -pow(abs(x), 1.0/3.0) 
	

	
""" 
Convert from polar (r,w) to rectangular (x,y) x = r cos(w) y = r sin(w)
""" 
def rect(r, w, deg=0): # radian if deg=0; degree if deg=1 
	from math import cos, sin, pi 
	if deg: 
		w = pi * w / 180.0 
		return r * cos(w), r * sin(w) 
		
"""
Convert from rectangular (x,y) to polar (r,w) r = sqrt(x^2 + y^2) w = arctan(y/x) = [-\pi,\pi] = [-180,180] 
""" 
def polar(x, y, deg=0): # radian if deg=0; degree if deg=1 
	from math import hypot, atan2, pi 
	if deg: 
		return hypot(x, y), 180.0 * atan2(y, x) / pi 
	else: 
		return hypot(x, y), atan2(y, x) 
		

""" 
x^2 + ax + b = 0 (or ax^2 + bx + c = 0) By substituting x = y-t and t = a/2,
the equation reduces to y^2 + (b-t^2) = 0 which has easy solution y = +/- sqrt(t^2-b) 
""" 
def quadratic(a, b, c=None): 
	import math, cmath 
	if c: # (ax^2 + bx + c = 0) 
		a, b = b / float(a), c / float(a)
	t = a / 2.0 
	r = t**2 - b 
	if r >= 0: # real roots 
		y1 = math.sqrt(r) 
	else: # complex roots 
		y1 = cmath.sqrt(r) 
	y2 = -y1 
	return y1 - t, y2 - t 
	
"""
x^3 + ax^2 + bx + c = 0  (or ax^3 + bx^2 + cx + d = 0)
With substitution x = y-t and t = a/3, the cubic equation reduces to    
    y^3 + py + q = 0,
where p = b-3t^2 and q = c-bt+2t^3.  Then, one real root y1 = u+v can
be determined by solving 
    w^2 + qw - (p/3)^3 = 0
where w = u^3, v^3.  From Vieta's theorem,
    y1 + y2 + y3 = 0
    y1 y2 + y1 y3 + y2 y3 = p
    y1 y2 y3 = -q,
the other two (real or complex) roots can be obtained by solving
    y^2 + (y1)y + (p+y1^2) = 0
"""
def cubic(a, b, c, d=None):
	from math import cos
	if d:			# (ax^3 + bx^2 + cx + d = 0)
		a, b, c = b / float(a), c / float(a), d / float(a)
	t = a / 3.0
	p, q = b - 3 * t**2, c - b * t + 2 * t**3
	u, v = quadratic(q, -(p/3.0)**3)
	if type(u) == type(0j):	# complex cubic root
		r, w = polar(u.real, u.imag)
		y1 = 2 * cbrt(r) * cos(w / 3.0)
	else:			# real root
		y1 = cbrt(u) + cbrt(v)
	y2, y3 = quadratic(y1, p + y1**2)
	return y1 - t, y2 - t, y3 - t
