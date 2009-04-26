from FRET import FRET
from math import cos,sqrt,pi , sin,fabs, tan, acos
from numpy import max, min,array,matrix,transpose
from emath import vectpermatrix,rotation,vectorangle,norm
 
 
 #~ for gnuplot  set pm3d depthorder
#~ set view map
#~ set ticslevel 0
#~ set xrange [0:1]
# Calculate sphere points for vector generation.

inca = 5.* pi /180.
incb = 5.* pi /180.
dmax = 1.
incd = 1.

vectors = []

d= 0.
while d < dmax:
	alpha = -pi/2.
	while alpha <pi/2.:
		beta = -pi/2.
		while beta < pi/2.:
			v = array([1.,0,0])*matrix(rotation(beta,alpha,0.))
			vectors.append((v[0,0], v[0,1],v[0,2]))
			beta+=incb
		alpha +=inca
		vectors.append((2, 2,2))
	d += incd	

#~ mCitrine:mCerulean
tauDA = 2.89
tauA = 3.6
tauD = 3.5
overlap = 959732990869504
qyD = 0.62
 
f = FRET( tauA, tauD, tauDA, overlap, qyD,44e3)

distance = [0.,60.,0.]
r = norm(distance)

#donor TM vector is (0,0,1)

others = []

for v in vectors:
	#donor-acceptor angle
	donor = [0.,0.,1.]
	acceptor = v

	if v == (2,2,2):
		others.reverse()
		for v in others:
			acceptor = (-v[0], v[1],v[2])
			
			f.calcK2( donor ,acceptor, distance)
			
			print -v[0],v[1],v[2],f.energyTransferEfficiency(r)[0]
			
		others =[]
		
		print
		
	else:
		f.calcK2( donor ,acceptor, distance)
		
		print v[0],v[1],v[2],f.energyTransferEfficiency(r)[0]
		
		others.append(( v[0],v[1],v[2],f.energyTransferEfficiency(r)[0]))

#~ f.k2  = 2./3.
#~ print f.calcR0()
#~ print f.energyTransferEfficiency(r)
#~ print 47.16**6 / (40**6+47**6)
