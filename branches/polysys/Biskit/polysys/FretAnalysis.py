from FRET import FRET
from math import cos,sqrt,pi , sin,fabs, tan, acos
from numpy import max, min,array,matrix,transpose
from emath import vectpermatrix,rotation,vectorangle,norm
 
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
			v = array([1,0,0])*matrix(rotation(beta,alpha,0))
			
			#~ print v[0,0], v[0,1],v[0,2]
			#~ print sqrt(v[0,0]**2 + v[0,1]**2+ v[0,2]**2)
			
			vectors.append((v[0,0], v[0,1],v[0,2]))
			#~ vectors.append((-v[0,0], v[0,1],v[0,2]))
			beta+=incb
		alpha +=inca
		vectors.append((2, 2,2))
	d += incd	

#~ mCitrine:mCerulean
tauDA = 2.89
tauA = 3.6
tauD = 3.5
overlap = 31.178
overlap = 14618e-13
qyD = 0.62
 
f = FRET( tauA, tauD, tauDA, overlap, qyD)

distance = [0.,4.,0.]
#donor TM vector is (0,0,1)
others = []
for v in vectors:
	#donor-acceptor angle
	donor = [0.,0.,1.]
	acceptor = v
	thetaT = vectorangle(donor, acceptor)

	#donor-distance angle, distance vector is along acceptor y axis (0,y,0)
	thetaD = vectorangle(donor, distance)

	#acceptor-distance angle
	thetaA = vectorangle(acceptor,distance)


	r = norm(distance)
	
	#~ print v[0],v[1],v[2],f.energyTransferRate(thetaA,thetaD,thetaT,r)
	
	if v == (2,2,2):
		others.reverse()
		for v in others:
			print v[0],v[1],v[2],v[3]
		others =[]
		#~ print
	else:
		f.energyTransferRate(thetaA,thetaD,thetaT,r)
		#~ print f.kappa2
		#~ print v[0],v[1],v[2],f.energyTransferRate(thetaA,thetaD,thetaT,r)
		#~ others.append(( v[0],v[1],v[2],f.energyTransferRate(thetaA,thetaD,thetaT,r)))

for v in vectors:
	#donor-acceptor angle
	donor = [0.,0.,1.]
	acceptor = (-v[0], v[1],v[2])
	thetaT = vectorangle(donor, acceptor)

	#donor-distance angle, distance vector is along acceptor y axis (0,y,0)
	thetaD = vectorangle(donor, distance)

	#acceptor-distance angle
	thetaA = vectorangle(acceptor,distance)


	r = norm(distance)
	
	#~ print v[0],v[1],v[2],f.energyTransferRate(thetaA,thetaD,thetaT,r)
	
	#~ if v == (2,2,2):
		#~ print
	#~ else:
		#~ print -v[0],v[1],v[2],f.energyTransferRate(thetaA,thetaD,thetaT,r)
	


print f.energyTransferEfficiency(thetaA,thetaD,thetaT,1.33)
exit()