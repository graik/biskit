from Biskit import *
from numpy import *
from math import atan2
import cPickle 
import os
from Biskit.mathUtils import eulerRotation,cartesianToPolar
from math import acos,sqrt

from quaternion import *

def extractFromPDB ( path = "" ):
	modelraw = PDBModel(path)

	model = modelraw.compress( modelraw.maskProtein() )

	try:
		f = open ('./residues_db/residues',"r")
		total_res = cPickle.load(f)

	except (cPickle.UnpicklingError, IOError,EOFError):
		f = open ('./residues_db/residues',"w")
		total_res = {}

	f.close()

	for c in range( model.lenChains() ):
		chain = model.takeChains( [c] )
		
		for i in  range(1,chain.lenResidues()+1):
			res = where(chain['residue_number'] == i)
			res_mask = zeros( len(chain), bool )
			put( res_mask, res, 1) 
			residue = chain.compress(res_mask)
			r = residue.sequence()
			
			if not r in total_res:
				total_res[r] = 1
			else:
				total_res[r] += 1
			
			# cazar N y C para reorientar el aa
			
			iNC= where(residue.maskFrom( 'name', ['N','CA','C'] ))
			
			if iNC[0] != []:
				
				if( not os.access("./residues_db/"+r, os.F_OK)):
					os.mkdir("./residues_db/"+r)
				
				# orientar haciendo NC un vector sobre X
				
				N = residue.xyz[iNC[0][0]]
				residue.xyz = residue.xyz - N
				
				
				C = residue.xyz[iNC[0][2]]
				Ca = residue.xyz[iNC[0][1]]
				nNCNCa = normalized(cross(C,Ca))	
				NC = array([[0.,0.,0.],[0.,0.,0.],nNCNCa])
				NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,0.,-1.]])	
				q= rotquat(NC,NC2)[2]
				R = transpose(matrix(rotmat(q)))				
				for i in range(len(residue.xyz)):
					residue.xyz[i] = residue.xyz[i] * R
				residue.update()
				
				C = residue.xyz[iNC[0][2]]
				Cnorm = normalized(C)
				NC = array([[0.,0.,0.],[0.,0.,0.],Cnorm])
				NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,1.,0.]])
				q= rotquat(NC,NC2)[2]
				R = transpose(matrix(rotmat(q)))
				for i in range(len(residue.xyz)):
					residue.xyz[i] = residue.xyz[i] * R
				residue.update()
				
				
				residue.renumberResidues()
				residue['serial_number'] = range(1,residue.lenAtoms()+1)
				residue.writePdb("./residues_db/"+r+"/"+r+str(total_res[r])+".pdb")
				


	f = open ('./residues_db/residues',"w")
	cPickle.dump(total_res,f,cPickle.HIGHEST_PROTOCOL)
	f.close()
	
def rotateAA ( paa='',v1 = [0.,0.,0.],v2 = [0.,0.,0.],writeIt = False,name=''):
	
	if not isinstance(paa,PDBModel):
		aa = PDBModel(paa)
	else:
		aa = paa
		
	# we assume it's in the YZ plane ( it must be from DB )
	
	NC = array([[0.,0.,0.],[0.,0.,0.],v1])
	NC2 = array([[0.,0.,0.],[0.,0.,0.],v2])
	q= rotquat(NC,NC2)[2]
	#~ print "q:",q,rotquat(NC,NC2)
	R = transpose(matrix(rotmat(q)))	
	#~ print R,v1,v2
	for i in range(len(aa.xyz)):
		aa.xyz[i] = aa.xyz[i] * R
	aa.update()
	
	if writeIt :
		aa.writePdb(name)
	
	return aa

def stickAAs( a,b , flip ):
	iC= where(a.maskFrom( 'name', ['C'] ))
	#~ print a
	#~ print b
	
	if flip :
		
		b = rotateAA(b,[1.,0.,0.],[0.,0.,1.],True,'lol.pdb')
		b = rotateAA(b,[0.,0.,1.],[-1.,0.,0.],True,'lol.pdb')
		
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
	
	aCnNn = normalized(aCn - aNn)
	
	bC1 = b.xyz[iCb[0][0]] # C1
	bN1 = b.xyz[iNb[0][0]] # N1
	
	bC1N1 = normalized(bC1 - bN1)
	
	#center
	a.xyz = a.xyz - aNn
	b.xyz = b.xyz - bN1
	#~ print "centered-",bN1,bC1,iCb[0][0],iCa[0][0]
	#Reorient b 
	b = rotateAA(b,bC1N1,aCnNn)
	
	c =stickAAs(a,b,False)
	#~ print a 
	#~ print b
	#~ print c
	return c
	
	
	


def calculateBulkRMS(aaname = 'A', writeIt = False):
	
	try:
		f = open ('./residues_db/residues',"r")
		total_res = cPickle.load(f)
		f.close()
	except (cPickle.UnpicklingError, IOError,EOFError):
		print "[ERROR calculateBulkRMS] no residue database."
		return
	
	residue = None
	residues = []
	
	if writeIt:
		f = open ("./rmsd_"+aaname,"w")
	
	#load all aa of name aaname
	for i in range(1,total_res[aaname]+1):
		print "./residues_db/"+aaname+"/"+aaname+str(i)
		residues.append(PDBModel("./residues_db/"+aaname+"/"+aaname+str(i)+".pdb"))
		
	
	
	for i in range(total_res[aaname]):
		for j in range(total_res[aaname]):
			if i!= j:
				if writeIt:
					f.write(str(i)+" "+str(j)+ " " + str (residues[i].rms(residues[j]))+"\n")
		f.write("\n")
	
	if writeIt:
		f.close()
	
	return 

				
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
			
	print "angles",aX*180./pi,aZ*180./pi
	return (aX ,aZ)


#~ extractFromPDB ( '1HUY.pdb')
#~ extractFromPDB ( '1CV7.pdb')
#~ extractFromPDB ( '1GZX.pdb')
#~ flip = True
#~ ab = stickAAs( PDBModel('E1.pdb'),PDBModel('E6.pdb'),flip )
#~ ac = stickAAs( ab,PDBModel('E8.pdb'),False )
#~ ad = stickAAs( ac,PDBModel('E8.pdb'),flip )
#~ ad.writePdb("concat.pdb")
#~ ad.report()

#~ print ad._chainIndex
 
