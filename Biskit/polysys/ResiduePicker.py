from Biskit import *
from numpy import *
import cPickle 
import os
from Biskit.mathUtils import eulerRotation,cartesianToPolar
from math import acos,sqrt

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
			
			iNC= where(residue.maskFrom( 'name', ['N','C','CA'] ))
			
			if iNC[0] != []:
				
				if( not os.access("./residues_db/"+r, os.F_OK)):
					os.mkdir("./residues_db/"+r)
				
				# orientar haciendo NC un vector sobre X
				C = residue.xyz[iNC[0][1]]
				N = residue.xyz[iNC[0][0]]
				NC = N-C
				pNC = cartesianToPolar(array([NC]))
				R = rotation(pNC[0,1],0,pNC[0,2])
				aNC = matrix(NC)
				
				print "original"
				#~ print residue.xyz
				print pNC[0,1]*180./pi,pNC[0,2]*180./pi
				
				for i in range(len(residue.xyz)):
					residue.xyz[i] = residue.xyz[i] * R
				
				C = residue.xyz[iNC[0][2]]
				N = residue.xyz[iNC[0][0]]
				NC = N-C
				pNC = cartesianToPolar(array([NC]))
				R = rotation(pNC[0,1],0,pNC[0,2])
				aNC = matrix(NC)
				for i in range(len(residue.xyz)):
					residue.xyz[i] = residue.xyz[i] * R
				
				C = residue.xyz[iNC[0][1]]
				N = residue.xyz[iNC[0][0]]
				NC = N-C
				pNC = cartesianToPolar(array([NC]))
				print "orientado"
				#~ print residue.xyz
				print pNC[0,1]*180./pi,pNC[0,2]*180./pi
				
				residue.update()
				
				residue.xyz = residue.xyz - N
				residue.update()
				print "centrado"
				#~ print residue.xyz
				residue.renumberResidues()
				residue['serial_number'] = range(1,residue.lenAtoms()+1)
				residue.writePdb("./residues_db/"+r+"/oriented_"+r+str(total_res[r])+".pdb")
				
				#~ if( not os.access("./residues_db/"+r, os.F_OK)):
					#~ os.mkdir("./residues_db/"+r)
				#~ f = open ("./residues_db/"+r+"/"+r+str(total_res[r]),"w")
				#~ cPickle.dump(residue,f,cPickle.HIGHEST_PROTOCOL)
				#~ f.close()


	f = open ('./residues_db/residues',"w")
	cPickle.dump(total_res,f,cPickle.HIGHEST_PROTOCOL)
	f.close()

def rotation (alpha,beta,gamma):
	#alpha -> z gamma -> x beta ->y
	cos_alpha = cos(alpha); sin_alpha = sin(alpha)
	cos_beta  = cos(beta);  sin_beta  = sin(beta)
	cos_gamma = cos(gamma); sin_gamma = sin(gamma)
	R = matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	
	R[0,0] = cos_alpha*cos_beta
	R[0,1] = cos_alpha*sin_beta*sin_gamma - sin_alpha*cos_gamma
	R[0,2] = cos_alpha*sin_beta*cos_gamma + sin_gamma * sin_alpha

	R[1,0] = sin_alpha*cos_beta
	R[1,1] =  sin_alpha*sin_beta*sin_gamma  + cos_gamma * cos_alpha
	R[1,2] =  cos_gamma * sin_beta*sin_alpha - cos_alpha*sin_gamma

	R[2,0] = -sin_beta 
	R[2,1] = cos_beta * sin_gamma
	R[2,2] = cos_beta*cos_gamma
	
	return R


#~ extractFromPDB('1huy.pdb')
extractFromPDB('oriented_E1.pdb')
extractFromPDB('oriented_E2.pdb')

#~ NC = [1,1,1]				
#~ pNC = cartesianToPolar(array([NC]))
#~ print (NC[1]/sqrt(NC[1]*NC[1] + NC[0]*NC[0]))*180/pi, (NC[1]/sqrt(NC[1]*NC[1] + NC[2]*NC[2]))*180/pi
#~ R = rotation(-pNC[0,1],0,-pNC[0,2])
#~ aNC = matrix(NC)
#~ print aNC, pNC[0,1]*180./pi,pNC[0,2]*180./pi
#~ p = aNC*matrix(R)
#~ print p
#~ print rotation (1,1,1)