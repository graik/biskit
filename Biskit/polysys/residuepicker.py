from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,matrix
import cPickle 
import os
from math import acos,sqrt,atan2
import emath 
from quaternion import rotquat, rotmat

class residuePicker:
	
	def __init__(self):
		try:
			f = open ('./residues_db/residues',"r")
			self.total_res = cPickle.load(f)

		except (cPickle.UnpicklingError, IOError,EOFError):
			if( not os.access("./residues_db/", os.F_OK)):
				os.mkdir("./residues_db/")
			f = open ('./residues_db/residues',"w")
			self.total_res = {}

		f.close()
		
		try:
			f = open ('./residues_db/apvect',"r")
			self.apvect = cPickle.load(f)

		except (cPickle.UnpicklingError, IOError,EOFError):
			f = open ('./residues_db/apvect',"w")
			self.apvect = {}

		f.close()

	def save(self):
		f = open ('./residues_db/residues',"w")
		cPickle.dump(self.total_res,f,cPickle.HIGHEST_PROTOCOL)
		f.close()
		
		f = open ('./residues_db/apvect',"w")
		cPickle.dump(self.apvect,f,cPickle.HIGHEST_PROTOCOL)
		f.close()

	def extractFromPDB ( self,path = "" ):
		modelraw = PDBModel(path)

		model = modelraw.compress( modelraw.maskProtein() )

		for c in range( model.lenChains() ):
			chain = model.takeChains( [c] )
			
			residuos = chain.resModels()
			
			
			for i in range(len(residuos)) :
				
				r = residuos[i].sequence()
				
				residue = residuos[i]
				
				if i<len(residuos)-1 :
					nextres = residuos[i+1]
				else:
					nextres = None

				if not r in self.total_res:
					self.total_res[r] = 1
				else:
					self.total_res[r] += 1
				
				# cazar N y C para reorientar el aa
				
				iNC= where(residue.maskFrom( 'name', ['N','CA','C'] ))
				
				if nextres != None:
					iNC2= where(residue.maskFrom( 'name', ['N','CA','C'] ))
				
				if iNC[0] != []:
					
					if( not os.access("./residues_db/"+r, os.F_OK)):
						os.mkdir("./residues_db/"+r)
					
					# orientar haciendo NC un vector sobre X
					N = residue.xyz[iNC[0][0]]
					residue.xyz = residue.xyz - N
					
					
					C = residue.xyz[iNC[0][2]]
					Ca = residue.xyz[iNC[0][1]]
					nNCNCa = emath.normalized(cross(C,Ca))	
					NC = array([[0.,0.,0.],[0.,0.,0.],nNCNCa])
					NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,0.,-1.]])	
					q= rotquat(NC,NC2)[2]
					R = transpose(matrix(rotmat(q)))				
					for i in range(len(residue.xyz)):
						residue.xyz[i] = residue.xyz[i] * R
					residue.update()
					
					C = residue.xyz[iNC[0][2]]
					Cnorm = emath.normalized(C)
					NC = array([[0.,0.,0.],[0.,0.,0.],Cnorm])
					NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,1.,0.]])
					q= rotquat(NC,NC2)[2]
					R = transpose(matrix(rotmat(q)))
					for i in range(len(residue.xyz)):
						residue.xyz[i] = residue.xyz[i] * R
					residue.update()
				
					if nextres  != None :
						N2 = residue.xyz[iNC2[0][0]]
						if not nextres.sequence() in self.apvect:
							self.apvect[nextres.sequence()] = [ ]
						
						self.apvect[nextres.sequence()] .append( ( r , N2 - C))
						
					
					residue.renumberResidues()
					residue['serial_number'] = range(1,residue.lenAtoms()+1)
					residue.writePdb("./residues_db/"+r+"/"+r+str(self.total_res[r])+".pdb")

		self.save()
		
	def rotateAA (self, paa='',v1 = [0.,0.,0.],v2 = [0.,0.,0.],writeIt = False,name=''):
		
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

	def stickAAs( self,a,b , flip ):
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
		
	def concatenateAAChains (self,a,b) :
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


	def calculateBulkRMS(self,aaname = 'A', writeIt = False):
		
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

					
	

r = residuePicker()

r.extractFromPDB ( '1HUY.pdb')
#~ extractFromPDB ( '1CV7.pdb')
#~ extractFromPDB ( '1GZX.pdb')
#~ flip = True
#~ ab = stickAAs( PDBModel('E1.pdb'),PDBModel('E6.pdb'),flip )
#~ ac = stickAAs( ab,PDBModel('E8.pdb'),False )
#~ ad = stickAAs( ac,PDBModel('E8.pdb'),flip )
#~ ad.writePdb("concat.pdb")
#~ ad.report()

#~ print ad._chainIndex
 
