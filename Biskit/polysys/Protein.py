from Biskit import *
from FRETProtein import *
import re
from Biskit import molUtils as mu
from numpy import *
from string import *
from Bio import pairwise2

class UndefProtein (BlockEntity):
	def __init__(self,sequence = ""):
		self.sequence = sequence
		
	def setSequence(self,sequence):
		self.sequence = sequence
		
	def run (self):
		BlockEntity.run(self)
		#create atomic data from sequence and return a PDBModel
		
		

class Protein (PDBModel, BlockEntity):
	
	def __init__(self,name ="Protein",path=None,method = 'correlation'):
		self.undef_sequence = ""
		BlockEntity.__init__(self,name)
		PDBModel.__init__(self,path)
		self.method = method
		
		self.__prot_chains = []
		self.__prot_chains_current = -1
		self.__prot_chains_sequences = []
		
		#Get SEQRES sequence
		
		if(path!=None):
			f = open (path,"r")
			self.lineas = f.readlines()
			f.close()
			
			for i in self.lineas:
				o = re.match(".*SEQRES\s*([\d]*)\s*([A-Z]+)\s*([\d]+)\s*([A-Z\s]+).+",i)
					
				if o!= None:
					
					if self.__prot_chains_current == -1 :
						self.__prot_chains_current+=1
						self.__prot_chains.append(o.group(2))
						self.__prot_chains_sequences.append("")
					
					if self.__prot_chains_current >= 0 and self.__prot_chains[self.__prot_chains_current]!=o.group(2) :
						self.__prot_chains_current+=1
						self.__prot_chains.append(o.group(2))
						self.__prot_chains_sequences.append("")
					
					self.__prot_chains_sequences[self.__prot_chains_current]+=self.getSeq(o.group(4))
		
	def onInsertion(self,myassembly=None):
		a = Assembly("dummy")
		
		BlockEntity.onInsertion(self)
		
		chains=[]
		
		#Extract the undefined sequences and create the undefined blocks
		chain_number=0
		
		for i in self.__prot_chains:
			chain_string = self.takeChains( [chain_number], breaks=1 ).sequence()
			
			if (self.method == 'correlation'):
				alignment =  self.sequenceCorrelation(chain_string,self.__prot_chains_sequences[chain_number])
				print alignment
				start = int(alignment[1][0])
				end = start + len(self.takeChains( [chain_number] ).sequence())
				print start,end, len(chain_string)
				
			elif (self.method =='alignment'):
				alignment = self.sequenceAlignment(chain_string,self.__prot_chains_sequences[chain_number])
				start = alignment[0]
				end = alignment[1]
				print start,end, len(chain_string)
			else:
				print '[ERROR Protein __Init__] Unsupported method name:',method
			
			if end > len(self.__prot_chains_sequences[chain_number])-1:
				print "[WARNING Protein:",self.name," onInsertion] SEQRES has undefined residues than the ones defined in ATOM"
			else:
				bstart = None; bend = None; bchain = None
				e1=None;e2=None;e3=None
				if start >0:
					#create the proximal undefined seq
					e1 = BlockEntity(self.name+"_undef_chain_start_"+self.__prot_chains[chain_number])
					bstart = Block(self.name+"_undef_chain_start_"+self.__prot_chains[chain_number],e1,self)
					bstart.addInterval(Interval(0,start-1))
					a.blocks.append(bstart)
				if end<len:
					#create distal undef seq
					e2 = BlockEntity(self.name+"_undef_chain_end_"+self.__prot_chains[chain_number])
					bend = Block(self.name+"_undef_chain_end_"+self.__prot_chains[chain_number],e2,self)
					bend.addInterval(Interval(end+1,len(self.__prot_chains_sequences[chain_number])-1))
					a.blocks.append(bend)
				#create chain seq
				e3 = BlockEntity(self.name+"_def_chain_"+self.__prot_chains[chain_number])
				bchain = Block(self.name+"_def_chain_"+self.__prot_chains[chain_number],e3,self)
				bchain.addInterval(Interval(start,start+len(chain_string)))
				a.blocks.append(bchain)
				chains.append(bchain)
				
				#add constraints
				if e1!=None:
					a.addConstraint(Constraint("linked",e1,e3))
					
				if e3!=None:
					a.addConstraint(Constraint("linked",e3,e2))
				
				if(chain_number>0):
					a.addConstraint(Constraint("glued",chains[0],chains[chain_number]))
				
				chain_number +=1
		
		if(myassembly!=None):
			myassembly.addAssembly(a)
	
	def getSeq(self,s=""):
		aas = s.split()
		s = ""
		for i in mu.singleAA(aas):
			s+= i
		return s
	
	def sequenceCorrelation	(self, a="",b = "",normalize = False):
		al = []
		bl = []
		
		if (b ==""):
			b= self.takeChains( [0] ).sequence()
		
		print a, b
		
		for i in a:
			al+= [float(ord(i)) or 0]
		
		for i in b:
			if( i in a) :
				bl+= [float(ord(i)) or 0]
			else:
				bl+=[0]
		
		a1 = array(al)
		b1 = array(bl)
		ab = correlate(a1,b1)
		
		if (normalize):
			c = ab /(correlate(a1,a1)*correlate(b1,b1))
			return (max(c) , where(c == max(c)))
		else:
			return (max(ab) , where(ab==max(ab)))	

	def sequenceAlignment (self, a= "",b= ""):
		c = b
		if (b ==""):
			c= self.takeChains( [0] ).sequence()
		
		(a,b,c,d,e) = pairwise2.align.localms(c, a, 3, -1,-2,-2)[0]
		
		return (d,e)
		
	def __str__(self):
		return BlockEntity.__str__(self)+PDBModel.__str__(self)
	

#~ p = Protein("testProtein",'2Q57.pdb')
#~ p.onInsertion()

#~ GGG 
#~ MSKGEELFTGVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTL?VQCFARYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNAISDNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGI

#~ a = "GLY CYS CYS GLY GLY CYS CYS GLY"
#~ b = "GLY VAL VAL"

#~ print p.sequenceCorrelation("GVVPILKELDG","")
#~ print p.sequenceAlignment("GVVPILKELDG","")
