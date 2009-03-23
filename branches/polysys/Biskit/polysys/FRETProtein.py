 
from Biskit import PDBModel
from numpy import*

class Interval:
	
	def __init__(self,s,e):
		self.start = s
		self.end  = e
		self._Interval__checkConsistence()
	
	def __checkConsistence(self):
		if( self.start<0 or self.end<0 ) :
			print "[ERROR Interval] Start or end can't be negative numbers when defining the range."
			self.start = 0
			self.end  = 0
		if(self.start>self.end):
			i = self.start
			self.start=self.end
			self.end=i
			print "[WARNING Interval",self.__str__(),"] Swapping start and end."
			
	def checkOverlap(self, others = []):
		for i in others:
			if (not ((self.end >= i.end and self.start >= i.end) or (self.start<=i.start and self.end <= i.start))):
				return True
		return False
		
	def __len__(self):
		return self.end - self.end
	
	def __str__(self):
		return '[Start: %d , End: %d]'%(self.start, self.end)
		
	def __cmp__(self,o):
		if o.start == self.start and o.end == self.end :
			return 0
		else :
			return -1

class BlockEntity:
	def __init__(self,name= "Unnamed",father=None):
		self.name = name
		self.sons = []
		self.father = father

	def __str__(self):
		return "[Entity: "+self.name+"]"
		
	def onInsertion(self,assembly= None):
		pass
	
	def run(self):
		pass
	
class Block(object):
	
	#
	#If no interval is given, the whole length of entity profiles or 1 is defaulted
	#If no entity is given, an interval of length 1 will be defaulted
	#
	
	def __init__(self,name,entity = None , father = None):
		self.name = name
		assert not isinstance(self,Interval) ,'[ERROR Block __init__] Entity is not a BlockEntity'
		
		self.entity = entity
		
		if(entity != None ):
			self.entity.father = father
			if(father!=None):
				self.entity.father.sons.append(entity)
	
		self.intervals = []
		
	
	def __dettachFromFather(self):
		if(father!=None):
			self.father.sons.remove ( self.entity )
			self.entity.father = None
	
	def setEntity (self,e):
		#reparent
		myfather = self.entity.father
		if(self.entity.father!=None):
			self.__dettachFromFather()
			myfather.sons.append(e)	
		self.entity = e;
		
		
		
	def addInterval(self, s,e=-1):
		if e==-1:
			assert not isinstance(self,Interval) ,'[ERROR Block addInterval] Please check parameter types.'
			i = s
			
		else:
			i = Interval(s,e)
			
		if i.checkOverlap(self.intervals) :
			print "[WARNING Block addInterval] Interval overlap for the same Block"
	
		self.intervals.append(i)	
		
	def remInterval(self, s,e = -1):
		if e==-1:
			assert not isinstance(self,Interval) ,'[ERROR Block addInterval] Please check parameter types.'
			i1 = s
		else:
			i1 = Interval(s,e)
			
		for i in self.intervals:
			if i == i1:
				self.intervals.remove(i1)
				return
		
	def __cmp__(self,o):
		if self.name == o.name :
			return 0
		else:
			return -1
	
	def __str__(self):
		mystring = "[Block "+self.name+"]"
		if self.intervals != []:
				for i in self.intervals:
					mystring += i.__str__()
		if self.entity != None :
			mystring += self.entity.__str__()+"\n"
			
			if self.entity.sons != [] :
				for i in self.entity.sons :
					mystring +="--->SubBlock : "+ i.__str__() +'\n'
		
		return mystring
	
class FRETProtein (PDBModel,BlockEntity):
	def __init__ (self,name,source, qy = 0, lt = 0,chromophore = None):
		PDBModel.__init__(self,source)
		BlockEntity.__init__(self,name)
		self.quantumYield=qy
		self.lifetime = lt
		self.name = name
		self.chromophore = chromophore
	
	def __str__(self):
		return BlockEntity.__str__(self)




class ConstraintAngle:
	def __init__(self,x,y,z):
		self.x = x;self.y=y;self.z =z

class ConstraintGlue:
	pass

class ConstraintClash(object):
	def __init__(self,typea,typeb):
		self.typea=typea
		self.typeb=typeb
	
	def __getType(self):
		if (typea,typeb) in c_clashes:
			return (typea,typeb)
		else:
			return (typeb,typea)
	
	type = property(fget=__getType)
	
	def __cmp__(self,o):
		if (self.typea==o.typea and self.typeb == o.typeb) or (self.typea==o.typeb and self.typeb == o.typea):
			return 0
		else:
			return -1
	
	def __str__(self):
		return "[Clash: "+self.typea+","+self.typeb+"]"

c_types = ("distance","orientation","glued","linked","max_distance","min_distance")

#distance: one is at certain distance of the other
#glued: when moving one,the other is always moved,so distance between them is never changed
#linked: ending and starting aminoacids are connected

c_clashes = ( ConstraintClash("distance","distance"),ConstraintClash("distance","glued"),ConstraintClash("orientation","orientation"))
	

class Constraint (object):
	
	def __init__(self,type,a = None,b = None,value =None):
		self.a = a
		self.b = b
		self.type = type
		self.value = value
	
	def sameOps (self, o):
		return ((self.a == o.a and self.b == o.b)or (self.a == o.b and self.b == o.a))
	
	def __cmp__(self,o):
		if self.value == o.value and self.type == o.type and ((self.a == o.a and self.b == o.b)or (self.a == o.b and self.b == o.a)) :
			return 0
		else:
			return -1

	def __str__(self):
		s=""
		s+="[Constraint: "+self.type+"] "
		if self.a != None:
			s+=self.a.name+" "
		if self.b!=None:
			s+=self.b.name
		return s

		
		
		
		
class Assembly :
	def __init__(self,name):
		self.blocks = []
		self.constraints = []
		self.name = name
		self.engines = []
		
	def addBlock(self, f):
		self.blocks.append(f)
		if self.blocks[-1].entity != None:
			self.blocks[-1].entity.onInsertion(self)
	
	def remBlock(self,f):
		if isinstance(f,Block):
			self.blocks.remove(f)
		else:
			self.blocks.remove(Block(f))
		
	def addConstraint(self,c):
		self.constraints.append(c)
	
	def remConstraint(self,c):
		self.constraints.remove(c)
	
	def addEngine(self,e):
		self.engines.append(e)
	
	def remEngine(self,e):
		self.engines.remove(e)
		
	def addAssembly(self,a,addEngines= False):
		self.blocks.extend (a.blocks)
		self.constraints.extend (a.constraints)
		if (addEngines ):
			self.engines.extend (a.engines)
	
	def run(self):
		#Check blocks
		reblocks =  []
		for i in range(len(self.blocks)-1):
			reblocks.append(self.blocks[i])
			
			if reblocks.count(self.blocks[i])>1 :
				print "[WARNING Polymer _build] Duplicated Block("+self.blocks[i].name+ ").Merging its intervals"
				pos = reblocks.index(self.blocks[i])
				for j in self.blocks[i].intervals:
					reblocks[pos].addBlock(blocks[i].intervals[j])
				
		#Check constraint coherence
		for i in range (len(self.constraints)-1):
			for j in range(i+1,len(self.constraints)-1):
				clash = ConstraintClash(self.constraints[i].type,self.constraints[j].type)
				if( clash in c_clashes ) :
					#Solve possible constraint clash
					if self.constraints[i].type == "distance":
						if  clash == ConstraintClash("distance","distance"):
							pass
							
						elif self.constraints[i].type == "distance"  and clash == ConstraintClash("distance","glued"):
							pass
							
						else:
							print "[ERROR Constraint _build] Unhandled clash descriptor (" +clash.__str__() +")." 
					
					elif self.constraints[i].type == "orientation":
						pass
						
					elif self.constraints[i].type == "glued":
						pass
						
					else:
						print "[ERROR Constraint _build] Unhandled constraint descriptor (" +self.constraints[i].type +")." 
		
		if( self.engines != []):
			for e in self.engines:
				e.run(self)
		
		
	def __str__(self):
		mystring = ""
		mystring+= "------------------------------------\n"
		mystring+= "Assembly:"+self.name+"\n"
		mystring+= "------------------------------------\n"
		
		if self.blocks != []:
			for i in self.blocks:
				mystring+= i.__str__()+"\n"
					
		if self.constraints != []:
			for i in self.constraints:
				mystring+= i.__str__()+"\n"
		mystring+= "------------------------------------\n"
		
		return mystring
			


class Engine(object):
	def __init__(self,a):
		self.myAssembly = a
	
	def run (self):
		pass
	
class CreateLinkEngine( Engine ):
	def __init__ (self,a):
		Engine.__init__(self,a)
		
	def run (self):
		Engine.run(self)
		#And then, let's do something



