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
	def __init__(self,name= "Empty",father=None):
		self.name = name
		self.sons = []
		self.father = father
		

	def __str__(self):
		return "[Entity: "+self.name+"]"
		
	def onInsertion(self,assembly= None):
		pass
	
	def run(self):
		# As it's empty it tries to return its father
		# An absolute empty entity has no sense otherwise.
		return self.father
	
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
		if not type in c_types:
			print "[WARNING Constraint _init_] "+type + " is not an allowed type name."
		if self.a == self.b:
			print "[ERROR Constraint _init_] Entities of a constraint can not be the same."
	
	def sameOps (self, o):
		return ((self.a == o.a and self.b == o.b)or (self.a == o.b and self.b == o.a))
	
	def triangularLink(self, c):
		if c.type != "linked" or self.type != "linked" :
			print "[WARNING Constraint triangularLink] "+str(c) + " is not a link constraint."
			return False
		
		if self.a == c.a or self.b == c.b:
			return True
		else:
			return False
		
	def linked (self,c):
		if c.type != 'linked' or self.type != "linked":
			print "[WARNING Constraint linked] (type)"+ c.type+"  "+str(c) + " is not a link constraint."
			return False
	
		if self.b == c.a:
			return True
		else:
			return False
		
	def ilinked (self,c):
		if c.type != 'linked' or self.type != "linked":
			print "[WARNING Constraint ilinked] (type)"+ c.type+"  "+str(c) + " is not a link constraint."
			return False
	
		if self.a == c.b:
			return True
		else:
			return False
	
	def __cmp__(self,o):
		if self.value == o.value and self.type == o.type and self.sameOps(o) :
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

	def __contains__(self, c):	
		if self.a == c or self.b == c :
			return True
		else:
			return False
		

