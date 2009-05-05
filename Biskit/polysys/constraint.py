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
		"""
		Interval comparison.
		
		@return: True if the two ntervals have the samestart and ending points.
		@rtype: bool
		"""
		
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
		

