from FRETUtils import dbPre
from Chromophore import Chromophore
from Protein import Protein

class FRETProtein (Protein):
	def __init__ (self,name="FRETProtein"):
		
		self.name = name
		
		if self.name != "FRETProtein":
			
			self._loadFromDB(name)
			
			self.chromophore = Chromophore(self.name+'_Cromophore',self._source,self)
			
			Protein.__init__(self,name,self._source+'.pdb')
			
			
			if not self.chromophore.defined:
				print "[WARNING FRETProtein __init__] Chromophore not defined. Defaulting." 
				self.name = "mCitrine"
				self.chromophore = Chromophore(self.name+'_Cromophore_Default',self._source,self)
				self.name = name
	
	def __str__(self):
		desc = Protein.__str__(self)
		desc += "\n["+str(self.name)+", lifetime: "+str(self.lifetime)+", abso. wavelenth: "+str(self.abswl)+", emis. wavelenth: "+str(self.emwl)+", quantum yield: "+str(self.quantumYield)+", PDB source file: "+self._source+"]"
		return desc
		
	def onInsertion(self,myassembly=None):
		if self.chromophore.defined:
			chromo = Block(self.chromophore.name,self.chromophore,self)
			chromo.addInterval(self.chromophore.atomrange[0],self.chromophore.atomrange[1])
			myassembly.addBlock(chromo)
		
	def _loadFromDB(self,name):
		f = open ('./fret_prots_db/single_parameters.db',"r")
		lineas = f.readlines()
		f.close()
		line = -1
		
		if len(lineas)<=1:
			print "[WARNING FRETProtein _loadFromDB (__init__)] "+name+" is not an available name in database (empty or bad format database?)." 
			return False
		else:
			for i in range(len(lineas)):
				if name in lineas[i]:
					line = i 
		if line == -1 :
			print "[WARNING FRETProtein _loadFromDB (__init__)] "+name+" is not an available name in database (name not found)." 
			return False
		
		parameters = lineas[line].split('\t',7)
		
		self.lifetime = dbPre(parameters[1],'float')
		self.abswl = dbPre(parameters[2],'float',-1)
		self.emwl = dbPre(parameters[3],'float',-1)
		self.epsilon = dbPre(parameters[4],'float')
		self.quantumYield = dbPre(parameters[5],'float')
		
		self._source = parameters[6]
		if self._source =='X' or self.abswl == -1 or self.emwl == -1:
			print "[ERROR FRETProtein _loadFromDB (__init__)] some mandatory parameters are not defined (X) in DB." 
			return False
		
		if len(parameters)>6:
			self.notes = parameters[7]
		
		return True