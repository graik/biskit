


from block import BlockEntity

class PolyChromophore ( BlockEntity , Chromophore):
	""""
	Wrapper for Chromophore. Please see Chromophore for further information.
	"""
	def __init__(self,name,source, prot):
	
		BlockEntity.__init__(self,name,prot)
		Chromophore.__init__(self,name,source)
		

class PolyFRETEntity(BlockEntity , FRETEntity):
	
	"""Path for """
	PROTEIN_STRUCTURE_SOURCE_DIR = "./"
	
	""""
	Wrapper for FRETEntity. Please see FRETEntity for further information.
	"""
	def __init__(self,name):
	
		BlockEntity.__init__(self,name,prot)
		
		structure = PDBModel(PROTEIN_STRUCTURE_SOURCE_DIR+source+.pdb")
		
		f.structure = structure
		FRETEntity.__init__name="FRETEntity",database="", structure = None,chromo_autodef = False,verbose = True
		
