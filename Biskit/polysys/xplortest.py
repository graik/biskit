from psfGen import pdbToPSF
from regularize import addUnknownAtoms
from pdbTool import PDBTool

pdbToPSF("aux_in.pdb")
addUnknownAtoms()
pdb = PDBTool("aux_out.pdb")
pdb.write()
exit()
