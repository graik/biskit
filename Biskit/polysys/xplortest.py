import protocol
protocol.initTopology(('protein'))
protocol.initParams(('protein'))

protocol.initStruct('res2.psf')
protocol.initCoords('res2.pdb')


#~ from protocol import loadPDB
#~ loadPDB("2PPZ.pdb")

#~ from protocol import addUnknownAtoms
#~ addUnknownAtoms()

from protocol import fixupCovalentGeom
fixupCovalentGeom('all')


#~ from ivm import IVM
#~ integrator = IVM()
#~ integrator.fix( AtomSel("") ) #resid 100:120
#~ integrator.group( AtomSel("resid 130:140") )
#~ from protocol import torsionTopology
#~ torsionTopology(integrator,oTensors=listOfVarTensors)


#~ from monteCarlo import *

#~ mc = MC(simulation.currentSimulation())
#~ randomizeTorsions(mc, sel='all')


import pdbTool
pdb = pdbTool.PDBTool("aux_out.pdb")
pdb.write()

exit()
