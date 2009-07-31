import protocol
import ivm
from xplorPot import XplorPot
import pdbTool

protocol.initTopology(('protein'))
protocol.initParams(('protein'))

protocol.loadPDB("03_min.pdb")
xplor.simulation.deleteAtoms("not known")

## protocol.fixupCovalentGeom(maxIters=100,useVDW=1)
protocol.massSetup()

dyn = ivm.IVM()

dyn.group( AtomSel( 'resid 1:107' ) )
dyn.group( AtomSel( 'resid 122:338') )

protocol.torsionTopology( dyn )

from potList import PotList

pot = PotList()
pot.append( XplorPot('CDIH') )

protocol.initRamaDatabase()
pot.append( XplorPot('RAMA') )

pot.append( XplorPot('VDW') )
pot.append( XplorPot("BOND") )
pot.append( XplorPot("ANGL") )
pot.append( XplorPot("IMPR") )

init_t = 500

from monteCarlo import randomizeTorsions
randomizeTorsions(dyn)

pdbTool.PDBTool("05_random.pdb").write()

from simulationTools import AnnealIVM

anneal= AnnealIVM( initTemp = 3000,
                   finalTemp= 25,
                   tempStep = 25,
                   ivm      = dyn,
                   rampedParams = coolParams )
anneal.run()

pdbTool.PDBTool("05_anneal.pdb").write()


protocol.initDynamics( ivm=dyn,
                       potList=pot,
                       bathTemp=init_t,
                       finalTime=1,
                       numSteps=500,
                       printInterval=10,
                       )
dyn.setETolerance( init_t/100 )
dyn.run()

pdbTool.PDBTool("05b_dyn.pdb").write()
