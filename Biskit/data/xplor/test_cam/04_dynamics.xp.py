import protocol
import ivm
from xplorPot import XplorPot
import pdbTool

protocol.initTopology(('protein'))
protocol.initParams(('protein'))

protocol.loadPDB("1up5_cleaned.pdb")

xplor.command(\
"""hbuild
   selection=( hydrogen ) 
   phistep=10.0
end""")

pdbTool.PDBTool("01_hbuild.pdb").write()
protocol.addUnknownAtoms()

pdbTool.PDBTool("02_addUnknown.pdb").write()

protocol.fixupCovalentGeom(maxIters=10,useVDW=1)

dyn = ivm.IVM()

dyn.group( AtomSel( 'resid 1:69' ) )
dyn.group( AtomSel( 'resid 84:145') )

dyn.setTrajectory( 'rb_traj.log' )

protocol.torsionTopology( dyn )

from potList import PotList

pot = PotList()
pot.append( XplorPot('CDIH') )
pot.append( XplorPot('DIHE') )
pot.append( XplorPot('VDW') )
pot.append( XplorPot("BOND") )
pot.append( XplorPot("ANGL") )
pot.append( XplorPot("IMPR") )

init_t = 500

protocol.initDynamics( ivm=dyn,
                       potList=pot,
                       bathTemp=init_t,
                       finalTime=5,
                       numSteps=1000,
                       printInterval=10,
                       )
dyn.setETolerance( init_t/100 )
dyn.run()

pdbTool.PDBTool("03_dyn.pdb").write()
