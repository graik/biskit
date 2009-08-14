import protocol
import ivm
from xplorPot import XplorPot
import pdbTool

protocol.initTopology(('protein'))
protocol.initParams(('protein'))

protocol.loadPDB( '1up5_cleaned.pdb' )
    
xplor.command(\
"""hbuild
   selection=( hydrogen ) 
   phistep=10.0
end""")


protocol.addUnknownAtoms()

protocol.fixupCovalentGeom(maxIters=10,useVDW=1)

pdbTool.PDBTool( 'input_protein.pdb' ).write()
xplor.command("write psf output=input_protein.psf end")


dyn = ivm.IVM()

## define fixed parts
dyn.group( AtomSel( 'resid 1:69' ) )
dyn.group( AtomSel( 'resid 84:145') )

import trajFile

t = trajFile.TrajFile( 'rb_traj.crd' )
dyn.setTrajectory( t )

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

## randomization

from monteCarlo import randomizeTorsions
randomizeTorsions(dyn)

pdbTool.PDBTool("05_random.pdb").write()


## try "normal" Torsion dynamics

protocol.initDynamics( ivm=dyn,
                       potList=pot,
                       bathTemp=init_t,
                       finalTime=10,
                       numSteps=5000,
                       printInterval=10,
                       )
dyn.setETolerance( init_t/100 )
dyn.run()

pdbTool.PDBTool("03_dyn.pdb").write()
