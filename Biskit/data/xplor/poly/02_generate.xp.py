import psfGen, protocol, pdbTool
import random as rand

simWorld.setRandomSeed(5521)

import protocol
protocol.initTopology(('protein'))
protocol.initParams(('protein'))

seq = open('01_sequence.txt').read()

psfGen.seqToPSF(seq,seqType='prot',startResid=1, segName='p1')

xplor.command( 'vector do (x = 0) (all)' )
xplor.command( 'vector do (y = %4.2f) (all)' % rand.random() )
xplor.command( 'vector do (z = %4.2f) (all)' % rand.random() )

## protocol.genExtendedStructure(maxFixupIters=50)

## read structures
protocol.initCoords( ['01_domain1.pdb', '01_domain2.pdb'] )

xplor.command(\
"""hbuild
   selection=( hydrogen and (resid 1:107 or resid 122:338) ) 
   phistep=10.0
end""")

protocol.fixupCovalentGeom(sel=AtomSel('(attr x=0.0) and resid 1:107'),
                           maxIters=20, verbose=1 )
protocol.fixupCovalentGeom(sel=AtomSel('(attr x=0.0) and resid 122:338'),
                           maxIters=20, verbose=1 )

## fix linker
xplor.command( 'vector do (x = %4.2f) (resid 108:121)' % rand.random())
xplor.command( 'vector do (y = %4.2f) (resid 108:121)' % rand.random())
xplor.command( 'vector do (z = %4.2f) (resid 108:121)' % rand.random())

protocol.fixupCovalentGeom(sel=AtomSel('resid 108:121'), maxIters=20)

## restrain domain 1 and domain 2 to themselves but allow them to move

xplor.command(r"""
flags 
   include ncs 
end
ncs restraints
  initialize
   group                          ! make another group for other segment
         equi ( (resid 1:107) )   !fix these two lines.
!         equi ( (resid 1:107) )
         weight 500.
   end
   group
         equi ((resid 122:338) )   !fix these two lines.
!         equi ((resid 122:338) )
         weight 500.
   end
   ?   {* print the NCS relations when starting *}
end
""")


xplor.command( r"""minimize powell nstep=500 nprint 5 end""")

xplor.command(r"""
flags 
   exclude ncs 
end
""")

xplor.command("write psf output=03_system.psf end")
pdbTool.PDBTool("03_min.pdb").write()
