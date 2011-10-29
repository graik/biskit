import Biskit.tools as T
import Biskit as B
import numpy as N

com = T.load( 'com.traj' )

# re-order frames into 4 parallel trajectories
frames = N.zeros( len(com), int )

for i in range( 11 ):
    N.put( frames, range(i*4,i*4+4), N.arange(i,44,11) )

etraj = B.EnsembleTraj( n_members=4 )
etraj.frames = com.takeFrames( frames ).frames
etraj.ref = com.ref
etraj.resetFrameNames()
etraj.ref.disconnect()

# separate protein and DNA into two chains
etraj.ref.chainIndex(breaks=True, force=True, cache=True)
etraj.ref.addChainId()

## extract only some residues for speed
t1 = etraj.takeAtoms( etraj.ref.res2atomIndices(range(10)) )
t2 = etraj.takeChains( [1] )

etraj = t1.concatAtoms( t2 )

T.dump( etraj, 'com_fake.etraj' )
