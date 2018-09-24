
from Biskit import *
import Biskit.tools as T
import numpy as N

###############################
## Creating a Trajectory
## Example 1:
## starting from several single structures
###############################

m = PDBModel( '3TGI' )

t = Trajectory( 20 * [ m ] )  ## trajectory of 20 identical structures

t.ref          ## reference PDBModel
t.lenFrames()  ## number of frames, same as len(t)
t.lenAtoms()   ## number of atoms,  same as len( t.ref )

## kick out non-protein atoms (ions, solvent, etc.)
t = t.compressAtoms( t.ref.maskProtein() )

## shift each frame by incrementing delta-x
for i in range( len( t ) ):
    t.frames[i] += [ i, 0, 0 ]

pm = Pymoler()
## t[i] returns a PDBModel instance for frame i
pm.addMovie( [ t[i] for i in range(len(t)) ], 'traj' )
pm.add( 'mplay' )
pm.show()

###################################
## Load and split a real trajectory
###################################

## converted from amber crd using:
## amber2traj.py -i sim.crd -r frame0.pdb
t = T.load( 'traj_0.dat' )

## kick out frames: take every 4th frame only
t_short = t.takeFrames( range(0, len(t), 4) )
len( t_short )

## split system into three trajectories ...
## ... containing *roughly* one spectrin repeat each
third = t.lenAtoms() / 3

t0 = t.takeAtoms( range(0, third) )
t1 = t.takeAtoms( range(third, 2*third) )
t2 = t.takeAtoms( range(2*third, 3*third) )

###################################
## RMS fits and plotting
###################################

## fit trajectory to average
t0.fit()
## fit to average using CA only
t0.fit( t0.ref.maskCA(), prof='rms_ca' )
## iterative fit to CA, kicks out outlier regions
t0.fit( t0.ref.maskCA(), prof='rms_ca_it', n_it=5 )

## plot the rms versus frame for the 3 above fits
## (uses ProfileCollection.plot -- also available for ...
## ... PDBModel.atoms, .residues & Trajectory.profiles)
p = t0.profiles.plot( 'rms', 'rms_ca', 'rms_ca_it' )
p.show()


##################################
## divide, fit, and re-connect MD
##################################

## t0 and t1 have same core sequence but different head and tail peptides... 
## align residue/atom content of first and second trajectory
i0, i1 = t0.ref.compareAtoms( t1.ref )

t0 = t0.takeAtoms( i0 )
t1 = t1.takeAtoms( i1 )
## now t0 and t1 have exactly the same atom content and order

t0.fit( n_it=3 )
t1.fit( ref=t0.avgModel(), n_it=3 )  ## fit t1 to average of t0

t0.ref.writePdb( 't0.pdb' )  ## write ref PDB
t0.writeCrd( 't0.crd' )      ## write Amber CRD
t1.writeCrd( 't1.crd' )

t_01 = t0.concat( t1 )       ## concat t0 and t1 in time

p = t_01.profiles.plot( 'rms' )
p.show()  ## repeat the show if the plot gets hidden behind other windows
## Note: do not close the xplot window during the python session!
p.title = 'RMS versus t'
p.write_eps( 'plot.eps' )  ## see biggles documentation

## t_01.ref.writePdb( 't_0_1.pdb' )
## t_01.writeCrd( 't_0_1.crd')
