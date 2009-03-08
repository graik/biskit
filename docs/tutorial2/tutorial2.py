import numpy as N

from Biskit import *
from Biskit.gnuplot import *


#########################
## Handling Structures I
## load & inspect
#########################

# loading and quick overview:
m = PDBModel( '1r7m.pdb' )
## m = PDBModel('1R7M')  ## fetch from database

m.sequence()
m.report()
m.plot()  ## quick and dirty plot

len( m )  ## == m.lenAtoms()

# coordinates

m.xyz

# infos

m.info

# atom profiles:

m.atoms
plot( m.atoms['temperature_factor'] )
plot( m['temperature_factor'] )

# residue profiles .. empty to start with

m.residues
m.residues['aa'] = m.sequence()  ## create a residue profile

# convert atom profile to residue profile, average over atoms in residues
m['resT']     = m.atom2resProfile( 'temperature_factor', N.average )
# convert again but keep maximum value instead
m['resT_max'] = m.atom2resProfile( 'temperature_factor', max )

m.residues
p = m.residues.plot( 'resT', 'resT_max' )
p.show()

#########################
## numpy recap
#########################

a = N.arange( 100 )

mask = a % 2 == 0
c = N.compress( mask, a )

indices = N.flatnonzero( mask )

b = N.take( a, indices )

print N.alltrue( b == c )


#########################
## Chain handling
#########################

m.report()

## kick out the second copy of protein/DNA
m = m.takeChains( [0, 2, 3 ] )

m.plot()


#########################
## Adding Data
#########################

d = PDBDope( m )

d.addSecondaryStructure()

d.addSurfaceRacer()


#########################
## Plotting / Visual
#########################

pm = Pymoler()
pm.addPdb( m, 'prot' )
pm.colorAtoms( 'prot', m['AS'] )
pm.run()


#########################
## Handling Structures II
## select & compress
#########################

m.report()

mask = m.maskProtein()
prot = m.compress( mask )

dna  = m.compress( m.maskDNA() )

# repeat AS calculation only on protein

## backup old AS values into a new profile
prot['AS_bound'] = prot['AS']

## recalculate
d = PDBDope( prot )
d.addSurfaceRacer()

prot['AS_free'] = prot['AS']
prot['AS_delta']= prot['AS_free'] - prot['AS_bound']

## color complex by delta-ASA (free->bound)
pm = Pymoler()
pm.addPdb( prot, 'prot' )
pm.addPdb( dna, 'dna')
pm.colorAtoms( 'prot', prot['AS_delta'] )
pm.run()


#########################
## Saving and Pickling
#########################
import Biskit.tools as T

prot.saveAs( 'prot.model' )
m2 = T.Load( 'prot.model' )

prot.writePdb( 'prot.pdb' )
m2 = PDBModel( 'prot.pdb' )


#########################
## Trajectories
#########################

## see md tutorial!

#########################
## Summary Example
#########################

m = PDBModel( '3TGI' )
m = m.compress( m.maskProtein() )

diff = m.xyz - m.centerOfMass()

dist = N.sqrt( N.sum( diff**2, axis=1 ) )

# plot

mask = dist < 10.

m_center = m.compress( mask )

pm = Pymoler()
pm.addPdb( m_center, 'center' )
pm.addPdb( m, '3tgi' )
pm.colorAtoms( '3tgi', dist )
pm.show()




