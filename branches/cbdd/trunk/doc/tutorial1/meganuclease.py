from Biskit import *

## load a file
m0 = PDBModel('1R7M.pdb')
m0.sequence()

m0.report()

## explore profiles
m0.atoms
atom = m0.atoms[-1]

## explore atoms
m0[-1]

## access profiles
t = m0['temperature_factor']

import Biskit.gnuplot as G
G.plot( t )

#-------------------------------------

## explore residue profiles
m0.residues

## create profiles
m0['aa'] = m0.sequence()


m0.residues['rname'] = m0.atom2resProfile( 'residue_name' )
# is same as...
import numpy as N
# m0['rname'] = N.take( m0['residue_name'], m0.resIndex() )

#-------------------------------------
#
### Now let's do something useful!
#

## extract single copy
m0.report()
m = m0.takeChains( [0, 2, 3 ] )

## quick and dirty 2D-view
# G.plot( zip( m.xyz[:,0], m.xyz[:,1] ) )

m.gnuplot()
m0.report( plot=1 )

#-------------------------------------
## distance/contact matrix
import Biskit.mathUtils as MU

pw = MU.pairwiseDistances( m.xyz, m.xyz )

## nicer graphics with Matrixplot

p = MatrixPlot( pw, step=25 )
p.show()

## simple binary contact plot with gnuplot

## cont = pw * N.less( pw, 10 )
cont = N.less( pw, 10 )
G.scatter( zip( *N.nonzero( cont ) ) )

#
## inter-molecular contacts only
#

# separate DNA and protein

mask_dna = [ name in ['A','C','T','G']
             for name in m['residue_name'] ]

dna = m.compress( mask_dna )
dna.writePdb( 'dna_only.pdb' )

prot = m.compress( m.maskProtein() )

## Contact analysis with Complex
from Biskit.Dock import Complex

com = Complex( prot, dna )

## what is Complex??
dir( com )
print com.atomContacts.__doc__

cont = com.atomContacts( cutoff=10 )
print cont

G.scatter( zip( *N.nonzero( cont ) ) )
# --------------------------
##
## Visualization
##

# how many contacs has each protein residue?
com.lig_model['contacts'] = N.sum( cont, axis=0 )
com.rec_model['contacts'] = N.sum( cont, axis=1 )

pm = Pymoler()
pm.addPdb( com.rec(), 'rec' )
pm.addPdb( com.lig(), 'lig' )

pm.colorAtoms( 'rec', com.rec()['contacts'] )
pm.colorAtoms( 'lig', com.lig(force=1)['contacts'] )

pm.add( 'show surface')
pm.run()

##
## Add external data
##

d = PDBDope( com.rec_model )

# d.addConservation()
prot = T.load('rec.model')

prot.residues

pm = Pymoler()
pm.addPdb( prot, 'rec' )
pm.addPdb( com.lig(force=1), 'lig')
pm.colorRes( 'rec', prot['cons_ent'] )
pm.add( 'show surface')
pm.run()

##
## Correlate
##

prot['rescontacts'] = prot.atom2resProfile( 'contacts', f=N.sum )

G.scatter( zip( prot['rescontacts'], prot['cons_ent'] ) )

p = prot.residues.plot( 'rescontacts', 'cons_ent' )
p.show()

_cont = prot['rescontacts'] / 1. / N.max( prot['rescontacts'] )
_cons = prot['cons_ent'] / 1. / N.max( prot['cons_ent'] )

G.plot( _cont, _cons )

