from Biskit import PDBModel
import Biskit.tools as T
import Biskit.molUtils as MU
import numpy as N

from Biskit import PDBCleaner
import Biskit.tools as T

linker = 'ATGSGSGSGSGSGS'
outfile= '01_sequence.txt'

m1 = PDBModel( 'fkbp12.pdb' )
m2 = PDBModel( '1zgq_single.pdb' )

## cleanup structures
m1 = m1.compress( m1.maskProtein() * m1.maskHeavy() )
m2 = m2.compress( m2.maskProtein() * m2.maskHeavy() )

m1.remove( m1.indicesFrom( 'name', ['OXT'] ) )
m2.remove( m2.indicesFrom( 'name', ['OXT'] ) )

m1 = PDBCleaner( m1 ).process()
m2 = PDBCleaner( m2 ).process()


## orient structures

from Biskit.polysys import restools as RT

diff = m1.xyz[-1] - m2.xyz[0]

m2.xyz = m2.xyz + diff - [10,10,10]


## construct sequence
seq = m1.sequence() + linker + m2.sequence()

seq = ' '.join( MU.single2longAA( seq ) )

f = open( outfile, 'w' )
f.write( seq )
f.close()

## adapt residue enumeration to sequence

m1.renumberResidues()
m2.renumberResidues( start= m1.lenResidues() + len(linker) + 1 )

m1['serial'] = range( len(m1) )
m2['serial'] = range( len(m2) )

## equalize segid
m1['segment_id'] = ['p1'] * len( m1 )
m2['segment_id'] = ['p1']  * len( m2 )

## equalize chain ID
m1['chain_id'] = [''] * len( m1 )
m2['chain_id'] = ['']  * len( m2 )


m1.writePdb( '01_domain1.pdb' )
m2.writePdb( '01_domain2.pdb' )

