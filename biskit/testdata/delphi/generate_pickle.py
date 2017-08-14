import biskit.tools as T
from biskit.exe import Reduce        
from biskit import PDBModel

m1 = PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )

rec = PDBModel( T.testRoot('com/raw/rec.pdb'))
lig = PDBModel( T.testRoot('com/raw/lig.pdb'))

m2 = rec.concat(lig)

x = Reduce( m1, debug=False, verbose=True,
                 autocap=True )
m1 = x.run()

x = Reduce( m2, debug=False, verbose=True,
                 autocap=True )
m2 = x.run()

m1.saveAs( T.testRoot('delphi/1A19_reduced.model'))
m2.saveAs( T.testRoot('delphi/complex_reduced.model'))
