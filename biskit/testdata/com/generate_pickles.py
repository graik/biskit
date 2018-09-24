#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T

m1 = B.PDBModel('raw/rec.pdb')
m1.saveAs('rec.model')

m2 = B.PDBModel('raw/lig.pdb')
m2.saveAs('lig.model')

## not tested
import biskit.dock as D

com = D.Complex(m1, m2)
T.dump(com, 'ref.complex')
