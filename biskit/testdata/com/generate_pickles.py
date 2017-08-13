#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T

m1 = B.PDBModel('rec.pdb')
m1.saveAs('rec.model')

m2 = B.PDBModel('lig.pdb')
m2.saveAs('lig.model')

