#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T

m1 = B.PDBModel('1huy_citrine.pdb')
m1.saveAs('1huy_citrine.model')

m2 = B.PDBModel('1zgp_dsred_dimer.pdb')
m2.saveAs('1zgp_dsred_dimer.model')
