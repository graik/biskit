#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T

m = B.PDBModel('1A2P.pdb')
m = m.takeChains([0])
m = m.sort()

m.saveAs('1A2P_dry.model')
