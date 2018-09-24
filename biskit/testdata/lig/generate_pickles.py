#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T

m = B.PDBModel('1A19.pdb')
m = m.compress( m.maskProtein() )
m.saveAs('1A19_dry.model')
