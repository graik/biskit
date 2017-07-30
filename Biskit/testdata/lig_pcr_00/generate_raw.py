#!/usr/bin/env python

import Biskit.tools as T

## pickle to raw data
t = T.load( 'traj.dat' )

t.ref.writePdb('ref.pdb')
t.writeCrd('traj.crd')

