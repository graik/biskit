#!/usr/bin/env python

import Biskit.tools as T

## pickle to raw data
t = T.load( 'traj.dat' )

t.ref.writePdb('traj_ref.pdb')
T.dump(t.frameNames, 'traj_framenames.list')  ## for sorting by time and member
t.writeCrd('traj.crd')



