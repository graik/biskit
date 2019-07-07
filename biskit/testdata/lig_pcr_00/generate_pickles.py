#!/usr/bin/env python
## re-generate binary test data in this folder

import biskit as B
import biskit.tools as T
from biskit.md import Trajectory, EnsembleTraj, traj2ensemble

## create standard trajectory object
t = Trajectory('raw/traj.crd', 'raw/traj_ref.pdb' )
t.frameNames = T.load('raw/traj_framenames.list')

te = traj2ensemble(t, members=10)
te.fit(fit=0)  ## re-calculate profile 'rms' (all-atom fit to average structure)

T.dump(te, 'traj.dat')


