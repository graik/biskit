## convert renamed classes
## this requires the classes to be available under both the old and the new name
## otherwise unpickling will break

import biskit

f = open('traj.dat', 'rb')
this = pickle.load(f, encoding='latin1')
this.__class__ = biskit.EnsembleTraj

## all the following is untested, safer to re-create from raw data
this.profiles.__class__ = biskit.ProfileCollection

this.ref.__class__ = biskit.PDBModel
this.ref.source.__class__ = biskit.LocalPath

this.atoms.__class__ = biskit.pdbModel.PDBProfiles
this.residues.__class__ = biskit.pdbModel.PDBResidueProfiles
