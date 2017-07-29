import biskit

f = open('traj.dat', 'rb')
this = pickle.load(f, encoding='latin1')
this.__class__ = biskit.EnsembleTraj.EnsembleTraj
