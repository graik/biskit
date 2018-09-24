import biskit.tools as T

t = T.load( 'com_fake.etraj' )

x = t.takeFrames( range(0, t.n_members * 5) )

x.ref.disconnect()

T.dump( x, 'extract.etraj' )
