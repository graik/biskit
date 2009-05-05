from numpy import *
import numpy.random as R
import numpy as N

from Biskit.gnuplot import *

########################################
## creating arrays
## [array, arange, zeros]
########################################

a = array( [ [ 1,2,3,4,5 ], [6,7,8,9,10] ] )

#b = array( [ '1', '2', '3' ], int )
b = array( range(10) )
b = arange( 10 )

d = zeros( (10,3), float )
e = ones ( (10,3), float )

#########################################
## Shape, axis, indexing
## [shape, reshape, ravel]
#########################################

print len( e )
print shape( e )

print 'all: ',    sum( e )
print 'axis 0: ', sum( e, axis=0 )
print 'axis 1: ', sum( e, axis=1 )

a = arange( 30 )
a = reshape( a, (10,3) )

print a

# flatten things out again:
print ravel( a )


#########################################
## slicing and indexing
#########################################

print arange( 0, 10, 2 )
print b[0:10:2]

print b[ arange(0,5) ]
print b[0:5]

print a
print a[0]
print a[-1]
print a[0 , 0]
print a[-1,-2]
print a[: , 1]
print a[0 , :]

#########################################
## assigning values
#########################################

a[0,0] = -1
a[:,2] = 0
a[:,1] *= 10
a[-1, :] = [ 1, 2, 3]

print a

#########################################
## Random arrays
## [numpy.random]
#########################################

a = R.poisson( 50., (500, 2) )
b = R.poisson( 30., (500, 2) )

scatter( a, b )

center_a = average( a, axis=0 )
center_b = average( b, axis=0 )

scatter( a, b, [center_a, center_b] )

#########################################
## Operations are element-wise
#########################################

b = b + [0,30]
scatter( a, b, [center_a, center_b] )

# center an array:
b = b - average( b, axis=0 )
b = b + [30, 0]

# clip outliers
scatter( a, clip( b, 0, 99) )


#########################################
## example -- detect values above/below
#########################################

underflow = False

for i in b:
    for j in i:
        if j < 0:
            underflow = True

if underflow:
    print 'underperformer detected!'

# better:

if any( b < 0 ): ## speedup of factor 740
    print 'underperformer detected!'

# see also: less, greater, greater_equal, less_equal

#######################
## pairwise distances
#######################

b += [0, 25]

pw = zeros( (len(a), len(b)), float )
pw.shape

for i in arange( len(pw) ):
    pw[i, :] = sqrt( sum( (a - b[i])**2, axis=1) )

h = histogram( ravel( pw ), 100 )
plot( h[0] )

#######################
## masks
#######################

mask_a = sum( pw < 10, axis=0 )

a_contact = compress( mask_a, a, axis=0 )

scatter( a, b, a_contact )


#########################
## indices & masks
## take, flatnonzero, put
#########################

a = arange( 100 )

mask = a % 2 == 0
c = compress( mask, a )


indices = flatnonzero( mask )
indices = where( mask )

b = take( a, indices )

print alltrue( b == c )

put( a, indices, -1 )

R.shuffle( a )
plot( a )

# sorting

i = argsort( a )
a = take( a, i )

scatter( a )

# mask to indices and back to mask:

i = where( a > 10 )

mask = zeros( len(a), bool )
put( mask, i, 1 ) 


####################
## Structure example
####################

from Biskit import *

m = PDBModel( '3TGI' )
m = m.compress( m.maskProtein() )

diff = m.xyz - m.centerOfMass()

# or by hand:
diff = m.xyz - average( m.xyz, axis=0 )
# diff = m.xyz - average( m.xyz, axis=0, weights=m.masses() )

dist = sqrt( sum( diff**2, axis=1 ) )

# plot

mask = dist < 10.  # or less( dist, 10 )

m_center = m.compress( mask )

pm = Pymoler()
pm.addPdb( m_center, 'center' )
pm.addPdb( m, '3tgi' )
pm.show()

