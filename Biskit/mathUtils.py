## general purpose math methods
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##

## $Revision$
## last $Author$
## last $Date$


import Numeric as N
import random
import RandomArray
import math

def accumulate( a ):
    """
    cumulative sum of a[0], a[0]+a[1], a[0]+a[1]+[a2], ...
    normalized by N.sum( a )
    """
    return add.accumulate( a ) / N.sum( a )


def variance(x, avg = None):
    if avg is None:
        avg = N.average(x)

    if len(x) == 1:
        return 0.0
    
    return N.sum(N.power(N.array(x) - avg, 2)) / (len(x) - 1.)


def SD(x, avg = None):
    """
    Standard deviation
    """
    return N.sqrt(variance(x, avg))


def wMean(x, w=None):
    """
    Mean of data (x) weighted by (w).
    x - array, X-D array with numbers
    w - array, 1-D array of same length as x with weight factors
    -> array('f') or float
    """
    if w is None:
        wx = x
    else:
        wx = [ x[i] * 1. * w[i] for i in range( len(x) ) ]

    return N.sum(wx)/N.sum(w)


def wVar(x, w):
    """
    Variance of weighted (w) data (x).
    """
    wm = wMean(x,w)
    return ( N.sum(w) / ( (N.sum(w)**2-N.sum(w**2)) ) ) * N.sum(w*(x-wm)**2)


def wSD(x, w):
    """
    Standard deviation of weighted data.
    """
    return N.sqrt( wVar(x, w) )

    
def aboveDiagonal( pw_m ):
    """
    pw_m - symmetric square matrix
    -> raveled list of 'unique' values without diagonal
    """
    r = []
    for i in range( 0, len( pw_m ) ):

        for j in range( i+1, len( pw_m ) ):

            r.append( pw_m[i,j] )

    return r

def arrayEqual( a, b ):
    """
    Compare 2 arrays or lists of numbers for equality.
    a - array / list, (multi-dimensional is supported)
    b -       ~
    -> 1|0, 1 if array/list a equals array/list b
    """
    if a is None or b is None:
        return a is b
        
    if len(a) != len(b):
        return 0

    if type(a) is list:  a = N.array( a )
    if type(b) is list:  b = N.array( b )

    a = N.ravel( a )
    b = N.ravel( b )

    return N.sum( a==b ) == len(a)

   
def pairwiseDistances(u, v):
    """
    -> N.array( len(u) x len(v) ) of double, pw distances between 2 arrays.
    """
    diag1= N.diagonal(N.dot(u,N.transpose(u)))
    diag2= N.diagonal(N.dot(v,N.transpose(v)))
    dist= -N.dot(v,N.transpose(u))-N.transpose(N.dot(u,N.transpose(v)))
    dist= N.transpose(N.asarray(map(lambda column,a:column+a, \
                                N.transpose(dist), diag1)))
    return N.transpose(N.sqrt(N.asarray(
        map(lambda row,a: row+a, dist, diag2))))


def randomMask( nOnes, length ):
    """
    Create random array of given lenght and number of ones.
    nOnes - int, number of ones
    lenght - int, lenght of array
    -> array with ones and zeros
    """
    r = N.zeros( length )
    pos = []
    
    ## add random ones
    for i in range( nOnes ):
        pos += [ int( random.random() * length ) ]      
    N.put( r, pos, 1 ) 
        
    ## if two ones ended up on the same position
    while nOnes != N.sum(r):
        pos = int( random.random() * length )
        N.put( r, pos, 1 )

    return r


def random2DArray( matrix, ranNr=1, mask=None):
    """
    Create randomized 2D array containing ones and zeros.
    matrix - matrix to randomize
    mask - mask
    ranNr - integer, number of matricies to add up
    -> 2D array or |ranNr| added contact matricies
    !! MathUtilError, if mask does not fit matrix
    """
    ## get shape of matrix
    a,b = N.shape( matrix )

    ## get array from matrix that is to be randomized
    if mask:
        if len(mask) == len( N.ravel(matrix) ):
            array = N.compress( mask, N.ravel(matrix) )
        
        if len(mask) != len( N.ravel(matrix) ):
            raise MathUtilError(
                'MatUtils.random2DArray - mask of incorrect length' +
                '\tMatrix length: %i Mask length: %i'\
                  %(len( N.ravel(matrix) ), len(mask)))
            
    if not mask:
        array = N.ravel(matrix)
        
    ## number of ones and length of array
    nOnes = int( N.sum( array ) )
    lenArray = len( array )
    ranArray = N.zeros( lenArray )
    
    ## create random array
    for n in range(ranNr):
        ranArray += randomMask( nOnes, lenArray )
        
    ## blow up to size of original matix
    if mask:
        r = N.zeros(a*b)
        N.put( r, N.nonzero(mask), ranArray)
        return N.reshape( r, (a,b) )
    
    if not mask:
        return  N.reshape( ranArray, (a,b) )


def runningAverage( x, interval=2, preserve_boundaries=0 ):
    """
    x        - list of int/float
    interval - window size (-(interval-1)/2 to +(interval-1)/2)
    preserve_boundaries - 0||1, shrink window at edges to keep original
                          start and end value
    -> [ float ]
    """

    if interval == 0:
        return x
    
    l = []

    interval = int((interval-1)/2)

    if not preserve_boundaries:

        for i in range(len(x)):

            left = max(0, i - interval)
            right = min(len(x), i + interval + 1)

            slice = x[left:right]

            l.append(N.average(slice))
    else:

        for i in range( len(x) ):

            left = i - interval
            right= i + interval + 1

            if left < 0:
                right = right + left
                left = 0
            if right > len(x):
                left = left + right - len(x)
                right = len(x)

            slice = x[left:right]

            l.append(N.average(slice))

    return N.array(l)


def packBinaryMatrix( cm ):
    """
    Compress sparse array of 0 and ones to list of one-positions.
    cm - X by Y array of int
    -> {'shape':(X,Y), 'nonzero':[int] }
    """
    if cm == None or type( cm ) == dict:
        return cm
    
    result = {}
    result['shape'] = N.shape( cm )
    result['nonzero'] = N.nonzero( N.ravel( cm ) )
    result['nonzero'] = result['nonzero'].tolist()
    return result


def unpackBinaryMatrix( pcm, raveled=0 ):
    """
    Uncompress array of 0 and 1 that was compressed with packBinaryMatrix.
    pcm - {'shape':(X,Y,..), 'nonzero':[int]}
    raveled - 1|0, return raveled (default 0)
    -> N.array(X by Y by ..) of int
    """
    if type( pcm ) != dict:
        return pcm

    s = pcm['shape']

    m = N.zeros( N.cumproduct( s )[-1], 'i')
    m.savespace( 1 )
    N.put( m, pcm['nonzero'], 1 )

    if raveled:
        return m

    return N.reshape( m, s )


def matrixToList( cm ):
    """
    Convert matrix into standard python list remembering the dimensions..
    Not used right now.
    -> {'shape':(int,..), 'lst':[..] }
    """
    if cm == None or type( cm ) == dict:
        return cm

    result = {}
    result['shape'] = N.shape( cm )
    result['lst'] = N.ravel( cm ).tolist()

    return result


def listToMatrix( lcm ):
    """
    Convert result of matrixToList back into Numeric array
    Not used right now.
    -> Numeric.array
    """
    if type( lcm ) != dict:
        return lcm

    s = lcm['shape']
    return N.reshape( lcm['lst'], s )
    

def eulerRotation(alpha, beta, gamma):
    """
    Builds a rotation matrix from successive rotations:
    1. rotation about z-axis by angle alpha
    2. rotation about x-axis by angle beta
    3. rotation about z-axis by angle gamma
    -> 3 x 3 array of float
    Michael Habeck
    """
    cos_alpha = N.cos(alpha); sin_alpha = N.sin(alpha)
    cos_beta  = N.cos(beta);  sin_beta  = N.sin(beta)
    cos_gamma = N.cos(gamma); sin_gamma = N.sin(gamma)

    R = N.zeros((3,3), N.Float)

    R[0][0] = cos_gamma * cos_alpha - sin_gamma * cos_beta * sin_alpha
    R[0][1] = cos_gamma * sin_alpha + sin_gamma * cos_beta * cos_alpha
    R[0][2] = sin_gamma * sin_beta

    R[1][0] = -sin_gamma * cos_alpha - cos_gamma * cos_beta * sin_alpha
    R[1][1] = -sin_gamma * sin_alpha + cos_gamma * cos_beta * cos_alpha
    R[1][2] =  cos_gamma * sin_beta
  
    R[2][0] =  sin_beta * sin_alpha
    R[2][1] = -sin_beta * cos_alpha
    R[2][2] =  cos_beta

    return R

def randomRotation():
    """
    Get random rotation matrix.
    -> 3 x 3 array of float
    Michael Habeck
    """
    alpha = RandomArray.random() * 2 * N.pi
    gamma = RandomArray.random() * 2 * N.pi
    beta  = N.arccos(2*(RandomArray.random() - 0.5))

    return eulerRotation(alpha, beta, gamma)

def intersection( a, b, optimize=1 ):
    """
    intersection( [any], [any] ) -> [any], intersection of the two lists
    optimize - 1|0, result is sorted like the shorter of the two lists [1]
    """
    if optimize and len(a) > len(b):
        a, b = b, a

    return [ x for x in a if x in b ]

def nonredundant( l ):
    """
    noredundant( [any] ) -> [any], all members of l without duplicates
    """
    r = []
    for x in l:
        if x not in r:
            r.append( x )
    return r

def union( a, b ):
    """
    union( [any], [any] ) -> [any], union of two lists (without duplicates)
    """
    if type( a ) is N.arraytype:
        a = a.tolist()
    if type( b ) is N.arraytype:
        b = b.tolist()
        
    return nonredundant( a + b )

def difference( a, b ):
    """
    difference([any], [any]) -> [any], all members of a that are not in b
    """
    return [ x for x in a if x not in b ]

def removeFromList( l, v, all=1 ):
    """
    removeFromList( l, v, [all=1] )
    remove all or first occurrence(s) of v from l.
    l   - [ any ]
    v   - any or [ any ]
    all - 0|1, remove all occurrences (1) or only first one (0) [1]
    """
    if type( v ) != type( [] ):
        v = [ v ]

    for x in v:
        
        while x in l:
            l.remove( x )
            if not all:
                break

def randomRange( start, stop, n ):
    """
    start - minimal index
    stop  - 1+maximal index
    n     - number of indices
    -> [int] set of unique integers evenly distributed between start and stop
    """
    r = []
    while len( r ) < n:
        x = int( round( random.uniform( start, stop ) ) )
        if x < stop and not x in r:
            r += [x]
    r.sort()
    return r
                                        
def linfit( x, y ):
    """
    Calculate linear least-square fit to the points given by x and y.
    see http://mathworld.wolfram.com/LeastSquaresFitting.html
    x - [ number ]
    y - [ number ]
    -> float, float, float - m, n, r^2 (slope, intersection, corr. coefficient)
    !! BiskitError, if x and y have different number of elements
    """
    x, y = N.array( x, 'd'), N.array( y, 'd')
    if len( x ) != len( y ):
        raise Exception, 'linfit: x and y must have same length'

    av_x = N.average( x )
    av_y = N.average( y )
    n = len( x )

    ss_xy = N.sum( x * y ) - n * av_x * av_y
    ss_xx = N.sum( x * x ) - n * av_x * av_x
    ss_yy = N.sum( y * y ) - n * av_y * av_y

    slope = ss_xy / ss_xx

    inter = av_y - slope * av_x

    corr  = ss_xy**2 / ( ss_xx * ss_yy )

    return slope, inter, corr
    

def cartesianToPolar( xyz ):
    """
    Convert cartesian coordinate array to polar coordinate array
    xyz - array of cartesian coordinates (x, y, z)
    -> array of polar coordinates (r, theta, phi)
    """
    r = N.sqrt( N.sum( xyz**2, 1 ) )
    p = N.arccos( xyz[:,2] / r )

    ## have to take care of that we end up in the correct quadrant
    t=[]
    for i in range(len(xyz)):
        ## for theta (arctan)
        t += [math.atan2( xyz[i,1], xyz[i,0] )]

    return N.transpose( N.concatenate( ([r],[t],[p]) ) )
    

def polarToCartesian( rtp ):
    """
    Convert polar coordinate array to cartesian coordinate array
    rtp - array of cartesian coordinates (r, theta, phi)
    -> array of polar coordinates (x, y, z)
    """
    x = rtp[:,0] * N.cos( rtp[:,1] ) * N.sin( rtp[:,2] )
    y = rtp[:,0] * N.sin( rtp[:,1] ) * N.sin( rtp[:,2] )
    z = rtp[:,0] * N.cos( rtp[:,2] )
    
    return N.transpose( N.concatenate( ([x],[y],[z]) ) )


def projectOnSphere( xyz, radius=None, center=None ):
    """
    xyz    - array N x 3 of float, cartesian coordinates
    radius - float, radius of target sphere [maximal distance to center]
    center - array 0 x 3 of float, center of the sphere [average(xyz)]
    """
    if center is None:
        center = N.average( xyz )

    if radius is None:
        radius = max( N.sqrt( N.sum( N.power( xyz - center, 2 ), 1 ) ) )

    rtp = cartesianToPolar( xyz - center )
    rtp[ :, 0 ] = radius

    return polarToCartesian( rtp ) + center
