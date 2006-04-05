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

"""
general purpose math methods
"""

import Numeric as N
import random
import RandomArray
import math


def accumulate( a ):
    """
    cumulative sum of C{ a[0], a[0]+a[1], a[0]+a[1]+[a2], ... }
    normalized by C{ N.sum( a ) }

    @param a: array('f') or float
    @type  a: array

    @return: float
    @rtype: float
    """
    return N.add.accumulate( a ) / N.sum( a )


def variance(x, avg = None):
    """
    Variance, S{sigma}^2
    
    @param x: data
    @type  x: array('f') or float
    @param avg: use this average, otherwise calculated from x
    @type  avg: float OR None
    
    @return: float
    @rtype: float    
    """
    if avg is None:
        avg = N.average(x)

    if len(x) == 1:
        return 0.0

    return N.sum(N.power(N.array(x) - avg, 2)) / (len(x) - 1.)


def SD(x, avg = None):
    """
    Standard deviation, S{sigma}

    @param x: data
    @type  x: array('f') or float
    @param avg: use this average, otherwise calculated from x
    @type  avg: float OR None
    
    @return: float
    @rtype: float        
    """
    return N.sqrt(variance(x, avg))


def wMean(x, w=None):
    """
    Weighted mean: Mean of data (x) weighted by (w).
    
    @param x: X-D array with numbers
    @type  x: array
    @param w: 1-D array of same length as x with weight factors
    @type  w: array
    
    @return: array('f') or float
    @rtype: array('f') or float
    """
    if w is None:
        wx = x
    else:
        wx = [ x[i] * 1. * w[i] for i in range( len(x) ) ]

    return N.sum(wx)/N.sum(w)


def wVar(x, w):
    """
    Variance of weighted (w) data (x).

    @param x: X-D array with numbers
    @type  x: array
    @param w: 1-D array of same length as x with weight factors
    @type  w: array
    
    @return: array('f') or float
    @rtype: array('f') or float    
    """
    wm = wMean(x,w)
    return ( N.sum(w) / ( (N.sum(w)**2-N.sum(w**2)) ) ) * N.sum(w*(x-wm)**2)


def wSD(x, w):
    """
    Standard deviation of weighted data.

    @param x: X-D array with numbers
    @type  x: array
    @param w: 1-D array of same length as x with weight factors
    @type  w: array
    
    @return: array('f') or float
    @rtype: array('f') or float     
    """
    return N.sqrt( wVar(x, w) )


def cross(v, w):
    """
    Cross product between two coordinates.

    @param v: coordinate vector
    @type  v: list or array
    @param w: coordinate vector
    @type  w: list or array
    
    @return: array('f') or float
    @rtype: array('f') or float  
    """
    x = v[1]*w[2] - v[2]*w[1]
    y = v[2]*w[0] - v[0]*w[2]
    z = v[0]*w[1] - v[1]*w[0]
 
    return N.array( (x, y, z ) )


def aboveDiagonal( pw_m ):
    """
    Collect att the values above the diagonal in a square
    matrix.
    
    @param pw_m: symmetric square matrix
    @type  pw_m: 2-D array
    
    @return: raveled list of 'unique' values without diagonal
    @rtype: list
    """
    r = []
    for i in range( 0, len( pw_m ) ):

        for j in range( i+1, len( pw_m ) ):

            r.append( pw_m[i,j] )

    return r


def arrayEqual( a, b ):
    """
    Compare 2 arrays or lists of numbers for equality.
    
    @param a: first array (multi-dimensional is supported)
    @type  a: array / list
    @param b: second array (multi-dimensional is supported)
    @type  b: array / list
    
    @return: 1 if array/list a equals array/list b
    @rtype: 1|0
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
    Pairwise distances between two arrays.
    
    @param u: first array 
    @type  u: array
    @param v: second array 
    @type  v: array
    
    @return: Numeric.array( len(u) x len(v) ) of double
    @rtype: array
    """
    diag1 = N.diagonal( N.dot( u, N.transpose(u) ) )
    diag2 = N.diagonal( N.dot( v, N.transpose(v) ) )
    dist = -N.dot( v,N.transpose(u) )\
           -N.transpose( N.dot( u, N.transpose(v) ) )
    dist = N.transpose( N.asarray( map( lambda column,a:column+a, \
                                        N.transpose(dist), diag1) ) )
    return N.transpose( N.sqrt( N.asarray(
        map( lambda row,a: row+a, dist, diag2 ) ) ))


def randomMask( nOnes, length ):
    """
    Create random array of given lenght and number of ones.
    
    @param nOnes: number of ones
    @type  nOnes: int
    @param length: lenght of array
    @type  length: int
    
    @return: array with ones and zeros
    @rtype: array( 1|0 )
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
    
    @param matrix: matrix to randomize
    @type  matrix: 2D array
    @param mask: mask OR None (default: None)
    @type  mask: list(1|0)
    @param ranNr: number of matricies to add up (default: 1)
    @type  ranNr: integer
    
    @return: 2D array or |ranNr| added contact matricies
    @rtype:2D array
    
    @raise MathUtilError: if mask does not fit matrix
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
    Running average (smoothing) over a given data window.
    
    @param x: data
    @type  x: list of int/float
    @param interval: window size C{ (-(interval-1)/2 to +(interval-1)/2) }
                     (default: 2)
    @type  interval: int
    @param preserve_boundaries: shrink window at edges to keep original
                                start and end value (default: 0)
    @type  preserve_boundaries: 0|1
    
    @return: list of floats
    @rtype: [ float ]
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
    Compress sparse array of 0 and ones to list of one-positions
    (space saving function, upack with L{unpackBinaryMatrix}).
    
    @param cm: X by Y array of int
    @type  cm: 2D array 
    
    @return: {'shape':(X,Y), 'nonzero':[int] }
    @rtype: dict
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
    Uncompress array of 0 and 1 that was compressed with L{packBinaryMatrix}.
    
    @param pcm: {'shape':(X,Y,..), 'nonzero':[int]}
    @type  pcm: dict
    @param raveled: return raveled (default: 0)
    @type  raveled: 1|0
    
    @return: N.array(X by Y by ..) of int
    @rtype: 2D array
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
    Convert matrix into standard python list remembering the dimensions.
    Unpack with L{listToMatrix}.
    
    @note: Not used right now.

    @param cm: array of int
    @type  cm: 2D array
    
    @return: {'shape':(int,..), 'lst':[..] }
    @rtype: dict
    """
    if cm == None or type( cm ) == dict:
        return cm

    result = {}
    result['shape'] = N.shape( cm )
    result['lst'] = N.ravel( cm ).tolist()

    return result


def listToMatrix( lcm ):
    """
    Convert result of L{matrixToList} back into Numeric array
    
    @note: Not used right now.

    @param lcm: {'shape':(int,..), 'lst':[..] }
    @type  lcm: dict    
    
    @return: Numeric.array
    @rtype: 
    """
    if type( lcm ) != dict:
        return lcm

    s = lcm['shape']
    return N.reshape( lcm['lst'], s )


def eulerRotation(alpha, beta, gamma):
    """
    Builds a rotation matrix from successive rotations:
      1. rotation about y-axis by angle alpha
      2. rotation about x-axis by angle beta
      3. rotation about z-axis by angle gamma

    @author: Michael Habeck

    @param alpha: euler angle S{alpha}
    @type  alpha: float
    @param beta: euler angle S{beta}
    @type  beta: float
    @param gamma: euler angle S{gamma}
    @type  gamma: float
    
    @return: 3 x 3 array of float
    @rtype: array
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

    @author: Michael Habeck
    
    @return: 3 x 3 array of float
    @rtype: array
    """
    alpha = RandomArray.random() * 2 * N.pi
    gamma = RandomArray.random() * 2 * N.pi
    beta  = N.arccos(2*(RandomArray.random() - 0.5))

    return eulerRotation(alpha, beta, gamma)


def intersection( a, b, optimize=1 ):
    """
    Intersection of the two lists (i.e. all values common to both two lists).
    C{ intersection( [any], [any] ) -> [any] }

    @param a: first list
    @type  a: [any]
    @param b: second list
    @type  b: [any]
    @param optimize: result is sorted like the shorter of the two lists
                     (default: 1)
    @type  optimize: 1|0

    @return: list
    @rtype: [any]
    """
    if optimize and len(a) > len(b):
        a, b = b, a

    return [ x for x in a if x in b ]


def nonredundant( l ):
    """
    All members of a list without duplicates
    C{ noredundant( [any] ) -> [any] }

    @param l: list
    @type  l: [any]

    @return: list
    @rtype: [any]    
    """
    r = []
    for x in l:
        if x not in r:
            r.append( x )
    return r


def union( a, b ):
    """
    Union of two lists (without duplicates)
    C{ union( [any], [any] ) -> [any] }
    
    @param a: first list
    @type  a: [any]
    @param b: second list
    @type  b: [any]

    @return: list
    @rtype: [any]    
    """
    if type( a ) is N.arraytype:
        a = a.tolist()
    if type( b ) is N.arraytype:
        b = b.tolist()

    return nonredundant( a + b )


def difference( a, b ):
    """
    All members of a that are not in b
    C{ difference([any], [any]) -> [any] }

    @param a: first list
    @type  a: [any]
    @param b: second list
    @type  b: [any]

    @return: list
    @rtype: [any]    
    """
    return [ x for x in a if x not in b ]


def removeFromList( l, v, all=1 ):
    """
    Remove all or first occurrence(s) of v from l.
    C{ removeFromList( l, v, [all=1] ) }
    
    @param l: list
    @type  l: [ any ]
    @param v: remove these values
    @type  v: any or [ any ]
    @param all: remove all occurrences (1) or only first one (0) (default: 1)
    @type  all: 0|1

    @return: list
    @rtype: [any]     
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
    Creates a set of n unique integers randomly distributed between
    start and stop. 
    
    @param start: minimal index
    @type  start: int
    @param stop: 1+maximal index
    @type  stop: int
    @param n: number of indices
    @type  n: int
    
    @return: set of unique integers evenly distributed between start and stop
    @rtype: [int]
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
    see U{http://mathworld.wolfram.com/LeastSquaresFitting.html}
    
    @param x: x-data
    @type  x: [ float ]
    @param y: y-data
    @type  y: [ float ]
    
    @return: float, float - m, n, r^2 (slope, intersection, corr. coefficient)
    @rtype: float, float, float
    
    @raise BiskitError: if x and y have different number of elements
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
    Convert cartesian coordinate array to polar coordinate array: 
    C{ x,y,z -> r, S{theta}, S{phi} }
    
    @param xyz: array of cartesian coordinates (x, y, z)
    @type  xyz: array
    
    @return: array of polar coordinates (r, theta, phi)
    @rtype: array
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
    Convert polar coordinate array to cartesian coordinate array: 
    C{ r, S{theta}, S{phi} -> x,y,z }
     
    @param rtp: array of cartesian coordinates (r, theta, phi)
    @type  rtp: array
    
    @return: array of cartesian coordinates (x, y, z)
    @rtype: array
    """
    x = rtp[:,0] * N.cos( rtp[:,1] ) * N.sin( rtp[:,2] )
    y = rtp[:,0] * N.sin( rtp[:,1] ) * N.sin( rtp[:,2] )
    z = rtp[:,0] * N.cos( rtp[:,2] )

    return N.transpose( N.concatenate( ([x],[y],[z]) ) )


def projectOnSphere( xyz, radius=None, center=None ):
    """
    Project the coordinates xyz on a sphere with a given radius around
    a given center.
    
    @param xyz: cartesian coordinates
    @type  xyz: array N x 3 of float
    @param radius: radius of target sphere, if not provided the maximal
                   distance to center will be used (default: None)
    @type  radius: float
    @param center: center of the sphere, if not given the average of xyz
                   will be assigned to the center (default: None)
    @type  center: array 0 x 3 of float

    @return: array of cartesian coordinates (x, y, z)
    @rtype: array    
    """
    if center is None:
        center = N.average( xyz )

    if radius is None:
        radius = max( N.sqrt( N.sum( N.power( xyz - center, 2 ), 1 ) ) )

    rtp = cartesianToPolar( xyz - center )
    rtp[ :, 0 ] = radius

    return polarToCartesian( rtp ) + center
