## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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

import numpy as N
import random
import math, cmath

class MathUtilError( Exception ):
    pass

def accumulate( a ):
    """
    cumulative sum of C{ a[0], a[0]+a[1], a[0]+a[1]+[a2], ... }
    normalized by C{ N.sum( a ) }

    @param a: array('f') or float
    @type  a: array

    @return: float
    @rtype: float
    """
    return N.cumsum( a ) / N.sum( a )


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


def aboveDiagonal( pw_m ):
    """
    Collect all the values above the diagonal in a square
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
    if mask is not None:
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
    if mask is not None:
        r = N.zeros(a*b)
        N.put( r, N.nonzero(mask), ranArray)
        return N.reshape( r, (a,b) )

    if not mask:
        return  N.reshape( ranArray, (a,b) )


def slidingAverage( y, window=2 ):
    if window == 0:
        return y

    assert window < len(y), 'window size too large for array'

    margin = int(round((window-1)/2.))

    return [ N.average( y[i-margin : i+margin] )
             for i in range(margin,len(y)-margin) ]


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

def area(curve, start=0.0, stop=1.0 ):
    """
    Numerically add up the area under the given curve.
    The curve is a 2-D array or list of tupples.
    The x-axis is the first column of this array (curve[:,0]).
    (originally taken from Biskit.Statistics.ROCalyzer)

    @param curve: a list of x,y coordinates
    @type  curve: [ (y,x), ] or N.array
    @param start: lower boundary (in x) (default: 0.0)
    @type  start: float
    @param stop: upper boundary (in x) (default: 1.0)
    @type  stop: float
    @return: the area underneath the curve between start and stop.
    @rtype: float
    """
    ## convert and swap axes
    curve = N.array( curve )
    c = N.zeros( N.shape(curve), curve.dtype )
    c[:,0] = curve[:,1]
    c[:,1] = curve[:,0]

    assert len( N.shape( c ) ) == 2

    ## apply boundaries  ## here we have a problem with flat curves
    mask = N.greater_equal( c[:,1], start )
    mask *= N.less_equal( c[:,1], stop )
    c = N.compress( mask, c, axis=0 )

    ## fill to boundaries -- not absolutely accurate: we actually should
    ## interpolate to the neighboring points instead
    c = N.concatenate((N.array([[c[0,0], start],]), c,
                       N.array([[c[-1,0],stop ],])) )
    x = c[:,1]
    y = c[:,0]

    dx = x[1:] - x[:-1] # distance on x between points 
    dy = y[1:] - y[:-1] # distance on y between points

    areas1 = y[:-1] * dx  # the rectangles between all points
    areas2 = dx * dy / 2.0 # the triangles between all points

    return N.sum(areas1) + N.sum(areas2)



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

    m = N.zeros( N.cumproduct( s )[-1], N.int)
    pass  ## m.savespace( 1 )
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

    R = N.zeros((3,3), N.float32)

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
    alpha = N.random.random() * 2 * N.pi
    gamma = N.random.random() * 2 * N.pi
    beta  = N.arccos(2*(N.random.random() - 0.5))

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
    if type( a ) is N.ndarray:
        a = a.tolist()
    if type( b ) is N.ndarray:
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

    @return: m, n, r^2 (slope, intersection, corr. coefficient)
    @rtype: float, float, float

    @raise BiskitError: if x and y have different number of elements
    """
    x, y = N.array( x, N.float64), N.array( y, N.float64)
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


def rotateAxis(theta, vector):
    """
    Calculate a left multiplying rotation matrix that rotates
    theta rad around vector.

    Taken from: http://osdir.com/ml/python.bio.devel/2008-09/msg00084.html

    Example:

    >>> m=rotateAxis(pi, N.array([1,0,0]))

    @type theta: float
    @param theta: the rotation angle


    @type vector: L{Vector}
    @param vector: the rotation axis

    @return: The rotation matrix, a 3x3 Numeric array.
    """
    vector = vector / N.linalg.norm(vector)
    x,y,z=vector
    c=N.cos(theta)
    s=N.sin(theta)
    t=1-c
    rot=N.zeros((3,3), "d")
    # 1st row
    rot[0,0]=t*x*x+c
    rot[0,1]=t*x*y-s*z
    rot[0,2]=t*x*z+s*y
    # 2nd row
    rot[1,0]=t*x*y+s*z
    rot[1,1]=t*y*y+c
    rot[1,2]=t*y*z-s*x
    # 3rd row
    rot[2,0]=t*x*z-s*y
    rot[2,1]=t*y*z+s*x
    rot[2,2]=t*z*z+c
    return rot


def cbrt(x):
    """
    cubic root.
    Author: Victor Gil.
    """
    if x >= 0: 
        return math.pow(x, 1.0/3.0) 
    else:
        return -math.pow(N.abs(x), 1.0/3.0) 


def cartesian2D(r, w, deg=0): # radian if deg=0; degree if deg=1 
    """ 
    Convert from polar (r,w) to rectangular (x,y) x = r cos(w) y = r sin(w)
    Author: Victor Gil
    """ 
    from math import cos, sin, pi 
    if deg: 
        w = pi * w / 180.0 
        return r * cos(w), r * sin(w) 

def polar2D(x, y, deg=0): # radian if deg=0; degree if deg=1 
    """
    Convert from rectangular (x,y) to polar (r,w) r = sqrt(x^2 + y^2) 
    w = arctan(y/x) = [-\pi,\pi] = [-180,180] 

    Author: Victor Gil
    """ 
    from math import hypot, atan2, pi 
    if deg: 
        return hypot(x, y), 180.0 * atan2(y, x) / pi 
    else: 
        return hypot(x, y), atan2(y, x) 


def quadratic(a, b, c=None): 
    """ 
    x^2 + ax + b = 0 (or ax^2 + bx + c = 0) By substituting x = y-t and t = a/2,
    the equation reduces to y^2 + (b-t^2) = 0 which has easy solution 
    y = +/- sqrt(t^2-b)

    Author: Victor Gil
    """ 
    if c: # (ax^2 + bx + c = 0) 
        a, b = b / float(a), c / float(a)
    t = a / 2.0 
    r = t**2 - b 
    if r >= 0: # real roots 
        y1 = math.sqrt(r) 
    else: # complex roots 
        y1 = cmath.sqrt(r) 
    y2 = -y1 
    return y1 - t, y2 - t 

def cubic(a, b, c, d=None):
    """
    x^3 + ax^2 + bx + c = 0  (or ax^3 + bx^2 + cx + d = 0)
    With substitution x = y-t and t = a/3, the cubic equation reduces to    
        y^3 + py + q = 0,
    where p = b-3t^2 and q = c-bt+2t^3.  Then, one real root y1 = u+v can
    be determined by solving 
        w^2 + qw - (p/3)^3 = 0
    where w = u^3, v^3.  From Vieta's theorem,
        y1 + y2 + y3 = 0
        y1 y2 + y1 y3 + y2 y3 = p
        y1 y2 y3 = -q,
    the other two (real or complex) roots can be obtained by solving
        y^2 + (y1)y + (p+y1^2) = 0

    Author: Victor Gil
    """
    cos = math.cos
    if d:       # (ax^3 + bx^2 + cx + d = 0)
        a, b, c = b / float(a), c / float(a), d / float(a)
    t = a / 3.0
    p, q = b - 3 * t**2, c - b * t + 2 * t**3
    u, v = quadratic(q, -(p/3.0)**3)
    if type(u) == type(0j):   # complex cubic root
        r, w = polar2D(u.real, u.imag)
        y1 = 2 * cbrt(r) * cos(w / 3.0)
    else:     # real root
        y1 = cbrt(u) + cbrt(v)
    y2, y3 = quadratic(y1, p + y1**2)
    return y1 - t, y2 - t, y3 - t

def outliers( a, z=5, it=5 ):
    """
    Iterative detection of outliers in a set of numeric values.
    Requirement: len(a) > 0; outlier detection is only performed if len(a)>2
    
    @param a: array or list of values
    @type  a: [ float ]
    @param z: z-score threshold for iterative refinement of median and SD
    @type  z: float
    @param it: maximum number of iterations
    @type  it: int
    
    @return: outlier mask, median and standard deviation of last iteration
    @rtype: N.array( int ), float, float
    """
    assert( len(a) > 0 )
    mask = N.ones( len(a) )
    out  = N.zeros( len(a) )
    
    if len(a) < 3:
        return out, N.median(a), N.std(a)
    
    for i in range( it ):
        b  = N.compress( N.logical_not(out), a )
        me = N.median( b )
        sd = N.std( b )
        
        bz = N.absolute((N.array( a ) - me) / sd)  # pseudo z-score of each value
        o  = bz > z
        ##        print 'iteration %i: <%5.2f> +- %5.2f -- %i outliers' % (i,me,sd,N.sum(o))

        ## stop if converged or reached bottom
        if (N.sum(o) == N.sum(out)) or (N.sum(o) > len(a) - 3):
            return o, me, sd
            
        out = o
    
    return out, me, sd
        
    
#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_mathUtils(self):
        """mathUtils.polar/euler test"""
        ## Calculating something ..
        self.d = N.array([[20.,30.,40.],[23., 31., 50.]])

        self.a = polarToCartesian( cartesianToPolar( self.d ) )

        self.t = eulerRotation( self.a[0][0], self.a[0][1], self.a[0][2]  )

        self.assertAlmostEqual( N.sum( SD(self.a) ), self.EXPECT )

    def test_area(self):
        """mathUtils.area test"""
        self.c = zip( N.arange(0,1.01,0.1), N.arange(0,1.01,0.1) )
        self.area = area( self.c )
        self.assertAlmostEqual( self.area, 0.5, 7 )

    EXPECT = N.sum( N.array([ 2.12132034,  0.70710678,  7.07106781]) )

if __name__ == '__main__':

    BT.localTest()
