## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Wolfgang Rieping
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

## last $Author$
## last $Date$
## $Revision$
"""
create a histogram from data
"""

import numpy as N


def histogram(data, nbins, range = None):
    """
    Create a histogram.
    Comes from Konrad Hinsen: Scientific Python

    @param data: data list or array
    @type  data: [any]
    @param nbins: number of bins
    @type  nbins: int
    @param range: data range to create histogram from (min val, max val)
    @type  range: (float, float) OR None

    @return: array (2 x len(data) ) with start of bin and witdh of bin. 
    @rtype: array
    """
    data = N.array(data, N.float)
    if range is None:
        min = data.min()
        max = data.max()
    else:
        min, max = range
        data = N.repeat(data,
                        N.logical_and(N.less_equal(data, max),
                                                  N.greater_equal(data, min)))
    bin_width = (max-min)/nbins
    data = N.floor((data - min)/bin_width).astype(N.int)
    histo = N.add.reduce(N.equal(
        N.arange(nbins)[:,N.newaxis], data), -1)
    histo[-1] = histo[-1] + N.add.reduce(N.equal(nbins, data))
    bins = min + bin_width*(N.arange(nbins)+0.5)
    return N.transpose(N.array([bins, histo]))


def density(x, nBins, range = None, steps = 1, hist = 0):
    """
    returns the normalized histogram of x::
      density( data, nBins [,range=None, steps=1, hist=0|1] ) -> array

    @param x: data list or array
    @type  x: [any]
    @param nBins: number of bins
    @type  nBins: int
    @param range: data range to create histogram from (min val, max val)
    @type  range: (float, float) OR None
    @param steps: 1: histogram appears as a discrete graph (default 1)
    @type  steps: 1|0
    @param hist: 0: normalize histogram (default 0)
    @type  hist: 1|0

    @return: array (2 x len(data) ) with start of bin and witdh of bin. 
    @rtype: array
    """
    h = histogram(x, nBins, range)
    binWidth = h[1,0] - h[0,0]

    if not hist:
        i = N.sum(h)[1]*binWidth
        h[:,1] = h[:,1]/i

    if steps:
        half = (h[1][0]-h[0][0])/2
        l = [(h[0][0]-half,0)]

        for row in h:
            l.append((row[0]-half,row[1]))
            l.append((row[0]+half,row[1]))

        l.append((h[-1][0]+half,0))

        h = l

    return N.array(h)


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def test_hist( self ):
        """hist test"""
        self.x = N.arange( 4, 12, 1.2 )
        self.data = density( self.x, 3, hist=1 )

        self.assert_( N.all( self.data == self.EXPECT) )

    EXPECT= N.array([[  4. ,   0. ],
                           [  4. ,   2. ],
                           [  6.4,   2. ],
                           [  6.4,   2. ],
                           [  8.8,   2. ],
                           [  8.8,   3. ],
                           [ 11.2,   3. ],
                           [ 11.2,   0. ]])


if __name__ == '__main__':

    BT.localTest()
