##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Wolfgang Rieping
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

## last $Author$
## last $Date$
## $Revision$
"""
create a histogram from data
"""

import Numeric


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
    data = Numeric.array(data, Numeric.Float)
    if range is None:
	min = Numeric.minimum.reduce(data)
	max = Numeric.maximum.reduce(data)
    else:
	min, max = range
	data = Numeric.repeat(data,
		   Numeric.logical_and(Numeric.less_equal(data, max),
                                       Numeric.greater_equal(data, min)))
    bin_width = (max-min)/nbins
    data = Numeric.floor((data - min)/bin_width).astype(Numeric.Int)
    histo = Numeric.add.reduce(Numeric.equal(
	Numeric.arange(nbins)[:,Numeric.NewAxis], data), -1)
    histo[-1] = histo[-1] + Numeric.add.reduce(Numeric.equal(nbins, data))
    bins = min + bin_width*(Numeric.arange(nbins)+0.5)
    return Numeric.transpose(Numeric.array([bins, histo]))


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
        i = Numeric.sum(h)[1]*binWidth
        h[:,1] = h[:,1]/i

    if steps:
        half = (h[1][0]-h[0][0])/2
        l = [(h[0][0]-half,0)]

        for row in h:
            l.append((row[0]-half,row[1]))
            l.append((row[0]+half,row[1]))

        l.append((h[-1][0]+half,0))

        h = l

    return Numeric.array(h)


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_hist( self ):
        """hist test"""
        self.x = Numeric.arange( 4, 12, 1.2 )
        self.data = density( self.x, 3, hist=1 )

        self.assert_( self.data, self.EXPECT )

    EXPECT= Numeric.array([[  4. ,   0. ],
			   [  4. ,   2. ],
			   [  6.4,   2. ],
			   [  6.4,   2. ],
			   [  8.8,   2. ],
			   [  8.8,   3. ],
			   [ 11.2,   3. ],
			   [ 11.2,   0. ]])


if __name__ == '__main__':

    BT.localTest()
