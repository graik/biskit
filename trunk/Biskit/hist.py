##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Wolfgang Rieping; All rights reserved
##
## create a histigram from data
## 
## Author: Wolfgang Rieping
##
## last $Author$
## last $Date$
## $Revision$


import Numeric


def histogram(data, nbins, range = None):
    """
    Comes from Konrad Hinsen: Scientific Python
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
    density( data, nBins [,range=None, steps=1, hist=0|1] ) -> array
    returns the normalized histogram of x
    steps = 1: histogram appears as a discrete graph
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



