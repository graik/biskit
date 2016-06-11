#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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

## 
## 

from Numeric import *
import sys
import os, re
import os.path
from string import *
import biggles

from Biskit.tools import *
from Biskit import mathUtils
from Biskit import Trajectory
from Biskit import PDBDope
from Biskit import EHandler
from Biskit.Dock import Complex as ProteinComplex

def _use():
    print """
a_multiDock  Visualize multidock results
Syntax       a_multiDock -cl |complexList.cl|

             cl - complexList, has to contain info dictionary data for key

             inv - 1||0 inverse data associated with key (i.e. for rmds plots)
             maxContour - scale contour circles to fit at most x solutions

Default values:
    """
    default = _defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def _defOptions():
    return { 'cl':'complexes_cont.cl',
             'ref':'../../com_wet/ref.complex',
             'key':'fnac_10',
             'inv':0}


def nameFromPath( path ):
    """
    Extract a nice name from the absolute complexes_grouped.cg path
    """
    path = string.split( path, '/' )
    name = ''
    for p in path:
        if re.match( '^[a,c][0-9]{2}', p ):
            name += re.findall( '^[a,c][0-9]{2}', p )[0] + '-'
        if re.match(  '^dock_[0-9,a-z,_]{4}', p ):
            try:
                name += re.findall( '^dock_[0-9,a-z,_]{4}.+', p )[0]
            except:
                name += re.findall( '^dock_[0-9,a-z,_]{4}', p )[0]
            
    return name


def data3DList( c, key='fnac_10', inverse=0,
                rm=range(1,12), soln=512 ):
    """
    Create an matrix: len(rec_model) * len((lig_model) * solutions
    containing the values of the info dic with given key.
    c - ComplexList
    """
    rm = range( 1, max(c.valuesOf( 'model1'))+1 )
    lm = range( 1, max(c.valuesOf( 'model2'))+1 )

    matrix = zeros( ( len( rm ), len( lm ), soln ), 'f' )
    
    try:

        for r in rm:
            rl = c.filter( 'model1', r )
            for l in lm:
                cl = rl.filter( 'model2', l )
                if inverse:
                    matrix[r-1][l-1] = (1./array(
                        cl.valuesOf( key, default=0 ))).tolist()
                else:
                    matrix[r-1][l-1] = cl.valuesOf( key, default=0 )

    except ValueError, why:
        try:
            lenM = len( matrix[r-1][l-1] )
            lenV = len( cl.valuesOf( key, default=0 ) )
        except:
            lenM = lenV = 0
        s =  '%i : %i len(matrix)=%i <> len(values)=%i' % (r, l, lenM, lenV)
        EHandler.error('Cannot extract fnac data. '+ s )
    
    return matrix


def leastNumberOfZeros( data ):
    """
    get the least number of zeros in any data array
    -> int
    """
    ## data size
    rec, lig, soln = shape( data )

    ## the lowest number of zero fields in an array
    minZeros = soln

    for r in range(rec):
        for l in range(lig):
            if sum( array(data[r][l]) == 0. ) < minZeros:
                minZeros = sum( array(data[r][l]) == 0. )
    return minZeros


def argmax_2d( a ):
    """
    a - array( X x Y ) of numbers,
    -> position( x,y ) of highest item
    """
    r = argmax( ravel( a ) )
    len_x, len_y = shape( a )

    return r / len_y, r % len_y
    
    
        
def plotLines( nr, style='dot', col='black', ref=None ):
    """
    Greate a plot grid
    -> list of biggles objects
    """
    l = []

    for i in range( 1, nr+1 ):
        l += [ biggles.LineY( i, type = style, color = col ) ]
        l += [ biggles.LineX( i, type = style, color = col ) ]

    if ref:
        l += [ biggles.LineY( ref[0], color = 'grey' ) ]
        l += [ biggles.LineX( ref[1], color = 'grey' ) ]
   
    return l


def moreDockingInfo( data, levels=[0.01, .1, .2, .3, .4, .5, .6, .7, .8, .9]):
    r, l, s = shape( data )
    result = zeros( (r, l, len(levels)) )
    for i in range(len(levels)):
        for r in range( shape(data)[0] ):
            for l in range( shape(data)[1] ):
                result[r][l][i] =len( nonzero( greater(data[r][l], levels[i])))

    return result


def contourPlot( data, inverse=0, symbol='filled circle', label=None,
                 maxContour=None ):
    plots = []
    
    ## get grid lines
    plots = plotLines( len(data) )

    ## get data points
    cutoffs = [ 0.01, .1, .2, .3, .4, .5, .6, .7, .8, .9 ]
    if inverse:
        cutoffs = 1./array([ 15, 10, 8, 6, 5, 4, 3, 2, 1 ])
    data = moreDockingInfo( data, levels=cutoffs )

    ## callculate symbol size
    if maxContour:
        max_data = maxContour*1.
    else:
        max_data = max( ravel(data) )*1.

    max_symbol_size = 15
    scale_factor = max_data/(max_symbol_size**2)

    def scaleCircle( value ):
        return int(round(sqrt(value/scale_factor)))
    
    ## colors
    col = colorSpectrum( len(cutoffs))
    col.reverse()
    col=['light grey'] + col[1:]
    label_col = ['dark grey'] + col[1:]

    range_max = []
    range_min = []

    for i in range(len(cutoffs)):
        range_max += [ max( ravel( data[:,:,i] ) ) ]
        range_min += [ min( ravel( data[:,:,i] ) ) ]

        for r in range( shape(data)[0] ):

            for l in range( shape(data)[1] ):
        
                plots += [ biggles.Point( r+1, l+1,
                                          color=col[i],
                                          size=scaleCircle(data[r][l][i]),
                                          type=symbol ) ]
        
    ## legend
    plots += [ biggles.PlotLabel( 1.04, .98, 'FNC>', size=2 ) ]         
    for i in range( len(cutoffs) ):
        if sum( take(data, (i,), 2) )>0:
##             plots += [ biggles.PlotLabel( 1.03, .95-i*0.03, str(i*.1),
##                                           color=label_col[i],size=2 ) ]
            plots += [ biggles.PlotLabel( 1.03, .95-i*0.03, '%.2f'%cutoffs[i],
                                          color=label_col[i],size=2 ) ]
            plots += [ biggles.PlotLabel( 1.03, .35-i*0.02,
                                          str(range_min[i])+'-'+str(range_max[i]),
                                          color=label_col[i],size=1 ) ]

    return plots


def getRmsdValues( rmsd_dic, mask ):
    
    rb = rmsd_dic[ mask ]['rec'][0]
    rf = rmsd_dic[ mask ]['rec'][1]
    r = [ ['%.2f'%rb[i], '%.2f'% rf[i] ] for i in range( len(rb) )]
    
    lb = rmsd_dic[ mask ]['lig'][0]
    lf = rmsd_dic[ mask ]['lig'][1]
    l = [ ['%.2f'%lb[i], '%.2f'% lf[i] ] for i in range( len(lb) )]

    m = [ str(i) for i in range( len(lb)) ]
    
    return m[1:], r[1:], l[1:]


def plotRmsLabels( nr, rmsd_dic, mask, pad=1 ):
    """
    nr - number of models (usually 1+10)
    mask - names of maska to use
    
    add rmsd values to plot
    -> biggles object
    """
    ## label position tweeking
    step = 1. / ( (nr -1) + (pad*2) )
    shift = step*pad + array( range(0, nr), 'f' )*step
    
    p = []
    for c in range( 0, nr ):
        m, r, l = getRmsdValues( rmsd_dic, mask )
        
        p += [ biggles.PlotLabel( shift[c], -0.02 , r[c][0], size=1 ) ]
        p += [ biggles.PlotLabel( shift[c], -0.04, r[c][1], size=1 ) ]
        p += [ biggles.PlotLabel( shift[c], -0.09,  m[c],    size=3 ) ]

        p += [ biggles.PlotLabel( -0.03 , shift[c]+0.01, l[c][0], size=1 ) ]
        p += [ biggles.PlotLabel( -0.03, shift[c]-0.01, l[c][1], size=1 ) ]
        p += [ biggles.PlotLabel( -0.09,  shift[c], m[c],    size=3 ) ]
      
    return p



###############
# main function

if len(sys.argv) < 2:
    _use()
    
options = cmdDict( _defOptions() )
#    options = test()

## paths ...
cFile = os.path.abspath( options['cl'] )
cName = nameFromPath( cFile ) + '-' + options['key']
dock_dir = os.path.split( absfile(options['cl']))[0]

## file names for output ...
contourPlotName = os.path.dirname(cFile) + '/' + cName + '_contour_plot.eps'


## Load data
flushPrint('Loading complex list \n')
cList = load( options['cl'] )

## extract data
flushPrint('Extracting %s data  \n'%options['key'])
inverse_data = int( options['inv'])

data = data3DList( cList, key=options['key'], inverse=inverse_data,
                   soln=512 )

#########
## Contour plot 
colorDic = {'x':'red', 'r':'cornflower blue', 'l':'pale green', 'a':'grey54'}
biggles.configure('fontsize_min', 0.7 )

p_cont = biggles.FramedPlot()
p_cont.x1.draw_ticklabels = 0
p_cont.y1.draw_ticklabels = 0
p_cont.title = cName

## labels
p_cont.x2.label_style['size'] = 2
p_cont.y2.label_style['size'] = 2
p_cont.x2.label = 'Receptor model'
p_cont.y2.label = 'Ligand model'
p_cont.x1.label = '  '
p_cont.y1.label = '  '

models = len( data )

p_cont.xrange = ( 0, models+1 )
p_cont.yrange = ( 0, models+1 )

plots = contourPlot( data, inverse=inverse_data, symbol='filled circle',
                     maxContour=toInt( options.get('maxContour',None)) )
#plots += plotRmsLabels( 11, rmsd_all_dic, 'CONT')
#plots += plotNrLabels( 11, pos=-0.05 ) 

for p in plots:
    p_cont.add( p )

p_cont.show()
p_cont.write_eps( contourPlotName )
