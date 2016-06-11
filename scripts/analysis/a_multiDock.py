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
from Biskit.Dock import Complex as ProteinComplex

def _use():
    print """
a_multiDock  Visualize multidock results
Note: interface rms values are for contact atoms not contact residues.

Syntax       a_multiDock -cl |complexList.cl|

             cl - complexList, has to contain info dictionary data for key
             r  - hex receptor pdbs (i.e rec/*_hex.pdb)
             l  - hex ligand pdbs (i.e lig/*_hex.pdb)
             ref - reference complex
             key - info dictionary key to plot (high values are considered good)
             inv - 1||0 inverse data associated with key (i.e. for rmds plots)
             maxContour - scale contour circles to fit at most x solutions
             additional_profile - add to profile plot (rec_model lig_model)
             
Result       5 plots, info txt file, dumped data

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
                rm=range(1,12), lm=range(1,12), soln=512 ):
    """
    Create an matrix: len(rec_model) * len((lig_model) * solutions
    containing the values of the info dic with given key.
    c - ComplexList
    """
    rm = toList( rm )
    lm = toList( lm )

    matrix = zeros( ( len( rm ), len( lm ), soln ), 'f' )

    for r in rm:
        rl = c.filter( 'model1', r )
        for l in lm:
            cl = rl.filter( 'model2', l )
            if inverse:
                matrix[r-1][l-1] = (1./array(cl.valuesOf( key ))).tolist()
            else:
                matrix[r-1][l-1] = cl.valuesOf( key )
                  
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


def getArea( x, y ):
    x = array(x)
    y = array(y)
    
    x = abs( x[1:]-x[:-1] )
    y = (y[1:]+y[:-1])/2
    
    return sum( y * x )


def docking_performance( data, cutoff, refCom ):
    """
    Get a single number docking performance for all 121 combinations of rec
    and lig model.
    data - array(11x11x512) of float, fnrc_4.5acts for all dockings
           each 512 values must be sorted
    cutoff - int, throw away all solutions with less than cutoff contacts
    refComplex - Complex, real rec-lig complex (needed for number of native
                 contacts )
    -> array( 11 x 11 ) of float, docking 'performance' of ech rec - lig pair
    """
    ## convert number of contacts into fraction of nat. contacts
    refContacts = ravel( refCom.resContacts() )
    ## allow for small discrepancy of fnrc_4.5
    fract_cutoff = 1.0 /  sum( refContacts  )
    fract_cutoff = fract_cutoff * (cutoff - 0.5)

    ## clip values below cutoff
    data = data * logical_not( less( data, fract_cutoff ) )

    result = sum( cumsum( data, 2 ), 2 )

    return result


def argmax_2d( a ):
    """
    a - array( X x Y ) of numbers,
    -> position( x,y ) of highest item
    """
    r = argmax( ravel( a ) )
    len_x, len_y = shape( a )

    return r / len_y, r % len_y
    
    
        
def diffPlot( data, lenZeros, colorDic, type='fnc', ref=(1,1) ):
    """
    Plot of accumulated better solutions at each level
    data - 3D matrix, sorted data
    lenZeros - integer, the lowest number of zeros on any rec-lig pair
    ref - tuple, position of reference
    -> plots - list of plot objects
       areas - matrix, area between ref and data
    """
    ## reference data array
    refData = data[ ref[0]-1 ][ ref[1]-1 ]

    ## data size
    rec, lig, soln = shape( data )

    ## area matrix
    areas = zeros( ( rec, lig ), 'f' )

    x = range( 1, soln-(lenZeros-1) )
    x.reverse()

    plots = []

    def plotFnc( x, y, color=colorDic['a'] ):
        return [ biggles.Curve( x, cumsum(y), color=color ) ]
        #return [ biggles.Curve( x, y, color=color ) ]
    
    def plotBetter( x, y, color=colorDic['a'] ):
        y = cumsum( less(0, y) - greater(0, y) )
        return [ biggles.Curve( x, y, color=color ) ]
    
    for r in range(rec):
        for l in range(lig):
            diff = ( data[r][l] - refData )[-(soln-lenZeros):]
            if type=='better':
                plots += plotBetter( x, diff )
                ## area between ref and current
                score = cumsum( less(0, diff) - greater(0, diff) )
                areas[r][l] = getArea( x, score )
                
            if type=='fnc':
                plots += plotFnc( x, diff )
                ## area between ref and current
                areas[r][l] = getArea( x, cumsum( diff ))
                #areas[r][l] = getArea( x, diff )
    
    ## contains x-ray rec
    for l in range(lig):
        diff = ( data[0][l] - refData )[-(soln-lenZeros):]       
        if type=='better':
            plots += plotBetter( x, diff, color=colorDic['r'] )
        if type=='fnc':
            plots += plotFnc( x, diff, color=colorDic['r'] )
  
    ## contains x-ray lig
    for r in range(rec):
        diff = ( data[r][0] - refData)[-(soln-lenZeros):]
        if type=='better':
            plots += plotBetter( x, diff, color=colorDic['l'] )
        if type=='fnc':
            plots += plotFnc( x, diff, color=colorDic['l'])

    ## ref line
    plots += [ biggles.Curve( x, zeros(len(x)), color=colorDic['x'], width=3 ) ]
                
    return plots, areas


def logPlot( data, lenZeros, colorDic, ref=None, ave=None ):
    """
    log(fraction of solutions) vs. fraction of native contacts (sorted)
    data - 3D matrix (rec * lig * soln)
    lenZeros - integer, the lowest number of zeros on any rec-lig pair
    ref - tuple, position of reference
    -> list of plot objects
    """
    ## data dimensions
    rec, lig, soln = shape( data )

    ## y data range (removing all zeros region)
    y = arange( soln-(lenZeros), 0, -1, 'f' )/soln

    plots = []

    if type(colorDic)==str:
        colorDic = {'a':colorDic}

    ## collect all
    for r in range(rec):
        for l in range(lig):
            dt = data[r][l][-(soln-lenZeros):]
            if ave:
                dt = mathUtils.runningAverage( dt, ave, 1 )
            plots += [ biggles.Curve( dt , y, color=colorDic['a'] ) ]

    if ref:
        nr = ref[0]-1
        nl = ref[1]-1
        ## contains x-ray rec
        for l in range(lig):
            dt = data[nr][l][-(soln-lenZeros):]
            if ave:
                dt = mathUtils.runningAverage( dt, ave, 1 )
            plots += [ biggles.Curve( dt , y, color=colorDic['r'] ) ]

        ## contains x-ray lig
        for r in range(rec):
            dt = data[r][nl][-(soln-lenZeros):]
            if ave:
                dt = mathUtils.runningAverage( dt, ave, 1 )
            plots += [ biggles.Curve( dt , y, color=colorDic['l'] ) ]

        ## add ref in other color
        dt = data[nr][nl][-(soln-lenZeros):]
        if ave:
            dt = mathUtils.runningAverage( dt, ave, 1 )
        plots += [ biggles.Curve( dt , y, color=colorDic['x'], width=4 ) ]

        ## labels
        plots += [ biggles.PlotLabel( .8, .95, 'Xray rec', color =colorDic['r'] ) ]
        plots += [ biggles.PlotLabel( .8, .90, 'Xray lig', color =colorDic['l'] ) ]
        plots += [ biggles.PlotLabel( .8, .85, 'Xray', color =colorDic['x'] ) ]

    return plots


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


def categorize( mat, posRange=None ):
    if posRange:
        range = posRange*2/20
    else:
        max_area = max( ravel( mat ) )
        min_area = min( ravel( mat ) )
        range = ( max_area - min_area ) / 10
    
    neg = greater( 0 , mat )
    pos = less( 0 ,mat )

    cat_pos = mat/range * pos
    cat_neg = mat/range * neg

    cat_pos = cat_pos.astype('i') + pos
    cat_neg = cat_neg.astype('i') - neg

    return cat_pos, cat_neg


def plotMatrixLabels( matrix, pad=1.0, color='black', offset = 0.03 ):
    """
    memberMat - frame cluster membership matrix
    offset - label position tweeking
    -> biggles object
    """
    a = []

    ## data dimensions
    rec, lig  = shape( matrix )

    ## label positions
    step = 1. / ( (rec-1) + (pad*2) )
    shift = step*pad + array( range(0, rec), 'f' )*step 

    ## iterate over matrix
    for r in range( rec ):
        for l in range( lig ):
            d = matrix[r][l]
            if type(d)==int or type(d)== float:
                a += [ biggles.PlotLabel( shift[r]+offset, shift[l]+offset,
                                          str(round( d, 1 )),
                                          color = color, size = 1 ) ]
            if type(d)==str:
                if matrix[r][l] != '0':
                    a += [ biggles.PlotLabel( shift[r]+offset, shift[l]+offset,
                                              d, color = color, size = 1 ) ]
    return a


def visualizePlot( areaMat, pos_symbol='circle', neg_symbol='circle',
                   pos_color='red', neg_color='blue', label=None, posRange=None ):
    plots = []
    
    ## get grid lines
    plots = plotLines( models, ref=(1,1) )

    ## get color points
    pos, neg = categorize( areaMat, posRange )
    
    for r in range( shape(areaMat)[0] ):
        for l in range( shape(areaMat)[1] ):
            plots += [ biggles.Point( r+1, l+1, color=pos_color, size=pos[r][l], \
                                       type=pos_symbol ) ]
            plots += [ biggles.Point( r+1, l+1, color=neg_color, size=neg[r][l], \
                                       type=neg_symbol ) ]

    ## get labels
    if label:
        plots += plotMatrixLabels( areaMat, color='dark slate grey' )

    return plots


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



def plotNrLabels( nr, pos=-0.05, pad=1 ):
    """
    nr - number of models (usually 1+10)
    
    add xy tick labels
    -> biggles object
    """
    ## label position tweeking
    step = 1. / ( (nr -1) + (pad*2) )
    shift = step*pad + array( range(0, nr), 'f' )*step
    
    p = []
    for c in range( 0, nr ):        
        p += [ biggles.PlotLabel( shift[c], pos,  str(c), size=3 ) ]
        p += [ biggles.PlotLabel( pos,  shift[c], str(c), size=3 ) ]
      
    return p


def moreDockingInfo( data, levels=[0.01, .1, .2, .3, .4, .5, .6, .7, .8, .9 ] ):
    r, l, s = shape( data )
    result = zeros( (r, l, len(levels)) )
    for i in range(len(levels)):
        for r in range(11):
            for l in range(11):
                result[r][l][i] = len( nonzero( greater( data[r][l], levels[i] ) ) )

    return result


def dockingInfo( areaMat , fnc3DMat, ref=(1,1) ):
    x_rec, x_lig = ref
    
    ## best models
    ligSum = sum( areaMat )
    recSum = sum( areaMat, 1 )
    lm = (ligSum / min(abs(ligSum))).astype('i')
    rm = (recSum / min(abs(recSum))).astype('i')
    lm = nonzero( lm == max(lm) )
    rm = nonzero( rm == max(rm) )

    ## best model combination
    maxMat = max(ravel(areaMat))
    maxData = nonzero(( equal( maxMat, ravel(areaMat) )))
    r, l = int( ( maxData / len(areaMat) )[0] )+1,\
           int( ( maxData % len(areaMat))[0])+1
    max_fnc_soln = max( fnc3DMat[r-1][l-1] )
    soln = nonzero( equal( max_fnc_soln, fnc3DMat[r-1][l-1] ))[0]
    
    ## overall best soln
    max_overall = max( ravel(fnc3DMat) )
    max_soln = nonzero( equal( max_overall, ravel(fnc3DMat) ))[0]
    mIdx =  ( max_soln / 512 )
    rec_max_soln, lig_max_soln  = ( mIdx / len(areaMat) ) +1,\
                                  ( mIdx % len(areaMat) ) +1
    max_soln = max_soln % 512
    
    max_fnc_ref = max( fnc3DMat[x_rec-1][x_lig-1] ) 
    
    msg = """
The best performing model(s) are:
    receptor: %s
    ligand:   %s

The best receptor/ligand model combination is %i / %i
    where solution %i is  best solution with a
    fraction of native contacts of %.3f

The single best solution is given receptor/ligand model
    combination %i / %i where solution %i has a
    fnrc_4.5 of %.3f (max in reference %.3f)""" \
    %(rm+1, lm+1,
      r, l, soln+1, max_fnc_soln, 
      rec_max_soln, lig_max_soln, max_soln+1, max_overall, max_fnc_ref)
    
    return msg, (r,l), (rec_max_soln, lig_max_soln)


def maxVal( data ):
    """
    go through array and keep the highest value so far
    """
    m = [ data[0] ]
    high = data[0]
    for i in range( 1, len(data) ):
        if high < data[i]:
            high = data[i]
            m += [ high ]
        else:
            m += [ high ]
    return m

    
def plotRanking( data, colorDic, ref=(1,1) ):
    """
    log(fraction of solutions) vs. fraction of native contacts (sorted)
    data - 3D matrix (rec * lig * soln)
    ref - tuple, position of reference
    -> list of plot objects
    """
    ## data dimensions
    rec, lig, soln = shape( data )

    ## y data range (removing all zeros region)
    x = arange( 1, soln+1, 1, 'f' )/soln
    area = zeros( ( rec, lig ), 'f' )
    plots = []

    ## collect all
    ref_area = getArea( x, (maxVal( data[ref[0]-1][ref[1]-1] )) )
    for r in range(rec):
        for l in range(lig):
            d = maxVal(data[r][l])
            area[r][l] = getArea( x, d ) - ref_area
            plots += [ biggles.Curve( x, d, color=colorDic['a'] )]

    if ref:
        nr = ref[0]-1
        nl = ref[1]-1
        ## contains x-ray rec
        for l in range(lig):
            plots += [ biggles.Curve( x, maxVal(data[nr][l]),
                                      color=colorDic['r'])]

        ## contains x-ray lig
        for r in range(rec):
            plots += [ biggles.Curve( x, maxVal(data[r][nl]),
                                      color=colorDic['l']) ]

        ## add ref in other color
        plots += [ biggles.Curve( x, maxVal(data[nr][nl]),
                                  color=colorDic['x'])]

        ## labels
        plots += [ biggles.PlotLabel( .8, .95, 'Xray rec', color =colorDic['r'] ) ]
        plots += [ biggles.PlotLabel( .8, .90, 'Xray lig', color =colorDic['l'] ) ]
        plots += [ biggles.PlotLabel( .8, .85, 'Xray', color =colorDic['x'] ) ]

    return plots, area

    
def analyzeRmsd( t, aMask=None ):
    """
    traj - Lst of Trajectories with equal atoms
    aMask - atomMask to apply
    -> 1xN arrays of pairwise rmsds
    """
    ## create copies with less atoms, if requested
    if aMask is not None:
        traj = t.compressAtoms( aMask )


    
    ## calculate pw rmsd
    pw = traj.pairwiseRmsd()

    return pw


def calcAllRmsds( rec, lig, masks ):
    """
    Calculate pairwise rmsds for all combinations of free and
    bound of the sets specified in masks .
    """
    result = {}
    
    for m in masks:
        mask_name = m['name']
        ## pairwise rmsd calculation
        flushPrint("analyzing rms ...\n")
        rec_pw  = analyzeRmsd( rec, aMask= m['rec'])
        lig_pw  = analyzeRmsd( lig, aMask= m['lig'])
        
        result.update( {mask_name:{'rec':rec_pw,
                                   'lig':lig_pw,}})
    return result


def allMasks( traj_rec, traj_lig, com ):
    """
    Return a list of masks 
    The number of atoms in the first mask is used to generate
           the % comment in the rest
    """
    masks =[]
    ## number of atoms in rec and lig (after casting)
    nRec = traj_rec.getRef().lenAtoms()
    nLig = traj_lig.getRef().lenAtoms()

    ## ALL - all heavy atom mask
    masks+= [{'rec':ones( nRec ),
              'lig':ones( nLig ),
              'name':'ALL'}]
    
    ## CONT - contact atom mask    
    aContMat = com.atomContacts() #heavy atom contacts
    rec_cont = greater( sum( aContMat, 1 ), 0 )
    lig_cont = greater( sum( aContMat, 0 ), 0 )
    masks+= [{'rec':rec_cont,
              'lig':lig_cont,
              'name':'CONT'}]
    
    ## NONE_CONT - none contact surface atom mask
    masks+= [ {'rec':logical_not( rec_cont ),
               'lig':logical_not( lig_cont ),
               'name':'NONE_CONT'}]

    ## BACKBONE - Back bone mask
    rec_bb = traj_rec.getRef().maskBB()
    lig_bb = traj_lig.getRef().maskBB()
    masks+= [ {'rec':rec_bb,
               'lig':lig_bb,
               'name':'BACKBONE'}]
    
    ## SIDECHAIN - Side chain mask
    rec_sc = logical_not(rec_bb)
    lig_sc = logical_not(lig_bb)
    masks+= [{'rec':rec_sc,
              'lig':lig_sc,
              'name':'SIDECHAIN'}]
    
    ## CARBON - hydrophobic - Element: C
    C_rec = traj_rec.getRef().maskF( lambda a: a['element'] == 'C')
    C_lig = traj_lig.getRef().maskF( lambda a: a['element'] == 'C')
    masks+= [{'rec':C_rec,
              'lig':C_lig,
              'name':'CARBON'}]
    
    ## POLAR - polar - Element: O, N, S
    elements = ['O', 'N', 'S']
    polar_rec = array( zeros(nRec) )
    polar_lig = array( zeros(nLig) )
    for e in elements:
        polar_rec += traj_rec.getRef().maskF( lambda a: a['element'] == e )
        polar_lig += traj_lig.getRef().maskF( lambda a: a['element'] == e )
    
    masks+= [{'rec':polar_rec,
              'lig':polar_lig,
              'name':'POLAR'}]
    
    ## AA_TYPE: CHARGED, POLAR or NONPOLAR
    charged = [ 'LYS', 'ARG', 'HIS','ASP', 'GLU' ]
    polar = [ 'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS' ]
    nonpolar = [ 'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP' ]

    aa_charged_rec =  array( zeros(nRec) )
    aa_polar_rec =  array( zeros(nRec) )
    aa_nonpolar_rec =  array( zeros(nRec) )
    aa_charged_lig = array( zeros(nLig) )
    aa_polar_lig = array( zeros(nLig) )
    aa_nonpolar_lig = array( zeros(nLig) )
    
    for c in charged:
        aa_charged_rec += traj_rec.getRef().maskF( lambda a: a['residue_name'] == c )
        aa_charged_lig += traj_lig.getRef().maskF( lambda a: a['residue_name'] == c )

    for p in polar:
        aa_polar_rec += traj_rec.getRef().maskF( lambda a: a['residue_name'] == p )
        aa_polar_lig += traj_lig.getRef().maskF( lambda a: a['residue_name'] == p )

    for n in nonpolar:
        aa_nonpolar_rec += traj_rec.getRef().maskF( lambda a: a['residue_name'] == n )
        aa_nonpolar_lig += traj_lig.getRef().maskF( lambda a: a['residue_name'] == n )

    masks+= [{'rec':aa_charged_rec,
              'lig':aa_charged_lig,
              'name':'AA_CHARGED'}]
    masks+= [{'rec':aa_polar_rec,
              'lig':aa_polar_lig,
              'name':'AA_POLAR'}]
    masks+= [{'rec':aa_nonpolar_rec,
              'lig':aa_nonpolar_lig,
              'name':'AA_NONPOLAR'}]

    ## add statistics
    masks[0]['rec_comment'] = nRec
    masks[0]['lig_comment'] = nLig
    for m in masks[1:]:
        m['rec_comment'] = 100.*sum(m['rec'])/nRec
        m['lig_comment'] = 100.*sum(m['lig'])/nLig

    mask_names = [ m['name'] for m in masks ]

    return masks, mask_names


def surfMasks( traj_rec, traj_lig, com ):
    """
    surface atom filtered masks
    """
    ## surface masks (atoms vith more that 0.01% exposed ASA)
    rec_surf = greater( traj_rec.getRef().atomProfile('relASA'), 0.01 )
    lig_surf = greater( traj_lig.getRef().atomProfile('relASA'), 0.01 )
    
    ## number of surface atoms in rec, lig
    nRec = sum( rec_surf )
    nLig = sum( lig_surf )

    ## collect masks
    surf_masks = []
    surf_mask_names = []

    ## get all atom masks (once more)
    masks, mask_names = allMasks( traj_rec, traj_lig, com )

    name = ''
    try:
    
        ## correct all but first mask 
        for i in range(len(masks)):
            name = 'SURF_' + mask_names[i]
            rec_mask = masks[i]['rec'] * rec_surf
            lig_mask = masks[i]['lig'] * lig_surf
            new_mask = {'rec':rec_mask,
                        'lig':lig_mask,
                        'name':name,
                        'rec_comment':100.*sum(rec_mask)/nRec,
                        'lig_comment':100.*sum(lig_mask)/nLig }
            surf_masks += [ new_mask ]
            surf_mask_names += [ name ]

    except Exception, why:
        print "ERROR: ", lastError()
        print "name ", name, " i ", i

    ## replace comment in first mask (reference)
    surf_masks[0]['rec_comment'] = nRec
    surf_masks[0]['lig_comment'] = nLig

    return surf_masks, surf_mask_names


def castTraj( traj_1, traj_2 ):
    """
    Make two similar trajectories equal with respect to their atoms.
    """
    aInd = traj_1.getRef().compareAtoms( traj_2.getRef() )

    traj_1 = traj_1.takeAtoms( aInd[0] )
    traj_2 = traj_2.takeAtoms( aInd[1] )

    return traj_1, traj_2


def test():
    options = _defOptions()
    dir = '/home/Bis/johan/interfaces/c06/dock_multi_0919/hex1008/'
    dir_rec = dir + 'rec/1BVL'
    dir_lig = dir + 'lig/3LZT'
    options['ref'] = '/home/Bis/johan/interfaces/c06/com_wet/ref.complex'
    options['cl'] = dir +'complexes_cont.cl'
    options['r'] = [ dir_rec +'_001_hex.pdb', dir_rec +'_002_hex.pdb',
                     dir_rec +'_003_hex.pdb', dir_rec +'_004_hex.pdb',
                     dir_rec +'_005_hex.pdb', dir_rec +'_006_hex.pdb',
                     dir_rec +'_007_hex.pdb', dir_rec +'_008_hex.pdb',
                     dir_rec +'_009_hex.pdb', dir_rec +'_010_hex.pdb',
                     dir_rec +'_011_hex.pdb' ]
    options['l'] = [ dir_lig + '_001_hex.pdb', dir_lig + '_002_hex.pdb',
                     dir_lig + '_003_hex.pdb', dir_lig + '_004_hex.pdb',
                     dir_lig + '_005_hex.pdb', dir_lig + '_006_hex.pdb',
                     dir_lig + '_007_hex.pdb', dir_lig + '_008_hex.pdb',
                     dir_lig + '_009_hex.pdb', dir_lig + '_010_hex.pdb',
                     dir_lig + '_011_hex.pdb' ]
    options['additional_profile'] = [3,2]
    return options
  

###############
# main function

if len(sys.argv) < 2:
    _use()
    
options = cmdDict( _defOptions() )
#    options = test()

## paths ...
cFile = os.path.abspath( options['cl'] )
dock_dir = os.path.split( absfile(options['cl']))[0]

rec_pdbs = toList( options['r'] )
lig_pdbs = toList( options['l'] )
ref_rec_file = os.path.split(rec_pdbs[0])[0] + '/ref_rec.pdb'
ref_lig_file = os.path.split(lig_pdbs[0])[0] + '/ref_lig.pdb'
if ref_rec_file in rec_pdbs:
    rec_pdbs.remove( ref_rec_file )
if ref_lig_file in lig_pdbs:
    lig_pdbs.remove( ref_lig_file )
    
## file names for output ...
cName = nameFromPath( cFile ) + '-' + options['key']
dotPlotName = os.path.dirname(cFile) + '/' + cName + '_dot_plot.eps'
logPlotName = os.path.dirname(cFile) + '/' + cName + '_log_plot.eps'
diffPlotName = os.path.dirname(cFile) + '/' + cName + '_diff_plot.eps'
rankPlotName = os.path.dirname(cFile) + '/' + cName + '_rank_plot.eps'
contourPlotName = os.path.dirname(cFile) + '/' + cName + '_contour_plot.eps'
profilePlotName = os.path.dirname(cFile) + '/' + cName + '_plot_profile.eps'
outFile = os.path.dirname(cFile) + '/' + cName + '_best.txt'
dumpFile = os.path.dirname(cFile) + '/' + cName + '.dat'


#############################################################
## RMSD analysis

flushPrint('Performing RMSD analysis...\n')

## Load Trajectories - read and sort bound ligand and recepror trajectories
traj_rec = Trajectory( rec_pdbs )
traj_lig = Trajectory( lig_pdbs )

## reference (bound) rec and lig
ref = load( absfile( options['ref'] ) )
ref.rec().remove( ref.rec().maskH() )
ref.lig_model.remove( ref.lig().maskH() )
ref.lig_transformed = None
ref.rec().writePdb( ref_rec_file )
ref.lig().writePdb( ref_lig_file )

## Load Trajectories - read and sort bound ligand and recepror trajectories
traj_rec_ref = Trajectory( [ ref_rec_file ] )
traj_lig_ref = Trajectory( [ ref_lig_file ] )

traj_rec, traj_rec_ref = castTraj( traj_rec, traj_rec_ref )
traj_lig, traj_lig_ref = castTraj( traj_lig, traj_lig_ref )

refCom = ProteinComplex( traj_rec_ref.ref, traj_lig_ref.ref )

traj_rec = traj_rec_ref.concat( traj_rec )
traj_lig = traj_lig_ref.concat( traj_lig )

## Add surface atom profile
flushPrint('Adding accessible surface profiles...\n')
rec_asa = PDBDope( traj_rec.getRef() )
rec_asa.addASA()

lig_asa = PDBDope( traj_lig.getRef() )   
lig_asa.addASA()

## Atom sets to use for analysis
## Sets are(only heavy atoms):
masks_all , mask_all_names = allMasks( traj_rec, traj_lig, refCom )

## Masks restricted to surface atoms
masks_surf , mask_surf_names = surfMasks( traj_rec, traj_lig, refCom )

## pairwise rmsd calculation
rmsd_all_dic = calcAllRmsds( traj_rec, traj_lig,  masks_all )
rmsd_surf_dic = calcAllRmsds( traj_rec, traj_lig, masks_surf )

## rmsd info to logfile
outf = open(outFile, 'w')
outf.writelines('----------------------- ' + cName + ' ------------------------\n\n')

msk=['ALL', 'CONT', 'NONE_CONT', 'BACKBONE', 'SIDECHAIN']
msk_surf = ['SURF_ALL', 'SURF_CONT', 'SURF_NONE_CONT', 'SURF_BACKBONE', 'SURF_SIDECHAIN']

outf.writelines( ['ALL ATOMS            '] + [ '%5s'%i for i in range(1, 12) ] + ['\n'] )
for m in msk:
    c, rec_rms, lig_rms = getRmsdValues( rmsd_all_dic, m )
    rec_txt = [ '%4s '%(rec_rms[i][0]) for i in range( len(rec_rms) ) ]
    outf.writelines(['%15s rec b: '%m] + rec_txt + ['\n'] )
    rec_txt = [ '%4s '%(rec_rms[i][1]) for i in range( len(rec_rms) ) ]
    outf.writelines(['%15s rec f: '%' '] + rec_txt + ['\n'] )
    lig_txt = [ '%4s '%(lig_rms[i][0]) for i in range( len(lig_rms) ) ]
    outf.writelines(['%15s lig b: '%m] + lig_txt + ['\n'] )
    lig_txt = [ '%4s '%(lig_rms[i][1]) for i in range( len(lig_rms) ) ]
    outf.writelines(['%15s lig f: '%' '] + lig_txt + ['\n'] )    

outf.writelines( ['\n\nSURF ATOMS           '] + [ '%5s'%i for i in range(1, 12) ] + ['\n'] )
for m in msk_surf:
    c, rec_rms, lig_rms = getRmsdValues( rmsd_surf_dic, m )
    rec_txt = [ '%4s '%(rec_rms[i][0]) for i in range( len(rec_rms) ) ]
    outf.writelines(['%15s rec b: '%m] + rec_txt + ['\n'] )
    rec_txt = [ '%4s '%(rec_rms[i][1]) for i in range( len(rec_rms) ) ]
    outf.writelines(['%15s rec f: '%' '] + rec_txt + ['\n'] )
    lig_txt = [ '%4s '%(lig_rms[i][0]) for i in range( len(lig_rms) ) ]
    outf.writelines(['%15s lig b: '%m] + lig_txt + ['\n'] )
    lig_txt = [ '%4s '%(lig_rms[i][1]) for i in range( len(lig_rms) ) ]
    outf.writelines(['%15s lig f: '%' '] + lig_txt + ['\n'] )  

##############################################################

flushPrint('Loading complex list \n')
cList = load( options['cl'] )

## extract data from dictionary
flushPrint('Extracting %s data  \n'%options['key'])
inverse_data = int( options['inv'])
#data = data3DList( cList, key=options['key'], inverse=inverse_data,
#                   rm=range(1,12), lm=range(1,12), soln=512 )
nr_rec = len(rec_pdbs)
nr_lig = len(lig_pdbs)
data = data3DList( cList, key=options['key'], inverse=inverse_data,
                   rm=range(1,nr_rec+1), lm=range(1,nr_lig+1),
                   soln=len(cList)/(nr_rec*nr_lig) )





## sort data
sorted_data = sort( data, 2 )

## get the lowest number of zeros in any array
lenZeros = leastNumberOfZeros( sorted_data )
 
colorDic = {'x':'red', 'r':'cornflower blue', 'l':'pale green', 'a':'grey54'}
biggles.configure('fontsize_min', 0.7 )

#################
## EXTRA PLOTS
hexE = cList.valuesOf( 'hex_etotal' )
order = argsort( hexE )
soln = 512

def singleDataPlot( data, inverse, sorted_ref, zeros, soln, color ):
    sorted_data = sort( data )
    if inverse:
        data = (1./array(data)).tolist()
        sorted_data = sort( data )
    # log plot
    log_p = logPlot( [[sorted_data]], zeros, color, ref=None, ave=5 )
    # diff plot
    x = range( 1, soln-(zeros-1) )
    x.reverse()
    diff = ( sorted_data - sorted_ref )[-(soln-zeros):]
    diff_p = biggles.Curve( x, cumsum(diff), color=color )
    # rank plot
    x = arange( 1, soln+1, 1, 'f' )/soln
    rank_p = biggles.Curve( x, maxVal(data), color=color)
    return log_p, diff_p, rank_p

## the best 512 solutions of all according to hex energy
data_512 = take( cList.valuesOf( options['key'] ), order[:512] )
#sorted_data_512 = sort( data_512 )

log_512, diff_512, rank_512 = singleDataPlot( data_512, 
                                              inverse_data,
                                              sorted_data[0][0], lenZeros,
                                              soln, 'medium purple')
    
## mimic hex multi model docking
data_5120 = take( cList.valuesOf( options['key'] ), order[:5120] )
take_range = range( len(data_5120)-1, 0, -10 )
take_range.reverse()
sorted_data_5120 = sort( data_5120 )
sorted_data_5120 = take( sorted_data_5120, take_range )
data_5120 = take( data_5120, take_range )

log_5120, diff_5120, rank_5120 = singleDataPlot( data_5120, 
                                                 inverse_data,
                                                 sorted_data[0][0], lenZeros,
                                                 soln, 'dark orange')

###########
## Log PLOT
p_fnc = biggles.FramedPlot()
p_fnc.title = cName
p_fnc.ylog = 1
p_fnc.xlabel = 'Fraction of native contacts'
p_fnc.ylabel = 'Fraction of solutions'

p_lst = logPlot( sorted_data, lenZeros, colorDic, ref=(1,1), ave=5 )

for p in p_lst:
    p_fnc.add( p )
    
p_fnc.add( log_512[0] )
p_fnc.add( log_5120[0] )

## plot and save
p_fnc.show()
p_fnc.write_eps( logPlotName, width="18cm", height="29cm" )


############
## DIFF PLOT
p_diff = biggles.FramedPlot()
p_diff.title = cName
p_diff.xlabel = 'Solution'
p_diff.ylabel = 'Quality compared to reference'

diff_plots, diff_area = diffPlot( sorted_data, lenZeros, colorDic,
                                  type='fnc', ref=(1,1) )

print '### SOLUTIONS (sorted data)'
msg, best_soln, best = dockingInfo( diff_area, sorted_data )
outf.writelines('\n\n### SOLUTIONS (sorted data)\n')
outf.writelines( msg )
print msg

for p in diff_plots:
    p_diff.add( p )

p_diff.add( diff_512 )
p_diff.add( diff_5120 )

## plot and save
p_diff.show()
p_diff.write_eps( diffPlotName, width="18cm", height="29cm" )


###########
## RANK PLOT 
p_rank = biggles.FramedPlot()
p_rank.title = cName
p_rank.ylabel = 'Maximal fraction of native contacts'
p_rank.xlabel = 'Solution coordinate'

rank_plots, rank_area = plotRanking( data, colorDic, ref=(1,1) )

print '### RANKING'
msg, best_rank, best = dockingInfo( rank_area, data )
outf.writelines('\n\n### RANKING\n')
outf.writelines( msg )
outf.close()
print msg

for p in rank_plots:
    p_rank.add( p )

p_rank.add( rank_512 )
p_rank.add( rank_5120 )

p_rank.show()
p_rank.write_eps( rankPlotName, width="18cm", height="29cm" )


#########
## Contour plot 

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
plots += plotRmsLabels( 11, rmsd_all_dic, 'CONT')
#plots += plotNrLabels( 11, pos=-0.05 ) 

for p in plots:
    p_cont.add( p )

p_cont.show()
p_cont.write_eps( contourPlotName )


#########
## Visualize plot (dot plot)

p_model = biggles.FramedPlot()
p_model.x1.draw_ticklabels = 0
p_model.y1.draw_ticklabels = 0
p_model.title = cName

## labels
p_model.x2.label_style['size'] = 2
p_model.y2.label_style['size'] = 2
p_model.x2.label = 'Receptor model'
p_model.y2.label = 'Ligand model'
p_model.x1.label = '  '
p_model.y1.label = '  '

models = len( diff_area )

p_model.xrange = ( 0, models+1 )
p_model.yrange = ( 0, models+1 )

plots = visualizePlot( diff_area,
                       pos_symbol='filled circle',
                       neg_symbol='filled circle',
                       pos_color='pink',
                       neg_color='light blue') +\
        visualizePlot( rank_area,                       pos_symbol='circle',
                       neg_symbol='circle',
                       pos_color='red',
                       neg_color='blue')

## add best label
rank_label = resize( ['0'],(models,models) )
rank_label[best_rank[0]-1][best_rank[1]-1] = 'r'

soln_label = resize( ['0'],(models,models) )
soln_label[best_soln[0]-1][best_soln[1]-1] = 's'
plots += plotMatrixLabels( rank_label, offset=0.02 )
plots += plotMatrixLabels( soln_label, offset=0.03 )

plots += plotRmsLabels( 11, rmsd_all_dic, 'CONT')

for p in plots:
    p_model.add( p )

p_model.show()
p_model.write_eps( dotPlotName )


#######################
## dump plot data
#dump( {'diff_area':diff_area, 'rank_area':rank_area,
#       'rmsd_all_dic':rmsd_all_dic, 'rmsd_surf_dic':rmsd_surf_dic} ,dumpFile )


#################################################
##### profile plots of best solutions
plots = 4

## add aditional plot as specified in options
additional = options.get('additional_profile',None)
if additional:
    aRec, aLig = toIntList( additional )
    plots = 5

b=biggles.FramedArray(plots,1)
b.title = cName
b.xlabel = 'solution'
b.ylabel = 'fraction of native contacts'

## xray
b[0,0].add( biggles.Curve(range(512),data[0][0]) )
b[0,0].add( biggles.PlotLabel( .5, .95 ,'Xray', size=4) )

## best solutions
r, l = best_soln[0], best_soln[1]
b[1,0].add( biggles.Curve(range(512),data[r-1][l-1]) )
label = 'the highest number of good scoring solutions (%i/%i)'%(r,l)
b[1,0].add( biggles.PlotLabel( .5, .95 , label, size=4) )

## best rank
r, l = best_rank[0], best_rank[1]
b[2,0].add( biggles.Curve(range(512),data[r-1][l-1]) )
label = 'the largest number of high rankng solutions (%i/%i)'%(r,l)
b[2,0].add( biggles.PlotLabel( .5, .95 , label, size=4) )

## best overall
r, l = best[0], best[1]
b[3,0].add( biggles.Curve(range(512),data[r-1][l-1]) )
label = 'contains the best overall scoring solutions (%i/%i)'%(r,l)
b[3,0].add( biggles.PlotLabel( .5, .95 ,label , size=4) )

## additional plot
if additional:
    b[4,0].add( biggles.Curve(range(512),data[aRec-1][aLig-1]) )
    label = 'additional plot (%i/%i)'%(aRec, aLig)
    b[4,0].add( biggles.PlotLabel( .5, .95 , label, size=4) )

b.show()
b.write_eps( profilePlotName, width="18cm", height="29cm" )




            
    

