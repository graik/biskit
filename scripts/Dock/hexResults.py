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
from string import *
import biggles
from Biskit.tools import *


def _use():
    print """
hexResults   Get info about docking results from one or more complexGroup files. 
Syntax       hexResult -cg |complexGroup.cg| [ -p |plot name| -o |file name| ]
            
Result       Plot and report file

Default values:
    """
    default = _defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def _defOptions():
    return { 'cg':'complexes_grouped.cg',
             'p':'compareHex.eps',
             'o':'compareHex.txt'}

def nameFromPath( path ):
    """
    Extract a nice name from the absolute complexes_grouped.cg path
    """
    path = string.split( path, '/' )
    name = ''
    for p in path:
        if re.match( '^[a,c][0-9]{2}', p ):
            name += re.findall( '^[a,c][0-9]{2}', p )[0] + '-'
        if re.match(  '^dock_[0-9,x,r,a,y]{4}', p ):
            try:
                name += re.findall( '^dock_[0-9,x,r,a,y]{4}.+', p )[0]
            except:
                name += re.findall( '^dock_[0-9,x,r,a,y]{4}', p )[0]
            
    return name


def clustInfo( complexesGrouped ):
    
    ## get cluter and data lists
    clst_list = []

    list = { 'rms':[], 'fnrc_4.5':[],
             'fractNatSurf_lig':[], 'fractNatSurf_rec':[],
             'cons_abs':[], 'cons_max':[], 'cons_ent':[],
             'hex_etotal':[], 'hex_eshape':[],
             'ePairScore':[], 'eProsa_surf':[], 'eProsa_pair':[] }
     
    for cg in complexesGrouped:
        clst_list += [ cg.info['cluster'] ]

        for key in [ 'rms', 'fnrc_4.5', 'fractNatSurf_lig',
                     'fractNatSurf_rec', 'cons_abs', 'cons_max',
                     'cons_ent', 'hex_etotal', 'hex_eshape', 'ePairScore' ]:
            list[key] += [ cg.info[ key ] ]

        list['eProsa_surf'] += [ cg.info['eProsa'][1] ]
        list['eProsa_pair'] += [ cg.info['eProsa'][0] ]
        
    ## get average, max, min and size of cluster
    data = { 'rms':[], 'fnrc_4.5':[], 'fractNatSurf_lig':[], 'fractNatSurf_rec':[],
             'cons_abs':[], 'cons_max':[], 'cons_ent':[], 'hex_etotal':[], 
             'hex_eshape':[], 'ePairScore':[], 'eProsa_surf':[], 'eProsa_pair':[] }
    
    clst_range = range( 1, max(clst_list)+1 )

    for key in data.keys():  
        for clst in clst_range:
            d = compress( equal(clst_list, clst), list[key] )
            data[key] += [ [ average( d ), max( d ), min ( d ), len( d ) ] ]      
        data[key] = transpose(data[key])

    return data, clst_range


def plotClusterSize( data ):
    """
    cluster vs. size inset
    """
    inset =  biggles.FramedPlot()

    inset.frame.draw_ticks = 0
    inset.x1.draw_ticklabels = 0
    inset.y1.draw_ticklabels = 0 
    inset.y2.draw_ticklabels = 1 
    inset.y2.draw_ticks = 1
    inset.y2.ticks_style['color'] = 'red'
 
    inset.add( biggles.Curve( clst_range, data['fnrc_4.5'][3], color='red' ) )

    return inset


def plotKey( data, name, key ):
    ## plot the cluster rmsd dist.
    plot = biggles.FramedPlot()
    plot.add( biggles.PlotLabel( 0.5, 0.90, name, size=4 ) )
    
    ## customize plot appearence
    #plot.x1.label = 'cluster'
    plot.y2.label = 'size'
    plot.y1.label = key

    plot.x2.draw_ticks = 0
    plot.y2.draw_ticks = 0

    ## Plot
    plot.add( biggles.ErrorBarsY(clst_range , data[key][1], data[key][2], width=1 ) )
    plot.add( biggles.Points( clst_range, data[key][0], type='cross', size=1 ) )
    plot.add( biggles.LineY( 0.25, type='dot' ) )

    ## add label with info about 'good' solutions
    good = []
    for clst in clst_range:
        if data[key][0][ clst -1 ] > 0.25:
            good += [ clst ]
    plot.add( biggles.PlotLabel( 0.5, 0.80, str(good), size=1 ) )

    return plot


def bestClust( data, key = 'fnrc_4.5' ):
    """
    the cluster with the highest average fraction of native contacts
    """
    max_ave = max( data[key][0] )
    mask = equal(  data[key][0], max_ave )

    return ravel(compress( mask, data[key] ))
    

def dataDist( data, key = 'fnrc_4.5' ):
    """
    the distribution of 'correct' clusters
    """
    fractions = [ .1, .2, .3, .4, .5, .6, .7, .8, .9 ]
    result = {}
    
    for f in fractions:
        name = str( int(f*100) )
        r = compress( greater( data[key][0], f ), data[key][3] )
        if r:
            result['c'+name] = len(r)
            result['a'+name] = average(r)
            result['n'+name] = min(r)
            result['x'+name] = max(r)
        else:
            result['c'+name] = 0
            result['a'+name] = 0
            result['n'+name] = 0
            result['x'+name] = 0
    return result


def clusterDist( data, name ):
    """
    create a statistics text block for plotting
    """
    best = bestClust( data )
    info = {'bestAve':best[0], 'bestMax':best[1], 'bestMin':best[2], 'bestSize':best[3] }
    
    info.update( { 'name':name } )

    dist = dataDist( data )
    info.update( dist )

    body = """
%(name)17s:  Best clust: average %(bestAve)4.2f - %(bestSize)i members
                                max/min %(bestMax)4.2f /%(bestMin)5.2f
      10   20   30   40   50   60   70   80   90
nr  %(c10)4i %(c20)4i %(c30)4i %(c40)4i %(c50)4i %(c60)4i %(c70)4i %(c80)4i %(c90)4i 
ave %(a10)4.1f %(a20)4.1f %(a30)4.1f %(a40)4.1f %(a50)4.1f %(a60)4.1f %(a70)4.1f %(a80)4.1f %(a90)4.1f 
min %(n10)4i %(n20)4i %(n30)4i %(n40)4i %(n50)4i %(n60)4i %(n70)4i %(n80)4i %(n90)4i 
max %(x10)4i %(x20)4i %(x30)4i %(x40)4i %(x50)4i %(x60)4i %(x70)4i %(x80)4i %(x90)4i 
""" %( info )

    return body


def clusterInfo( data, name ):
    """
    create a statistics text block for plotting
    """
    info = {'name': name,
            'nrClust':len(data['fnrc_4.5'][0]),
            'bestSoln':max( data['fnrc_4.5'][1] ),
            'avgSize':average( data['fnrc_4.5'][3]),
            'singletons':sum(equal(data['fnrc_4.5'][3], 1))}

    body = """
%(name)17s: Nr clusters: %(nrClust)4i
                   Best solution: %(bestSoln)4.2f
                   Average size: %(avgSize)4.1f
                   Singletons: %(singletons)i
 """ %( info )

    return body


def plotInfo( info ):
    """
    biggles FramedArray with only information labels
    """
    stat = biggles.FramedPlot()
    #stat.title = 'grouping info'
    
    ## turn off drawing of all axis related objects
    ## note: if all is turned off plotting doesnt work, so make one white
    stat.x.draw_ticks = 0
    stat.y.draw_ticks = 0
    stat.x.draw_ticklabels = 0
    stat.y.draw_ticklabels = 0
    stat.y.draw_spine = 0
    stat.x.spine_style['color'] = 'white'

    ## have to make it a plot - add a white line
    stat.add( biggles.LineY( 0.01, type='dot', color='white' ) )

    ## add info lines one by one (from bottom and up!)
    info = string.split( info , '\n' )
    for l in range(len(info)):
        yPos = 0.95 - 0.04*l
        stat.add( biggles.PlotLabel( .05, yPos, info[l], halign='left', \
                                     size=1, face='courier' ) )
    return stat


def test():
    options = _defOptions()
    dir = '/home/Bis/johan/interfaces/'
    options['cg'] = [ dir +'c01/dock_0519_pc2/hex01/complexes_grouped.cg',
                      dir + 'c01/dock_0423/hex01/complexes_grouped.cg',
                      dir + 'c01/dock_xray/hex01/complexes_grouped.cg' ]
    return options


################
# main function

if __name__ == '__main__':
    if len(sys.argv) < 2:
        _use()
    
#    options = test()
    options = cmdDict( _defOptions() )
    
## options and variables
cgList = toList( options['cg'] )

plotList = []
dist = ''
info = ''

## table
t = biggles.Table( len(cgList)+1 ,1)
t.cellspacing = 1.2 #default 2.0

biggles.configure('fontsize_min', 0.7 )

for cg in cgList:
    ## path and name
    cgFile = os.path.abspath( cg )
    cgName = nameFromPath(cgFile)
    print 'Working on ', cgName
    
    ## collect data
    group = load( cgFile )
    data, clst_range = clustInfo( group )

    ## create framedArray
    plot = plotKey( data, cgName, 'fnrc_4.5' )
    plot.add( biggles.Inset( (0.0,0.0), (1.0,1.0), plotClusterSize( data ) ) )
    plotList += [ plot ]

    ## gather cluster information
    dist += clusterDist( data, cgName )
    info += clusterInfo( data, cgName )
    
## grouping info
t2 = biggles.Table( 1, 2 )
stat1 = plotInfo( info )
stat2 = plotInfo( dist )
t2[0,0] = stat2
t2[0,1] = stat1

## add rfamedArrays to table
i = 0
for p in plotList:
    t[i,0] = p
    i += 1
    
t[i,0] = t2

## show and plot
t.show()
t.write_eps( options['p'], width="18cm", height="29cm" )


## Info to text file
outFile = open( options['o'], 'w' )
outFile.write('# --------------- ' +options['o']  + ' -------------------\n')
outFile.write(info)
outFile.write(dist)
outFile.close()

print 'Done \n'



