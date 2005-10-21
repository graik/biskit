#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## a_compare_rms_vs_fnc.py
##
## Compare rmsd vs. fraction nativ residue and/or atom contacts
##
## last $Author$
## $Date$

import Biskit.tools as T
from Numeric import *
import biggles
import re
import sys
import string
import Biskit.mathUtils

def _use():
    print """
a_compare_rms_vs_fnc.py: Plot interface rmsd (heavy and/or backbone) vs.
                           fraction native atom/residue contacts at
                           different cutoffs.

  creates up to 4 plots: rms_if_vs_cont.eps
                         rms_if_bb_vs_cont.eps
                         rms_if_bb_vs_rms_if.eps
                         rms_hex_vs_rms_if.eps
                           
                           
Syntax:		-i  complexes_cont.cl
                -o  str, output directory
                -v  [str], list of keys to plot
                -if     1||0 create plot of key vs. interface rmsd
                -if_bb  1||0 create plot of key vs. interface backbone rmsd

"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def defOptions():
    return {'o':'.',
            'i':'complexes_cont.cl',
            'v':['fnac_4.5', 'fnac_7.5', 'fnac_10', 'fnrc_4.5'],
            'if': 1,
            'if_bb':1 }  


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


def subPlot( data, data_keys, x_key, y_key  ):

    ## data
    x_dat = data[ keys.index( x_key ) ]
    y_dat = data[ keys.index( y_key ) ]

    ## outline data
    y_val = []
    for v in y_dat:
        if not v in y_val:
            y_val += [v]
            
    y_val = sort( y_val )

    x_max = []
    x_min = []
    for v in y_val:
        x_sub_dat = compress( v == y_dat, x_dat )
        x_max += [ max( x_sub_dat ) ]
        x_min += [ min( x_sub_dat ) ]

##     y_val = Biskit.mathUtils.runningAverage( y_val, 3 )
##     x_max = Biskit.mathUtils.runningAverage( x_max, 3 )
##     x_min = Biskit.mathUtils.runningAverage( x_min, 3 )

    ## collect plots
    sp = []
    
    sp += [ biggles.Curve( x_max, y_val, color='grey') ]
    sp += [ biggles.Curve( x_min, y_val, color='grey') ]

    sp += [ biggles.Points( x_dat, y_dat, type='dot', size=1 ) ]
    sp += [ biggles.PlotLabel( .8, .9, y_key, size=3 ) ]
    
    return sp


def test():
    options = defOptions()
    dir = '/home/Bis/raik/data/tb/interfaces/c23/dock_multi_0919/hex1008/'
    options['i'] = dir + 'complexes_cont.cl'
    options['o'] = dir
    options['v'] = ['fnac_4.5', 'fnac_7.5', 'fnac_10', 'fnrc_4.5']
    return options


    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        _use()
    
#    options = test()
    options = T.cmdDict( defOptions() )

    ##### COMPLEX LIST ####
    options['if'] = int( options['if'] )
    options['if_bb'] = int( options['if_bb'] )
    
    ## Load complex dictionary
    c_lst = T.Load( options['i'] )

    keys = T.toList(options['v']) + [ 'soln', 'rms', 'rms_if', 'rms_if_bb' ]
              
    ## get complex data
    data = [ transpose(c_lst.valuesOf('soln')) ]
    
    for key in keys[1:]:
        data +=  [ transpose(c_lst.valuesOf(key)) ]

    for i in range( len(data) ):
        if sum( [ j==None for j in data[i] ] ) > 0:
            print '\n################################'
            print 'Info data for %s == None for at least one complex'%keys[i]
            print 'Plotting will fail for this key.'
            print '################################\n'
        else:
            print 'Data for %s is OK' %keys[i]

    if options['if']:
        ## Array plot - rmsd vs. fnc
        ap = biggles.FramedArray( 2, 2, title=nameFromPath( options['i'] ) )
        ap.uniform_limits = 1

        ## labels
        ap.ylabel = 'fraction of contacts'
        ap.xlabel = 'interface rmsd'

        plot_keys = T.toList( options['v'] )
        x_key = 'rms_if'

        for plot in subPlot( data, keys, x_key, 'fnac_4.5' ):
            ap[0,0].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnac_7.5' ):
            ap[0,1].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnac_10' ):
            ap[1,0].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnrc_4.5' ):
            ap[1,1].add( plot )

        ap.show()
        ap.write_eps( options['o']+'/rms_if_vs_cont.eps',
                      width="18cm", height="29cm" )

        
    if options['if_bb']:
        ## Array plot - rmsd vs. fnc
        ap2 = biggles.FramedArray( 2, 2, title=nameFromPath( options['i'] ) )
        ap2.uniform_limits = 1

        ## labels
        ap2.ylabel = 'fraction of contacts'
        ap2.xlabel = 'interface bb rmsd'

        plot_keys = [ 'fnac_10_7.5', 'fnac_7.5', 'fnac_10', 'fnrc_4.5' ]
        x_key = 'rms_if_bb'

        for plot in subPlot( data, keys, x_key, 'fnac_4.5' ):
            ap2[0,0].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnac_7.5' ):
            ap2[0,1].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnac_10' ):
            ap2[1,0].add( plot )

        for plot in subPlot( data, keys, x_key, 'fnrc_4.5' ):
            ap2[1,1].add( plot )

        ap2.show()
        ap2.write_eps( options['o']+'/rms_if_bb_vs_cont.eps',
                       width="18cm", height="29cm" )


    if options['if'] and options['if_bb']:
        ## hex rmsd vs interface rmsd
        p_rms = biggles.FramedPlot()
        p_rms.title = nameFromPath( options['i'] )
        p_rms.ylabel = 'interface rmsd'
        p_rms.xlabel = 'interface bb rmsd '

        x_dat = data[ keys.index( 'rms_if_bb' ) ]
        y_dat = data[ keys.index( 'rms_if' ) ]

        p_rms.add( biggles.Points( x_dat, y_dat, type='dot', size=1 ) )
        p_rms.show()
        p_rms.write_eps( options['o']+'/rms_if_bb_vs_rms_if.eps' )


    if options['if']:
        ## hex rmsd vs interface rmsd
        p_rms2 = biggles.FramedPlot()
        p_rms2.title = nameFromPath( options['i'] )
        p_rms2.ylabel = 'interface rmsd'
        p_rms2.xlabel = 'Hex rmsd '

        x_dat = data[ keys.index( 'rms' ) ]
        y_dat = data[ keys.index( 'rms_if' ) ]

        p_rms2.add( biggles.Points( x_dat, y_dat, type='dot', size=1 ) )
        p_rms2.show()
        p_rms2.write_eps( options['o']+'/rms_hex_vs_rms_if.eps' )
