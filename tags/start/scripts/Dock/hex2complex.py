#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## hex2complex
##
## last $Author$
## $Date$

from Biskit.Dock import HexParser
from Biskit.tools import *

from Numeric import *
import string
import biggles

def _use():
    print """
hex2complex:	Parse output file from hex docking run, create dictionary of
		Complex(es), and pickle it to a file.
                Creates a plot of the cluster distribution (using hex rmsd).

                rec, lig  - receptor and ligang model dictionary
                hex       - output file from hex
                o         - name of resulting complex list
                p         - create biggles plot of rmsd vs. solution
                mac       - force rec and lig 'model 1' to be used for all

Syntax:		hex2complex -rec |models_rec| -lig |models_lig| -hex |hex.out|
		-o |output name| -p |create plot|

Example:	hex2complex -rec 1BZY_models.dic -lig 2AKZ_models.dic
		-hex 1BZY_2AKZ_hex.out -o complexes.cl -p
"""
    sys.exit(0)

def main(options):
    ## load pickeled model dictionaries with PCRModels indexed by hex model
    ## number
    rec_lst = Load( options['rec'] )
    lig_lst = Load( options['lig'] )

    ## open hex output file
    if options.has_key('mac'):
        parser = HexParser( options['hex'], rec_lst, lig_lst, forceModel=(1,1) )
    else:
        parser = HexParser( options['hex'], rec_lst, lig_lst )
        
    ## generate dictionary of Complex objects from hex output
    complex_lst = parser.parseHex()

    ## pickle list to file
    Dump( complex_lst, options['o'] )

    return complex_lst

#################################

def defOptions():
    return {'o':'complexes.cl'}  

def test():
    options = defOptions()
    dir = testRoot() + '/dock/'
    options['lig'] = dir + 'lig/1A19_model.dic'
    options['rec'] = dir + 'rec/1A2P_model.dic'
    options['hex'] = dir + 'hex/1A2P-1A19_hex.out'
    options['o'] =   dir + 'hex/complexes.cl'
    options['p'] = 1

    return options
    

if __name__ == '__main__':
 #   if len(sys.argv) < 4:
 #       _use()
    
    options = test()
#    options = cmdDict( defOptions() )
    
    complex_lst = main( options )
    
    if options.has_key('p'):
        ## plot the cluster rmsd dist.
        plot = biggles.FramedPlot()
        inset =  biggles.FramedPlot()

        ## plot title
        plot.title = '/'.join(string.split(absfile(options['o']), '/')[-5:-1])
        
        ## customize plot appearence
        plot.x1.label = 'cluster'
        plot.y2.label = 'size'
        plot.y1.label = 'rmsd'

        plot.x2.draw_ticks = 0
        plot.y2.draw_ticks = 0

        inset.frame.draw_ticks = 0
        inset.x1.draw_ticklabels = 0
        inset.y1.draw_ticklabels = 0 
        inset.y2.draw_ticklabels = 1 
        inset.y2.draw_ticks = 1
        inset.y2.ticks_style['color'] = 'red'


        ## get cluter and rmsd lists
        clst_list = []
        rms_list = []
        for compl in complex_lst:
            clst_list += [ compl.info['hex_clst'] ]
            rms_list += [ compl.info['rms'] ]
    
        ## get average, max, min and size of cluster
        data = []
        clst_range = range( 1, max(clst_list)+1 )
        for clst in clst_range:
            rms = compress( equal(clst_list, clst), rms_list)
            data += [ [ average( rms ), max( rms ), min ( rms ), len( rms ) ] ]
        data = transpose(data)

        ## Inset
        inset.add( biggles.Curve( clst_range, data[3], color='red' ) )

        ## Plot
        plot.add( biggles.ErrorBarsY(clst_range , data[1], data[2] ) )
        plot.add( biggles.Points( clst_range, data[0], type='cross', size=1 ) )
        plot.add(biggles.Inset( (0.0,0.0), (1.0,1.0), inset ) )
        plot.add( biggles.LineY( 10, type='dot' ) )

        ## add label with info about 'good' solutions (average rmsd < 10A)
        good = []
        for clst in clst_range:
            if data[0][clst -1] < 10:
                good += [ clst ]
        plot.add( biggles.PlotLabel( 0.5, 0.98,
                                     'Solutions with rmsd < 10A', size=1 ) )
        plot.add( biggles.PlotLabel( 0.5, 0.95, str(good), size=1 ) )

        ## plot and save
        plot.show()
        plot.write_eps( string.split(options['o'], '.')[0] +'.eps' )