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
## hex2complex
##

from Biskit.Dock import HexParser
from Biskit import CommandLine
import Biskit.settings as S

from Biskit.oldnumeric import *
import string
import biggles

from Biskit.tools import *
import os.path as osp


def testSuccess():
    import os
    return os.path.exists( options.test['o'] )

def setOptions():

    o = CommandLine( rec=str, lig=str, hex=str,
                     o='complexes.cl',
                     p=True,
                     mac=False )

    o.setDescription("""
Parse output file from hex docking run, create ComplexList and pickle it to
a file. Also creates a plot of the cluster distribution (using hex rmsd).

Note: Input model dictionaries are created by selectModels.py.

Syntax:    hex2complex -rec |models_rec| -lig |models_lig| -hex |hex.out|
-o |output name| -p |create plot|

Example:    hex2complex -rec 1BZY_models.dic -lig 2AKZ_models.dic
-hex 1BZY_2AKZ_hex.out -o complexes.cl -p
""")

    o.setRequired( [ 'rec', 'lig', 'hex' ] )

    o.setDocs( rec='pickled receptor model dictionary',
               lig='pickled ligand model dictionary',
               hex='output file from hex',
               o=  'name of resulting complex list',
               p=  'create biggles plot of rmsd vs. solution',
               mac="force rec and lig 'model 1' to be used for all" )

    f = testRoot()
    ftemp = S.tempDirLocal

    o.setTest( lig= osp.join(f, 'dock/lig/1A19_model.dic'),
               rec= osp.join(f, 'dock/rec/1A2P_model.dic'),
               hex= osp.join(f, 'dock/hex/1A2P-1A19_hex.out'),
               o=   osp.join(ftemp, 'complexes.cl'),
               p= 0,
               fsuccess=testSuccess )
    

    o.setTestCleanup( [ osp.join(ftemp, 'complexes.cl'),
                        osp.join(ftemp, 'complexes.eps') ] )

    o.parse()

    return o

## publish options record globally so that Test module can find it
options = setOptions()


def plot( complex_lst ):
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


def main( o ):
    """@param o: options dictionary"""

    ## load pickeled model dictionaries with PCRModels indexed by hex model
    ## number
    rec_lst = load( o['rec'] )
    lig_lst = load( o['lig'] )

    ## open hex output file
    if o['mac']:
        parser = HexParser( o['hex'], rec_lst, lig_lst, forceModel=(1,1) )
    else:
        parser = HexParser( o['hex'], rec_lst, lig_lst )

    ## generate dictionary of Complex objects from hex output
    complex_lst = parser.parseHex()

    ## pickle list to file
    dump( complex_lst, o['o'] )

    return complex_lst

#################################

if __name__ == '__main__':

    complex_lst = main(options)

    if options['p']:
        plot( complex_lst )

    if 'test' in options:
        print 'Test worked? ', options.testSuccess()
        options.testCleanup()
