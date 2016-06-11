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

## Shorter Replacement of a_rmsd_vs_docking_performance
## creates one plot and two data files for external plotting

from Numeric import *
from Biskit import *
from Biskit.tools import *
import biggles as B
from Biskit.hist import density

def use():
    print """
    Create figure for multidock paper with the change in docking performance
    plotted against the change in rms_to_bound (both relative to the free-free
    docking).

    k .. PDBModel.info key to be plotted (rms_to_bound, rmsCA_to_bound, ..)
    o .. file name of output eps
    i .. one or many complexists with 
    Default options:
    """
    o = default_options()
    for k in o:
        print "\t",k,"\t",o[k]
        
    sys.exit(0)

def default_options():
    return { 'i':'complexes_cont.cl', 'k':'rmsIF_to_bound',
             'eps':'rms_vs_performance.eps',
             'md':'rms_vs_performance_MD.txt',
             'pcr':'rms_vs_performance_PCR.txt'}


def docking_scores( cl ):
    """
    Calculate performance score for each combination of models.
    cl - ComplexList
    -> array( N_model1 x N_model2 ) of float, model1_numbers, model2_numbers
    """
    rec_nr = cl.valuesOf( 'model1', unique=1 )
    lig_nr = cl.valuesOf( 'model2', unique=1 )

    rec_nr.sort()
    lig_nr.sort()

    r = zeros( (len(rec_nr), len(lig_nr)), 'f')

    for i in range( len(rec_nr) ):

        l1 = cl.filter('model1', rec_nr[i] )

        for j in range( len(lig_nr) ):

            l = l1.filter('model2', lig_nr[j] )
            v = l.valuesOf( 'fnac_10' )

            ## calculate score
            r[i][j] = sum( compress( greater( v, 0.1), v)**2 )

    return r, rec_nr, lig_nr


def __model_dicts( cl ):
    """
    Extract PDBModel belonging to each model number
    -> { model_nr : PDBModel }, { model_nr : PDBModel } - rec_models,lig_models
    """
    ## first solution from every docking -> 11 * 11
    first = cl.filter('soln',1)

    rec_models = {}
    for c in first.filter('model2',1):
        rec_models[ c['model1'] ] = c.rec_model

    lig_models = {}
    for c in first.filter('model1',1):
        lig_models[ c['model2'] ] = c.lig_model

    return rec_models, lig_models


def model_infos( cl ):
    """
    Collect all available rms_to_bound values from the rec and lig_models in cl
    -> { 'rms_to_bound':[float], 'rmsIF_to_bound':[float],..}, {..} - rec , lig
    """
    rec_models, lig_models = __model_dicts( cl )
    rec_keys = rec_models.keys()
    lig_keys = lig_models.keys()
    rec_keys.sort()
    lig_keys.sort()

    rms_rec = {}
    rms_lig = {}

    for infKey in rec_models.values()[0].info.keys():
        try:
            rms_rec[infKey] = [ rec_models[k].info[ infKey ] for k in rec_keys]
            rms_lig[infKey] = [ lig_models[k].info[ infKey ] for k in lig_keys]
        except:
            errWriteln('Error extracting %s:\n'%infKey + lastError() )

    return rms_rec, rms_lig


def dScore_vs_dRms( scores, rec_rms, lig_rms ):
    """
    Prepare the plotting of delta_score vs. delta_rms.
    Assumes scores[0][0] is from free-free docking.
    Assumes rec_rms[0] is rms_to_bound for free rec.
    Assumes lig_rms[0] is rms_to_bound for free lig.
    scores  - array( n_models x n_models ) of float, docking perf. scores
    rec_rms - [float], n_models rms values of rec models
    lig_rms - [float], n_models rms values of lig models
    -> 2 arrays( 121x1, 'f' ), raveled delta scores, delta rms values
    """
##     d_scores = ravel( scores - scores[0][0] )
    d_scores = ravel( (scores - scores[0][0])/ max([scores[0][0],0.03]) )

    d_rms_rec = array( rec_rms ) - rec_rms[0]
    d_rms_lig = array( lig_rms ) - lig_rms[0]

    d_rms = []
    for drec in d_rms_rec:

        for dlig in d_rms_lig:

            d_rms += [ average( [drec, dlig] ) ]

    return d_scores, d_rms


def plot( delta_scores_1, delta_rms_1, delta_scores_2, delta_rms_2,
          feps ):
    """
    """
    p = B.FramedPlot()
    p.xlabel = r"Interface rmsd to bound relative to free"
    p.ylabel = r"Docking performance relative to free"
    
    points_1 = B.Points( delta_rms_1, delta_scores_1, type='circle',
                         symbolsize=1)
    points_1.label = 'MD'
    
    points_2 = B.Points( delta_rms_2, delta_scores_2, type='diamond',
                         symbolsize=1)
    points_2.label = 'PCR-MD'

    a = concatenate( (delta_rms_1, delta_rms_2) )
    h = density( a, 20 )
    scale = max( h[:,1] ) / max( a )
    histogram = B.Curve( h[:,0], h[:,1] * scale )
    histogram.label = "distribution"
    
    p.add( points_1 )
    p.add( points_2 )
    p.add( histogram )
    p.add( B.PlotKey( 0.73, 0.95, [ points_1, points_2, histogram ] ) )

    p.show()
    p.write_eps( feps )


def report( delta_scores, delta_rms, fout ):
    """
    """
    f = open( fout, 'w')

    f.write('dScore    <dRMS>\n')
    for i in range( len( delta_scores ) ):
        f.write('%5.2f   %5.2f\n' % (delta_scores[i], delta_rms[i]) )

    f.close()

##########
## MAIN ##

options = cmdDict( default_options() )
options['i'] = toList( options['i'] )
rmstype = options['k']

if len( sys.argv ) < 2:
    use()


all_scores_md = all_rmsd_md = []
all_scores_pcr= all_rmsd_pcr= []


for f in options['i']:

    try:
        errWriteln('Loading %s...' % f)
        cl = load( absfile( f ) )

        errWriteln('Calculating docking scores...')
        scores, rec_nr, lig_nr = docking_scores( cl )

        infos_rec, infos_lig = model_infos( cl )
        d_scores, d_rms = dScore_vs_dRms( scores, infos_rec[rmstype],
                                          infos_lig[rmstype] )

        if f.find('_pcr') != -1:
            all_scores_pcr = concatenate( (all_scores_pcr, d_scores) )
            all_rmsd_pcr = concatenate( (all_rmsd_pcr, d_rms) )
        else:
            all_scores_md = concatenate( (all_scores_md, d_scores) )
            all_rmsd_md = concatenate( (all_rmsd_md, d_rms) )

    except:
        errWriteln( 'skipping %s:\n' % f + lastError() )

errWriteln('Dumping data to text...')
report( all_scores_md, all_rmsd_md, absfile( options['md'] ) )

report( all_scores_pcr, all_rmsd_pcr, absfile(options['pcr'] ) )

errWriteln('Preparing plot...')
plot( all_scores_md, all_rmsd_md, all_scores_pcr, all_rmsd_pcr,
      absfile( options['eps']) )

