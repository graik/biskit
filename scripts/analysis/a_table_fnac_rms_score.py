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

from Numeric import *
from Biskit import *
from Biskit.tools import *

def use():
    print """
    a_table_fnac_rms_score.py [ -i complexes_cont.cl -ref ref.complex -t ]
    
    Creates one line of table 3 of multidock paper, i.e. fnac, score, rms, ..
    of free vs. free docking, the docking with highest fnac, and the docking
    with highest score.
    t .. print table header

    Default options:
    """
    o = default_options()
    for k in o:
        print "\t",k,"\t",o[k]
        
    sys.exit(0)

def default_options():
    return { 'i':'complexes_cont.cl',
             'ref':'../../com_wet/ref.complex' }


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


def highest_fnac_models( cl ):
    """
    -> (int, int) model1 and model2 of complex with highest fnac_10
    """
    c = cl[ cl.argmax( 'fnac_10' ) ]
    return c['model1'], c['model2']


def highest_score_models( scores,
                          rec_numbers=range(1,12),
                          lig_numbers=range(1,12) ):
    """
    scores - array from score()
    rec_numbers - lst of rec model numbers (default 1..11)
    lig_numbers - lst of lig model numbers
    -> (int, int) model1 and model2 of docking with highest score.
    """
    ## index of highest value in each row
    best_of_row = argmax( scores, 0 )
    ## highest value in each row
    best_values = [ scores[best_of_row[i]][i] for i in range(len( scores )) ]

    best_axis_1 = argmax( best_values )
    best_axis_0 = best_of_row[ best_axis_1 ]

    return rec_numbers[ best_axis_0 ], lig_numbers[ best_axis_1 ]

def model_dicts( cl ):
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


def rms_to_bound( rec, lig, refCom ):
    """
    Calculate all-heavy, CA, and (heavy) interface rms between the given rec
    and lig model and the bound conformation in refCom. The values are put into
    the info dict of rec and lig (but the models are otherwise not changed).
    rec - PDBModel
    lig - PDBModel
    refCom - Complex, bound reference complex
    -> int, int - if_rms_rec, if_rms_lig
    """
    ## casting of models and contact matrix
    i_rec, i_ref_rec = rec.compareAtoms( refCom.rec_model )
    i_lig, i_ref_lig = lig.compareAtoms( refCom.lig_model )
    
    ires_ref_rec = refCom.rec_model.atom2resIndices( i_ref_rec )
    ires_ref_lig = refCom.lig_model.atom2resIndices( i_ref_lig )

    cont = refCom.resContacts()  ## will give the cached matrix
    cont = take( take( cont, ires_ref_lig, 1), ires_ref_rec, 0 )

    mrec = rec.take( i_rec )
    mlig = lig.take( i_lig )

    ref_rec = refCom.rec_model.take( i_ref_rec )
    ref_lig = refCom.lig_model.take( i_ref_lig )

    ## calculate total rms
    rec.info['rms_to_bound'] = mrec.rms( ref_rec, mrec.maskHeavy() )
    lig.info['rms_to_bound'] = mlig.rms( ref_lig, mlig.maskHeavy() )
    
    rec.info['rmsCA_to_bound'] = mrec.rms( ref_rec, mrec.maskCA() )
    lig.info['rmsCA_to_bound'] = mlig.rms( ref_lig, mlig.maskCA() )

    ## get interface mask
    if_mask_rec = mrec.res2atomMask( sum( cont, 1) )
    if_mask_lig = mlig.res2atomMask( sum( cont, 0) )

    if_mask_rec = if_mask_rec * mrec.maskHeavy()
    if_mask_lig = if_mask_lig * mlig.maskHeavy()

    ## calculate interface RMS
    rms_rec = mrec.rms( ref_rec, if_mask_rec )
    rms_lig = mlig.rms( ref_lig, if_mask_lig )

    rec.info['rmsIF_to_bound'] = rms_rec
    lig.info['rmsIF_to_bound'] = rms_lig

    return rms_rec, rms_lig


def report( cl, header=0, scoreMatrix=None, rec_nr=None, lig_nr=None ):
    """
    """
    if not scoreMatrix:
        scoreMatrix, rec_nr, lig_nr = docking_scores( cl )
    
    if header:
        print "        FREE            ",
        print "         Highest Fnac             ",
        print "         Highest Score            "
        s =   "fnac #rec #lig RMS_r RMS_l score"
        print "fnac RMS_r RMS_l score  " + s + "  " + s + "  <score>"

    ## free - free: I assume the lowest model number is for the free
    ff = cl.filter('model1',rec_nr[0]).filter('model2',lig_nr[0])

    v = ff.valuesOf('fnac_10')
    ff_fnac = max( v )
    ff_score= sum( compress( greater( v, 0.1), v)**2 )
    ff_rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    ff_rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    print "%4.2f %4.2f  %4.2f  %4.2f  " % \
          (ff_fnac,ff_rms_rec,ff_rms_lig,ff_score),

    ## best fnac
    m1, m2 = highest_fnac_models( cl )
    ff = cl.filter('model1',m1).filter('model2',m2)

    v = ff.valuesOf('fnac_10')
    ff_fnac = max( v )
    ff_score= sum( compress( greater( v, 0.1), v)**2 )
    ff_rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    ff_rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    print "%4.2f %2i  %2i    %4.2f  %4.2f  %4.2f  " % \
          (ff_fnac,m1-1,m2-1,ff_rms_rec,ff_rms_lig,ff_score),

    ## best score
    m1, m2 = highest_score_models( scoreMatrix, rec_nr, lig_nr )
    ff = cl.filter('model1',m1).filter('model2',m2)

    v = ff.valuesOf('fnac_10')
    ff_fnac = max( v )
    ff_score= sum( compress( greater( v, 0.1), v)**2 )
    ff_rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    ff_rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    print "%4.2f %2i  %2i    %4.2f  %4.2f  %4.2f  " % \
          (ff_fnac,m1-1,m2-1,ff_rms_rec,ff_rms_lig, ff_score),

    ## average score
    v = ravel( scoreMatrix )
    avg = average( v[1:] )
    print "%4.2f" % avg
    

##########
## MAIN ##

options = cmdDict( default_options() )
options['title'] = 't' in options

try:
    errWrite('Loading list...')

    cl = load( absfile( options['i'] ) )
    ref= load( absfile( options['ref']))

except:
    errWrite('error.\n')
    use()
    
## calculate dock performance scores for all model combinations
errWrite('\nCalculate docking scores...')
scores, rec_model_numbers, lig_model_numbers = docking_scores( cl )

## extract rec and lig PDBModels into dicts indexed by model number
rec_models, lig_models = model_dicts( cl )

## add RMS_TO_BOUND values to all
## I try to make it work also with unequal numbers of rec and lig models
errWrite('\nCalculating RMSD to bound for all models...\n')

n_models = len( rec_models )
if len( lig_models ) > n_models:
    n_models = len( lig_models )

for i in range( n_models ):
    try: rec = rec_models[ rec_model_numbers[i] ]
    except: rec = rec_models.values()[-1]
    
    try: lig = lig_models[ lig_model_numbers[i] ]
    except: lig = lig_models.values()[-1]

    rms_to_bound( rec, lig, ref )

report( cl, options['title'], scores, rec_model_numbers, lig_model_numbers )

## Dump ComplexList with PDBModels, that have RMS values in info dict
errWrite('\nDumping ComplexList with RMS values...\n')
dump( cl, absfile( options['i'] ) )
