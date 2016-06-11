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
from Biskit.mathUtils import SD
from Biskit.Statistics.lognormal import logConfidence

def use():
    print """
    a_random_contacting.py -i 1.cl 2.cl 3.cl .. -ref ref.complex
                           -nat natively_contacted.cl [ -t ]
                           [ -nout summary_output_file
                             -rout random_output_file ]
    
    Get confidence of native scores from several scores to random reference.

    Prints table 3 of multidock paper, i.e. fnac, score, rms, ..
    of free vs. free docking, the docking with highest fnac, and the docking
    with highest score - FOR EACH complex list. The ref.complex is used to
    calculate the interface rmsd to the bound.
    Calculates averages and confidence of highest score, rms..
    The 'real' table 3 line for the native is appended to a separate file.
    
    t   .. print table header
    nout .. append line with native fnac, scores,.. and confidence to this file
    rout .. append lines with random fnac, scores,.. to this file [STDOUT]
    
    Default options:
    """
    o = default_options()
    for k in o:
        print "\t",k,"\t",o[k]
        
    sys.exit(0)

def default_options():
    return { 'i':'complexes_cont.cl',
             'ref':'../../com_wet/ref.complex'}


def docking_scores( cl, soln=None ):
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

    if soln:
        cl = cl.filter('soln', (0, soln) )

    for i in range( len(rec_nr) ):

        l1 = cl.filter('model1', rec_nr[i] )

        for j in range( len(lig_nr) ):

            l = l1.filter('model2', lig_nr[j] )

            v = l.valuesOf( 'fnac_10', 0. )

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

def model_dicts( cl, old_recs=None, old_ligs=None ):
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

    if old_recs:
        copy_rms_to_bound( old=old_recs, new=rec_models )

    if old_ligs:
        copy_rms_to_bound( old=old_ligs, new=lig_models )

    return rec_models, lig_models

def copy_rms_to_bound( old, new ):
    """
    -> 0 | 1, 1 if rms_to_bound could be copied from model_dic entry
    """
    for k in new:
        new[ k ].info.update( old[ k ].info )

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


def getReport( cl, scoreMatrix=None, rec_nr=None, lig_nr=None ):
    """
    Collect highest fnac, score, IF rms for docking of free-free (ff),
    docking with highest fnac (fn), and docking with highest score (sc).
    -> {'ff':{..}, 'fn':{..}, 'sc':{..}, '<sc>':float }
    """
    if not scoreMatrix:
        scoreMatrix, rec_nr, lig_nr = docking_scores( cl )
    
    ## result
    r = {'ff':None, 'fn':None, 'sc':None, '<sc>':None }

    ## free - free: I assume the lowest model number is for the free
    ff = cl.filter('model1',rec_nr[0]).filter('model2',lig_nr[0])
    v = ff.valuesOf('fnac_10')
    
    fnac = max( v )
    score= sum( compress( greater( v, 0.1), v)**2 )
    rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    r['ff'] = {'fnac':fnac,'score':score,'rms_rec':rms_rec,'rms_lig':rms_lig}

    ## best fnac
    m1, m2 = highest_fnac_models( cl )
    ff = cl.filter('model1',m1).filter('model2',m2)
    v = ff.valuesOf('fnac_10')
    
    fnac = max( v )
    score= sum( compress( greater( v, 0.1), v)**2 )
    rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    r['fn'] = {'fnac':fnac,'score':score,'rms_rec':rms_rec,'rms_lig':rms_lig,
               'm1':m1-1, 'm2':m2-1}

    ## best score
    m1, m2 = highest_score_models( scoreMatrix, rec_nr, lig_nr )
    ff = cl.filter('model1',m1).filter('model2',m2)
    v = ff.valuesOf('fnac_10')
    
    fnac = max( v )
    score= sum( compress( greater( v, 0.1), v)**2 )
    rms_rec = ff[0].rec_model.info['rmsIF_to_bound']
    rms_lig = ff[0].lig_model.info['rmsIF_to_bound']

    r['sc'] = {'fnac':fnac,'score':score,'rms_rec':rms_rec,'rms_lig':rms_lig,
               'm1':m1-1, 'm2':m2-1}

    ## average score
    v = ravel( scoreMatrix )
    avg = average( v[1:] )

    r['<sc>'] = avg

    return r

def report( d, f=sys.stdout, header=0, prefix='' ):
    """
    d - dict of dicts with highest fnac, score, rms for 3 different dockings
    f - open file handle
    """
    if header:
        if prefix:
            f.write(" "*(len(prefix)+1) )
        f.write("        FREE            ")
        f.write("         Highest Fnac             ")
        f.write("         Highest Score            \n")
        s =   "fnac #r #l RMS_r P  RMS_l P  score P  "
        if prefix:
            f.write(" "*(len(prefix)+1) )
        f.write("fnac RMS_r RMS_l score P  " + 2*s + "<sc> P\n")

    ## fill in fake confidence values
    if not 'score_p' in d['ff']:
        for ki in ['ff','sc','fn']:
            for kj in ['score','rms_rec','rms_lig']:
                d[ ki ][ kj + '_p' ] = 0.
                d[ ki ][ kj + '_+' ] = ''
        d['<sc>_p'] = 0.
        d['<sc>_+'] = ''

    if prefix:
        f.write( prefix + ' ' )
        
    f.write("%(fnac)4.2f %(rms_rec)4.2f  %(rms_lig)4.2f %(score)5.2f%(score_+)1s%(score_p)02i  "%d['ff'] )
    
    f.write("%(fnac)4.2f %(m1)2i %(m2)2i %(rms_rec)4.2f%(rms_rec_+)1s%(rms_rec_p)02i  %(rms_lig)4.2f%(rms_lig_+)1s%(rms_lig_p)02i %(score)5.2f%(score_+)1s%(score_p)02i  " % d['fn'] )

    f.write("%(fnac)4.2f %(m1)2i %(m2)2i %(rms_rec)4.2f%(rms_rec_+)1s%(rms_rec_p)02i  %(rms_lig)4.2f%(rms_lig_+)1s%(rms_lig_p)02i %(score)5.2f%(score_+)1s%(score_p)02i  " % d['sc'] )

    f.write("%(<sc>)4.2f%(<sc>_+)1s%(<sc>_p)02i\n" % d )
    f.flush()
    

def __getConfidence( nat, random, key, avgsc=1 ):
    """
    nat    - dict,  native fnac, score, rms for 3 docking runs
    random - [dict], several random fnac, score, rms for 3 docking runs
    key    - str, which of the 3 docking runs to consider ('ff','fn','sc')
    avgsc  - 0|1, also get confidence for average score of all dockings
    Adds confidence values and whether ('') or not ('-') the native is above
    mean to the native dict. New keys are:
    'score_p', 'score_+', 'rms_rec_p'/'..+', 'rms_lig_p'/'..+'
    """
    vKeys = ['score'] + (key != 'ff') * ['rms_rec', 'rms_lig']

    for k in vKeys:

        r = [ d[key][k]   for d in random ]
        x = nat[key][k]

        conf, median = logConfidence( x, r )
        nat[key][k+'_p'] = int( round( conf*100 ) )
        nat[key][k+'_+'] = ''+'-'*(x < median)

    if avgsc:
        r = [ d['<sc>'] for d in random ]
        x = nat['<sc>']

        conf, median = logConfidence( x, r )
        nat['<sc>_p'] = int( round( conf*100) )
        nat['<sc>_+'] = ''+'-'*(x < median)


def reportNative( native, R, fout, prefix='' ):
    """
    native - dict with fnacs, scores, rms .. of natively contacted list
    R      - list of dict with collected random fnacs, scores, rms..
    fout   - str, outfile name
    """
    print 70 * '-'
    av = 1

    errWrite('calculating confidence...')
    for k in ['ff', 'fn', 'sc']:

        __getConfidence( native, R, k, av )
        av = 0

        errWrite('#')

    errWrite('done\n')

    report( native, f=fout, header=0, prefix=prefix )

##########
## MAIN ##

if len( sys.argv ) < 3:
    use()

options = cmdDict( default_options() )
options['title'] = 't' in options

cl_files  = toList( options['i'] )

cl_file_native = options['nat']

maxSoln = int( options.get('max', 0 ) )

## output files
f_native_out = sys.stdout
if options.get('nout',0):
    f_native_out = open( absfile( options['nout'] ), 'a' )

f_random_out = sys.stdout
if options.get('rout',0):
    f_random_out = open( absfile( options['rout'] ), 'a' )

prefix = options.get('prefix', '')

## native reference complex
ref_nat = load( absfile( options['ref'] ) )

## collect dict with fnac, score, rms for three docking runs from each list
result = []
rec_models = lig_models = None

try:

    errWrite('\nLoading natively contacted list...')
    cl = load( absfile( cl_file_native ) )
    if maxSoln:
        cl = cl.filter('soln', (0, maxSoln) )
        
    ## calculate dock performance scores for all model combinations
    errWrite('\nCalculate docking scores...')
    scores, rec_model_numbers, lig_model_numbers = docking_scores( cl )

    ## extract rec and lig PDBModels into dicts indexed by model number
    ## copy over rms to bound
    rec_models, lig_models = model_dicts( cl, rec_models, lig_models )

    ## add RMS_TO_BOUND values to all
    ## make it work also with unequal numbers of rec and lig models
    errWrite('\nCalculating RMSD to bound for all models...\n')

    n_models = len( rec_models )
    if len( lig_models ) > n_models:
        n_models = len( lig_models )

    for i in range( n_models ):
        try: rec = rec_models[ rec_model_numbers[i] ]
        except: rec = rec_models.values()[-1]

        try: lig = lig_models[ lig_model_numbers[i] ]
        except: lig = lig_models.values()[-1]

        rms_to_bound( rec, lig, ref_nat )

    result_native= getReport(cl, scores, rec_model_numbers, lig_model_numbers )

    print "Similarity to random reference complexes"

    for cl_f in  cl_files:

        errWrite('\nLoading list...')

        cl = load( absfile( cl_f ) )
        if maxSoln:
            cl = cl.filter('soln', (0, maxSoln) )

        ## calculate dock performance scores for all model combinations
        errWrite('\nCalculate docking scores...')
        scores, rec_model_numbers, lig_model_numbers=docking_scores(cl)

        ## extract rec and lig PDBModels into dicts indexed by model number
        ## copy over rms to bound
        rec_models, lig_models = model_dicts( cl, rec_models, lig_models )

        result += [ getReport( cl, scores,
                            rec_model_numbers, lig_model_numbers ) ]

        report( result[-1],f_random_out,options['title'],prefix=prefix  )
        options['title'] = 0 ## only print title once (if at all)
    

except:
    EHandler.fatal("Error")

## calculate averages and SDs
reportNative( result_native, result, f_native_out, prefix=prefix )

for f in [ f_random_out, f_native_out ]:
    if not f is sys.stdout:
        f.close()

