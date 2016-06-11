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
from Biskit.tools import *
from string import *
from Biskit.Dock.Complex import Complex as ProteinComplex
from Biskit.mathUtils import *

def _use():
    print """
Collect, calculate and save to disc (both as text files and pickled
dictionaries) various data about a complex list. This script is written
to collect data from multidocking runs and assums that the first
ligand and receptor model is the free xray structure.

Syntax:   a_rmsd_vs_dock_performance.py -cl |complexList.cl|
                         -ref |ref.complex| [-key |str| -inv [int|]

          cl - complexList, has to contain info dictionary data for key
          ref - reference complex
          key - info dictionary key to plot (high values are considered good)
          inv - 1||0 inverse data associated with key (i.e. for rmds plots)

Output:   An output directory 'model_vs_rmsd' is created and various text
          files and corresponding dictionaries are written to it.

Default values:
    """
    default = _defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def _defOptions():
    return { 'cl':'complexes_cont.cl',
             'ref':'../../ref.complex',
             'key':'fnac_10',
             'inv':0}


def nameFromPath( path ):
    """
    Extract a nice name for the plot from the absolute complexesList path
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
    # if all fails
    if name=='':
        name='hex'
            
    return name


def data3DList( c, key='fnac_10', inverse=0,
                rm=range(1,12), lm=range(1,12), soln=512 ):
    """
    Create an matrix: len(rec_model) * len(lig_model) * solutions
    containing the values of the info dic with given key.
    inverse - 1||0 inverse data (e.g. for energies)
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


def docking_performance( data, cutoff ):
    """
    Get a single number docking performance for all 121 combinations of rec
    and lig model.
    data - array(11x11x512) of float, fnc acts for all dockings
           each 512 values must be sorted
    cutoff - int, throw away all solutions with less than cutoff contacts
    -> array( 11 x 11 ) of float, docking 'performance' of ech rec - lig pair
    """
    ## sort data
    data = sort( data, 2 )

    ## clip values below cutoff
    data = data * logical_not( less( data, cutoff ) )

    data = data**2
    result = sum( data, 2 )

    return result


def docking_rank( performance ):
    """
    Rank the docking results from best (1) to worst (121).
    performance - array(11x11) performance measure
    ->  array(11x11)
    """
    r, l = shape( performance )
    performance = argsort( ravel( performance ) )
    
    rank = zeros( r*l )
    scale = range( 121, 0, -1 )
    
    for i in range( r*l ):
        rank[ performance[i] ] = scale[i]

    return reshape( rank, (r,l) )

    
def get_models( clst , ref_com ):
    """
    Get all rec and lig models including bound reference.
    Make sure that the residue and atom content is the same between
    the models and the bound. If not: make equal.
    clst - complex list (where first model is the free structures)
    ref_com - bound reference complex
    """
    rec_models = {}
    lig_models = {}
    for c in clst:
        if not c.info['model1'] in rec_models.keys():
            c.rec().remove( c.rec().maskH() )
            rec_models[c.info['model1']]= c.rec()
        if not c.info['model2'] in lig_models.keys():
            c.lig().remove( c.lig().maskH() )
            lig_models[c.info['model2']]= c.lig()

    rec_model_f_ref = rec_models[1]
    lig_model_f_ref = lig_models[1]

    ref_com.rec_model.remove( ref_com.rec_model.maskH() )
    rec_model_b_ref = ref_com.rec_model
    ref_com.lig_model.remove( ref_com.lig_model.maskH() )
    lig_model_b_ref = ref_com.lig_model

    ## make sure ref and models have same atoms and atom order
    ## this should always be the true, if not multidocking - skipp
    if len(rec_models)>1:
        rec_m_ind, rec_f_ind = rec_models[2].compareAtoms( rec_model_f_ref )

        if rec_m_ind != rec_f_ind:
            print 'WARNING! Free receptor reference (xray) and'
            print 'receptor models have different atom order!!'
            rec_model_f_ref = rec_model_f_ref.take( rec_f_ind )

            for k in rec_models.keys()[1:]:
                rec_models[k] = rec_models[k].take( rec_m_ind )
            lig_models[k] = lig_models[k].take( lig_m_ind )
            
    if len(lig_models)>1:
        lig_m_ind, lig_f_ind = lig_models[2].compareAtoms( lig_model_f_ref )
            
        if lig_m_ind != lig_f_ind:
            print 'Free ligand reference (xray) and '
            print 'ligand models have different atom order!!'
            lig_model_f_ref = lig_model_f_ref.take( lig_f_ind )

            for k in rec_models.keys()[1:]:
                rec_models[k] = rec_models[k].take( rec_m_ind )
            lig_models[k] = lig_models[k].take( lig_m_ind )

    ## make sure models and bound have same atoms and atom order
    rec_b2f_ind, rec_f2b_ind = rec_model_b_ref.compareAtoms( rec_model_f_ref )
    lig_b2f_ind, lig_f2b_ind = lig_model_b_ref.compareAtoms( lig_model_f_ref )

    for k in rec_models.keys():
        rec_models[k] = rec_models[k].take( rec_f2b_ind )
        lig_models[k] = lig_models[k].take( lig_f2b_ind )

    rec_model_b_ref = rec_model_b_ref.take( rec_b2f_ind )
    lig_model_b_ref = lig_model_b_ref.take( lig_b2f_ind )
            
    return rec_models, lig_models, rec_model_b_ref, lig_model_b_ref


def pw_rmsd( models_dic, ref_b, contacts ):
    """
    Pairwise rmsd calculation vs. bound and free structures.
    models_dic - model dictionary, where the free xray structure has key 1
    ref_b      - model, bound reference
    contacts   - atom contact mask
    -> dictionary, {1:[rmsd_b, rmsd_f], 2: ... }
    """
    rmsd = {}
    for k in models_dic.keys():
        rmsd_b = ref_b.rms( models_dic[k], contacts )
        rmsd_f = models_dic[1].rms( models_dic[k], contacts  )
        rmsd[k] = [ rmsd_b, rmsd_f ]
    return rmsd


def extract_info_about_best_solution( fnc_data, performance_data, rec_rms, lig_rms ):
    """
    Extract information about the docking with the highest fnac (BEST OVERALL)
    -> info dictionary
    """
    i={}
    
    ## fnc and score info
    n = fnc_data == max( ravel( fnc_data ) )
    rec_nr = nonzero( sum(sum( swapaxes( n, 0,2  ))) )[0]
    lig_nr = nonzero( sum(sum( swapaxes( n, 1,2  ))) )[0]
    best_soln = argmax( fnc_data[ rec_nr, lig_nr ] )
    best_10_fnc = sort( fnc_data[ rec_nr, lig_nr ] )[-10:]

    i['single_best_rec'] = rec_nr+1
    i['single_best_lig'] = lig_nr+1
    i['single_best_soln'] = best_soln+1
    i['single_best_fnc'] = fnc_data[ rec_nr, lig_nr, best_soln ]
    i['single_best_fnc10'] = best_10_fnc
    i['single_best_score'] = performance_data[ rec_nr, lig_nr ]

    i['single_best_rec_rms'] = rec_rms[rec_nr+1][0]
    i['single_best_lig_rms'] = lig_rms[lig_nr+1][0]
    i['all_rec_avg_rms'] = average( array( rec_rms.values() )[:,0] )
    i['all_lig_avg_rms'] = average( array( lig_rms.values() )[:,0] )
    i['all_rec_avg_sd'] = SD( array( rec_rms.values() )[:,0] )
    i['all_lig_avg_sd'] = SD( array( lig_rms.values() )[:,0] )
    i['all_avg_score'] = average( ravel( performance ) )
    i['all_sd_score'] = SD( ravel( performance ) )

    return i

    
def extract_info_about_xray( fnc_data, performance_data, rec_rms, lig_rms ):
    """
    Extract information about free xray docking
    -> info dictionary
    """
    i={}
    
    ## fnc and score info
    best_soln = argmax( fnc_data[ 0, 0 ] )
    best_10_fnc = sort( fnc_data[ 0, 0 ] )[-10:]

    i['xray_best_soln'] = int( best_soln )+1
    i['xray_best_fnc'] = fnc_data[ 0, 0, best_soln ]
    i['xray_best_fnc10'] = best_10_fnc
    i['xray_score'] = performance_data[ 0, 0 ]

    i['xray_rec_rms'] = rec_rms[1][0]
    i['xray_lig_rms'] = lig_rms[1][0]

    return i


def extract_info_about_best_models( fnc_data, performance_data, rec_rms, lig_rms ):
    """
    Extract information about the best performing docking (BEST COMBINATION)
    -> info dictionary
    """
    i={}
    
    ## fnc and score info
    n = performance_data == max( ravel( performance_data ) )
    lig_nr = nonzero( sum( n, 0 ) )[0]
    rec_nr = nonzero( sum( n, 1 ) )[0]
    best_soln = argmax( fnc_data[ rec_nr, lig_nr ] )
    best_10_fnc = sort( fnc_data[ rec_nr, lig_nr ] )[-10:]

    i['model_best_rec'] = rec_nr+1
    i['model_best_lig'] = lig_nr+1
    i['model_best_soln'] = int(best_soln)+1
    i['model_best_fnc'] = fnc_data[ rec_nr, lig_nr, best_soln ]
    i['model_best_fnc10'] = best_10_fnc
    i['model_best_score'] = performance_data[ rec_nr, lig_nr ]

    i['model_rec_rms'] = rec_rms[rec_nr+1][0]
    i['model_lig_rms'] = lig_rms[lig_nr+1][0]

    return i


def get_info_dic(  fnc_data, performance_data, rec_rms, lig_rms ):
    """
    Collect the info dictionaries into one (XRAY, BEST OVERALL, BEST COMBINATION)
    -> info dictionary
    """
    info = extract_info_about_best_models(fnc_data, performance_data, rec_rms, lig_rms)
    info.update( extract_info_about_xray( fnc_data, performance_data, rec_rms, lig_rms))
    info.update( extract_info_about_best_solution( fnc_data, performance_data, rec_rms, lig_rms))
    return info


def score_at_population_level(score, p):
    """
    the score at 'p' percent population
    """
    score = sort( ravel( score ) )
    pos = int( round( len( score ) * ( 1.0 - p )) )
    
    return sort(score)[pos]


def population_level_above_score(score, s):
    """
    the percent population above score 's'
    """
    score = ravel( score )
    return ( 1.*sum( greater_equal( score, s ) ) )/ len( score )


def get_score_dic( data ):
    r = {}
    levels = arange( 0.01, 1.01, 0.01)
    score = [ score_at_population_level( data, i ) for i in levels ]
    percent = [ population_level_above_score( data, i ) for i in levels ]
    r['levels'] = levels
    r['score_at_percent'] = score
    r['percent_above_score'] = percent
    return r


def test():
    options = _defOptions()
    dir = '/home/Bis/johan/interfaces/c24/dock_multi_0919/hex1008/'
    options['ref'] = '/home/Bis/johan/interfaces/c24/com_wet/ref.complex'
    options['cl'] = dir +'complexes_cont.cl'
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

## file names for output ...
outDir = os.path.dirname(cFile) + '/model_vs_rmsd'
if not os.path.exists( outDir ):
    os.mkdir( outDir )
outDir = outDir + '/free'
if not os.path.exists( outDir ):
    os.mkdir( outDir )
    
cName = nameFromPath( cFile )
outName = outDir + '/' + cName + '_model-rms.txt'
outName_best = outDir + '/' + cName + '_model-rms_best.txt'
outName_best_dic = outDir + '/' + cName + '_model-rms_best.dic'
outName_score = outDir + '/' + cName + '_model-rms_score.txt'
outName_score_dic = outDir + '/' + cName + '_model-rms_score.dic'

## load
flushPrint('Loading complex list \n')
cList = load( options['cl'] )
ref_com = load( options['ref'] )

## get models
rec_models, lig_models, rec_b_ref, lig_b_ref = get_models( cList, ref_com )
nr_rec = len(rec_models)
nr_lig = len(lig_models)

## guess number of solutions per single docking
single_soln = len(cList) / (nr_rec*nr_lig)

## extract data from dictionary
flushPrint('Extracting %s data  \n'%options['key'])
inverse_data = int( options['inv'])
dta = data3DList( cList, key=options['key'], inverse=inverse_data,
                  rm=range(1,nr_rec+1), lm=range(1,nr_lig+1),
                  soln=single_soln )

## contacts
ref_com_b_cast = ProteinComplex( rec_b_ref, lig_b_ref )
ref_b_resCont = ref_com_b_cast.resContacts()

rec_aContacts = rec_b_ref.res2atomMask( sum( ref_b_resCont, 1 ) )
lig_aContacts = lig_b_ref.res2atomMask( sum( ref_b_resCont, 0 ) )

## get pw_rmsds
rec_rmsd_cont = pw_rmsd( rec_models, rec_b_ref, rec_aContacts )
lig_rmsd_cont = pw_rmsd( lig_models, lig_b_ref, lig_aContacts )

rec_rmsd_bb = pw_rmsd( rec_models, rec_b_ref, rec_b_ref.maskBB() )
lig_rmsd_bb = pw_rmsd( lig_models, lig_b_ref, lig_b_ref.maskBB() )

## docking performance
performance = docking_performance( dta, 0.1 )
rank = docking_rank( performance )

## get docking info
info_dic = get_info_dic(dta, performance, rec_rmsd_cont, lig_rmsd_cont)
dump( info_dic, outName_best_dic)

out_file = open( outName_best, 'w')

head_0 = '#\t' + cName
head_1 = """
# COLUMN    DESCRIPTION
#   1     - X-ray, best fraction of native cotacts
#   2     -        solution number of best 
#   3     -        average of the 10 best fnc
#   4     -        score (performance measure)
#   5     -        receptor interface rmds to bound 
#   6     -        ligand interface rmds to bound 
#   7     - Best overall solution, best fraction of native cotacts
#   8     -                        solution number of best
#   9     -                        receptor model number
#  10     -                        ligand model number
#  11     -                        average of the 10 best fnc
#  12     -                        score (performance measure)
#  13     -                        receptor interface rmds to bound 
#  14     -                        ligand  interface rmds to bound
#  15     - Best performing models, best fraction of native cotacts
#  16     -                         solution number of best
#  17     -                         receptor model number
#  18     -                         ligand model number
#  19     -                         average of the 10 best fnc
#  20     -                         score (performance measure)
#  21     -                         receptor interface rmds to bound 
#  22     -                         ligand  interface rmds to bound
#  23     - All models, average performance score
#  24     -             standard deviation of the average performance score
#  25     -             average receptor interface rmds to bound 
#  26     -             standard deviation of receptor interface rmds
#  27     -             average ligand interface rmds to bound
#  28     -             standard deviation of ligand interface rmds"""

collected_info = [ info_dic['xray_best_fnc'],
                   info_dic['xray_best_soln'],
                   average(info_dic['xray_best_fnc10']),
                   info_dic['xray_score'],
                   info_dic['xray_rec_rms'],
                   info_dic['xray_lig_rms'],
                   
                   info_dic['single_best_fnc'],
                   info_dic['single_best_soln'],
                   info_dic['single_best_rec'],
                   info_dic['single_best_lig'],
                   average(info_dic['single_best_fnc10']),
                   info_dic['single_best_score'],
                   info_dic['single_best_rec_rms'],
                   info_dic['single_best_lig_rms'],
                   
                   info_dic['model_best_fnc'],
                   info_dic['model_best_soln'],
                   info_dic['model_best_rec'],
                   info_dic['model_best_lig'],
                   average(info_dic['model_best_fnc10']),
                   info_dic['model_best_score'],
                   info_dic['model_rec_rms'],
                   info_dic['model_lig_rms'],
                   
                   info_dic['all_avg_score'],
                   info_dic['all_sd_score'],
                   info_dic['all_rec_avg_rms'],
                   info_dic['all_rec_avg_sd'],
                   info_dic['all_lig_avg_rms'],
                   info_dic['all_lig_avg_sd'] ]


line_0 = [ '%6s'%i for i in range( len(collected_info)+1 )]
line_1 = [ '%6.2f'%i for i in collected_info ]

out_file.write( head_0 + '\n' )
out_file.write( head_1 + '\n' )
out_file.write( string.join( line_0, '\t' ) + '\n' )
out_file.write( string.join( line_1, '\t' ) + '\n' )
out_file.close()

## print data per docking
out_file = open( outName, 'w')

head_0 = '#\t' + cName
head_1 = '# %11s %24s %20s'%('rec','lig','dock')
head_2 = '#%3s %8s %10s %7s %7s %8s %5s %5s'%\
         ('nr','rms_bb','rms_cont','nr','rms_bb','rms_cont','E','rank')

out_file.write( head_0 + '\n' )
out_file.write( head_1 + '\n' )
out_file.write( head_2 + '\n' )

rec_nr, lig_nr = shape( performance )
for r in range( 1, rec_nr+1 ):
    for l in range( 1, lig_nr+1 ):
        rec_rms_b_bb = rec_rmsd_bb[r][0]
        rec_rms_f_bb = rec_rmsd_bb[r][1]
        rec_rms_b_cont = rec_rmsd_cont[r][0]
        rec_rms_f_cont = rec_rmsd_cont[r][1]
        lig_rms_b_bb = lig_rmsd_bb[l][0]
        lig_rms_f_bb = lig_rmsd_bb[l][1]
        lig_rms_b_cont = lig_rmsd_cont[l][0]
        lig_rms_f_cont = lig_rmsd_cont[l][1]
        dock_value = performance[r-1][l-1]
        dock_rank = rank[r-1][l-1]

        line = '%4i %4.2f %4.2f %4.2f %4.2f %4i %4.2f %4.2f %4.2f %4.2f %7.1f %4i'\
               %(r, rec_rms_b_bb, rec_rms_f_bb, rec_rms_b_cont, rec_rms_f_cont,
                 l, lig_rms_b_bb, lig_rms_f_bb, lig_rms_b_cont, lig_rms_f_cont,
                 dock_value, dock_rank)
        out_file.write( line + '\n')

out_file.close()


#### get score levels, write dictionary and textfile
score_dic = {}
score_dic['xray'] = get_score_dic( dta[0][0] )
score_dic['best'] = get_score_dic( dta[ info_dic['model_best_rec']-1 ]
                                   [ info_dic['model_best_lig']-1 ] )
score_dic['all'] = get_score_dic( dta )

## dump dictionary
dump( score_dic, outName_score_dic )

## write text file    
out_file = open( outName_score, 'w')

head_0 = '#\t' + cName
head_1 = '# L = level, X=xray (512), A = all (11x11x512), B = best performing model combination (512)'
head_2 = '# first line = score at level; second line = percent population above level'
line_0 = [ '%6.4f'%i for i in score_dic['xray']['levels'] ]
line_1 = [ '%6.4f'%i for i in score_dic['xray']['score_at_percent'] ]
line_2 = [ '%6.4f'%i for i in score_dic['xray']['percent_above_score'] ]
line_3 = [ '%6.4f'%i for i in score_dic['all']['score_at_percent'] ]
line_4 = [ '%6.4f'%i for i in score_dic['all']['percent_above_score'] ]
line_5 = [ '%6.4f'%i for i in score_dic['best']['score_at_percent'] ]
line_6 = [ '%6.4f'%i for i in score_dic['best']['percent_above_score'] ]
           
out_file.write( head_0 + '\n' )
out_file.write( head_1 + '\n' )
out_file.write( head_2 + '\n' )
out_file.write( 'L:  ' + string.join( line_0, '\t' ) + '\n' )
out_file.write( 'X1: ' + string.join( line_1, '\t' ) + '\n' )
out_file.write( 'X2: ' + string.join( line_2, '\t' ) + '\n' )
out_file.write( 'A1: ' + string.join( line_3, '\t' ) + '\n' )
out_file.write( 'A2: ' + string.join( line_4, '\t' ) + '\n' )
out_file.write( 'B1: ' + string.join( line_5, '\t' ) + '\n' )
out_file.write( 'B2: ' + string.join( line_6, '\t' ) + '\n' )
out_file.close()
