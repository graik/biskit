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

import Biskit.tools as T
from Biskit.Dock import ComplexList

import sys
import os.path

def syntax():
    print """
    concat_complexLists.py -i complexes1.cl complexes2.cl -o complexes_out.cl
                           -mo out_folder_for_changed_models
                           -rdic correct_rec.dic -ldic correct_lig_models.dic
    """
    sys.exit()


def pairwise_cast( models ):
    """
    atom cast all models in list with each other
    -> modified list
    """
    for i in range( len( models ) ):

        for j in range( i+1, len( models) ):

            T.flushPrint('.')

            m1 = models[i]
            m2 = models[j]

            eq_res, eq_atm = m1.equals( m2 )

            if not (eq_res and eq_atm):

                i1, i2 = m1.compareAtoms( m2 )

                delta_1 = len( m1 ) - len( i1 )
                delta_2 = len( m2 ) - len( i2 )

                models[i].keep( i1 )
                models[j].keep( i2 )

                f1 = T.stripFilename( m1.sourceFile() )
                f2 = T.stripFilename( m2.sourceFile() )

                print "Removed %i atoms from %s" % (delta_1, f1)
                print "Removed %i atoms from %s" % (delta_2, f2)

    return models


def reduceComplexList( cl ):
    """
    If cl comes from a multiple model docking, and HEX macrodocking was
    switched on, the list contains more than 512 solutions per model
    combination. Keep only the first 512 solutions.
    -> shortened or unchanged ComplexList
    """
    if len( cl ) != len(cl.models.recModels()) * len(cl.models.ligModels()) * 512:
        T.flushPrint('\nRemoving HEX solutions greater than 512.\n')
        len_old = len( cl )
        cl = cl.filter('soln', (0,512) )
        T.flushPrint('%i solutions removed.\n' % (len_old - len(cl) ) )
        T.flushPrint('%i solutions remain.\n' % len( cl ) )

    return cl

def key_for_name( reverse_model_dic, name ):
    if not name in reverse_model_dic:
        reverse_model_dic[ name ] = len( reverse_model_dic ) + 1
    return reverse_model_dic[ name ]

def reverse_dict( d ):
    r = {}
    for key, model  in d.items():
        if str(model.source) in r:
            raise Exception('cannot reverse ambiguous dictionary')
        r[ str(model.source) ] = key
    return r
    

def correct_model_numbers( cl, model_dic_1, model_dic_2 ):
    """
    Ensure that info['model1'] and info['model2'] are in accord with
    the given model dictionary.
    """
    d1 = reverse_dict( model_dic_1 )
    d2 = reverse_dict( model_dic_2 )
    counter = 0
    for c in cl:
        m1 = key_for_name( d1, str( c.rec_model.source ) )
        m2 = key_for_name( d2, str( c.lig_model.source ) )
        if m1 != c['model1']:
            c['model1'] = m1
            counter += 1
        if m2 != c['model2']:
            c['model2'] = m2
            counter += 1
    if counter:
        print "changed %i model numbers." % counter


def save_changed_models( cl, out='./changed_models' ):
    """
    To limit PVM traffic, PDBModels should always be sourced to a pickled
    object with identical coordinates/atoms/profiles. This way only the empty
    model is transfered between master and slaves.
    """
    count = 0

    for m in cl.recModels() + cl.ligModels():

##         if m.isChangedFromDisc() != (0,0):
        dir, fname = os.path.split( m.source.local() )
        m.saveAs( out + '/' + fname )
        count += 1

    if count:
        print '%i changed models have been pickled.' % count
            

def stray_models( cl ):
    """
    DEBUGGING
    Look for models that are not in rec- or lig_models dict.
    """
    stray_ligs = {}
    stray_recs = {}
    known_recs = cl.recModels()
    known_ligs = cl.ligModels()
    for c in cl:
        if c.rec_model not in known_recs:
            stray_recs[ c['model1']] = stray_recs.get( c['model1'], []) \
                       + [ c.rec_model ]
        if c.lig_model not in known_ligs:
            stray_ligs[ c['model2']] = stray_ligs.get( c['model2'], []) \
                       + [ c.lig_model ]
    return stray_recs, stray_ligs


##############
## MAIN
##############

if __name__ == '__main__':

    options = T.cmdDict( {'o':'complexes_pooled.cl',
                          'mo':'changed_models',
                          'ldic':'pcr_lig/models.dic',
                          'rdic':'pcr_rec/models.dic'} )

    if len( sys.argv ) < 2:
        syntax()

    fs = [ T.absfile( f ) for f in T.toList( options['i'] ) ]

    result = ComplexList()
    
    rec_dic = T.load( T.absfile( options['rdic'] ) )
    lig_dic = T.load( T.absfile( options['ldic'] ) )

    for f in fs:

        T.flushPrint('Loading %s ...' % f )

        cl = T.load( f )
        
        cl = reduceComplexList( cl )

        result += cl

    T.flushPrint('done\n')

    T.flushPrint('correct model numbers...')
    correct_model_numbers( result, rec_dic, lig_dic )

    T.flushPrint( '\ncasting all rec models...' )
    pairwise_cast( result.models.recModels() )

    T.flushPrint( '\ncasting all lig models...')
    pairwise_cast( result.models.ligModels() )
    T.flushPrint('done\n')

    save_changed_models( result, T.absfile( options['mo']) )

    T.flushPrint('Dumping result to %s' % T.absfile( options['o'] ) )
    T.dump( result, T.absfile( options['o'] ) )
