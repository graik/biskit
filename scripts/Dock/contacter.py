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
## contacter
##

import Biskit.tools as t
import sys
import Biskit.mathUtils as MU
from Biskit.PVM.hosts import cpus_all, nice_dic

from Biskit.Dock import ContactMaster, ComplexList, ComplexEvolvingList

import tempfile, os
import os.path, time, copy

def _use():
    print """
contacter: Take ComplexList, calculate contactMatrix and other stuff
           for all complexes on several nodes. Pickle ComplexList to a file.
           The result values are put into the info dict of each complex.

Syntax:    contacter [-i |complex_lst| -o |file_complex_lst|
                      -c |chunk_value| -ref |ref_complex| -v |complex_version|
                      -a -h |n_hosts| -u
                      -f |name| -s | -n |min_nice| -all -e |host1 host2..|]

Options:   -i     pickeled list of Complex objects (file name)
           -o     file name for pickled complex dictionary
           -c     chunk size (number of complexes passed to each node)
           -a     add hosts to pvm
           -h     number of nodes to use (default: all available)
           -e     exclude hosts
           -ref   pickled reference Complex for fraction of native Contacts
           -w     show xterm for each node (default: off)
           -u     only fill empty info fields, or missing keys from -f 
           -f     force calculation on sub-set of measures, current measures:
                     'fnrc_4.5', 'fnac_10', 'fnac_4.5',
                     'fnarc_9', 'fnarc_10', 'c_ratom_9', 'c_ratom_10',
                     'eProsa', 'ePairScore', 'foldX',
                     'cons_ent', 'cons_max', 'cons_abs'
                     'rms_if', 'rms_if_bb', 'xplorEnergy'
           -v     work on a previous version of each complex, only valid if
                  input is ComplexEvolvingList (e.g. status before and after
                  refinement). 0..use oldest version, -1..use latest version
           -n     renice calc to, at least, this nice value
           -s     splits complex list into sublists of this size, dumps temporary
                  contacted lists, collects result (can be resumed)
           -all   allow more than 512 solutions per model pair (keep all)


Default options:
"""
    o = defaultOptions()
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


def defaultOptions():
    return {'i':'complexes.cl', 'o':'complexes_cont.cl',
            'c':'5', 's':None, 'v':-1 }


def reduceComplexList( cl ):
    """
    If cl comes from a multiple model docking, and HEX macrodocking was
    switched on, the list contains more than 512 solutions per model
    combination. Keep only the first 512 solutions.
    -> shortened or unchanged ComplexList
    """
    if len( cl ) == len(cl.recModels()) * len(cl.ligModels()) * 512:
        return cl

    if len( cl ) < 512:
        print '\nNOTE: THE COMPLEX LIST IS SHORTER THAN 512 COMPLEXES'
        print '   COMPLEX LIST CONTAINS %i ONLY COMPLEXES\n'%len(cl)
        return cl

    t.flushPrint('\nRemoving HEX solutions greater than 512.\n')
    r = cl.filter('soln', (0,512) )
    t.flushPrint('%i solutions removed.\n' % (len( cl ) - len( r ) ) )
    t.flushPrint('%i solutions remain.\n' % len( r ) )

    return r


def checkListStatus( cl, update=0, force_keys=[], version=-1 ):
    """
    If the list is to be only updated then check that the list/sublist has
    entries that can be updated. If an unupdatable list is sent to
    the contactMaster the calculation will halt.

    update     - 1||0, are we in update mode? if not return whole list
    force_keys - 1||0, calculation is restricted to certain keys
    """
    if not update:
        return cl

    ## all complexes where force_keys are missing or where values are None
    if version == -1:
        todo = [ c for c in cl
                 if None in c.values(keys=force_keys, default=None) ]
    else:
        todo = [ c[version] for c in cl
                 if None in c.values(keys=force_keys, default=None) ]

    ## add missing force keys with value None
    if force_keys:
        for c in todo:
            missing = [ k for k in force_keys if k not in c ]

            for k in missing:
                c[ k ] = None

    if len( todo ) > 0:
        msg = "\n%i out of %i complexes are to be updated\n\n"
        t.flushPrint( msg%(len(todo), len(cl) ) )
        return cl
    else:
        t.flushPrint( "\nList contains no data that can be updated\n" )


###########################
# MAIN
###########################

if len(sys.argv) < 3:
    _use()

options = t.cmdDict( defaultOptions() )

## ## current keys used for scoring
## scoreKeys = ['eProsa', 'ePairScore', 'foldX', 'cons_ent', 'cons_max']

## load docking solutions
t.flushPrint( "\nLoading complex list %s ... " % t.absfile( options['i'] ) )
complex_lst = t.load( options['i'] )
t.flushPrint( "done\n" )

## validate and expand list of keys to be calculated
force = []
if options.has_key('f'):
    raw_force = t.toList( options['f'] )

    ## check that the key is valid
    validKeys = ['fnac_4.5', 'fnac_10', 'fnrc_4.5',
                 'fnarc_9', 'fnarc_10', 'c_ratom_9', 'c_ratom_10',
                 'eProsa', 'ePairScore', 'foldX',
                 'cons_ent', 'cons_max', 'cons_abs',
                 'rms_if', 'rms_if_bb']

    for key in raw_force:
        if key in validKeys:
            force += [ key ]
        else:
            matching = [ v for v in validKeys if key == v[:len(key)] ]
            if matching:
                force += matching
            else:    
                print '\nThe force calculate key %s is not valid.'% key
                print '### CALCULATION ABORTED ###'
                sys.exit(0)

## reduce list to standard size (if macrodocking was used)
if not 'all' in options and not isinstance( complex_lst, ComplexEvolvingList):
    complex_lst = reduceComplexList( complex_lst )

## load reference complex if given
refComplex = None
if options.has_key('ref'):
    refComplex = t.load( t.absfile( options['ref'] ) )

## remove excluded hosts
if 'e' in options:
    MU.removeFromList( cpus_all, t.toList( options['e'] ) )

host_number = int( options.get('h', len( cpus_all )) )

## add hosts to PVM ?
add_hosts = options.has_key( 'a' )

show_x = options.has_key('w')

update = options.has_key('u')

forceUpdate = options.has_key('fu')

version = int( options.get('v', -1) )

## enforce minimum nice level
if 'n' in options:
    for k in nice_dic:
        if nice_dic[ k ] < int( options['n'] ):
            nice_dic[k] = int( options['n'] )

try:
    if options['s']:
        ## determine into how many parts to split the complex list
        target_size = int(options['s'])
        n_chunks = len(complex_lst)/target_size
        if len(complex_lst) & target_size != 0:
            n_chunks += 1

        subFile_names=[]

        for n in range(1, n_chunks+1):

            subLst = complex_lst[target_size*(n-1): target_size*(n)]
            subFile = options['o'] + '_subLst_' + str(n)
            subFile_names += [ subFile ]
            if os.path.exists( subFile ):
                print '\nSub list %i (complexes %i to %i) alredy exists!'\
                      %( n, target_size*(n-1), target_size*(n) )
            else:
                print '\nWorking on sub list %i of %i (complexes %i to %i)'\
                      %( n, n_chunks+1, target_size*(n-1), target_size*(n) )

                sLst = checkListStatus(subLst, update, force, version )

                ## subList contains fields that needs to be updated
                if sLst:
                    ## initialize nodes, and start distributed calculation
                    master = ContactMaster(sLst, int( options['c'] ),
                                           cpus_all[:host_number],
                                           refComplex = refComplex,
                                           updateOnly = update,
                                           force = force,
                                           niceness = nice_dic,
                                           outFile = subFile,
                                           com_version = version,
                                           show_output = show_x,
                                           add_hosts=add_hosts)

                    t.flushPrint('Start job processing .. ')
                    master.start()

                    ## wait until master is finished
                    while not master.isFinished():
                        time.sleep( 5 )

                ## subList contains no info to be updated   
                else:
                    t.dump( subLst, subFile )

        print '\nCollecting final results...'
        complex_lst = ComplexList()
        for f in subFile_names:
            sub = t.load( f )
            complex_lst += sub
            os.unlink( f )

        t.dump( complex_lst, options['o'] )

    else:

        subLst = checkListStatus(complex_lst, update, force, version )

        if subLst:
            ## initialize nodes, and start distributed calculation
            master = ContactMaster(complex_lst, int( options['c'] ),
                                   cpus_all[:host_number],
                                   refComplex = refComplex,
                                   updateOnly = update,
                                   force = force,
                                   niceness = nice_dic,
                                   outFile = options['o'],
                                   com_version = version,
                                   show_output = show_x,
                                   add_hosts = add_hosts)
            master.start()

        else:
            t.flushPrint( "\n #### Nothing to update! #### " )

except IOError, why:
    t.errWriteln("IOError while working on %s:" % t.absfile(options['i']) \
                 + str(why) )
    t.errWriteln( t.lastErrorTrace() )
