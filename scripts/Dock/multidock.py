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


import time

from Biskit.tools import *
from Biskit.Dock import Docker, HexParser
from Biskit.PVM.hosts import *


o = {'rdic':'pcr_rec/rec_model.dic',
     'rpdb':None,
     'ldic':'pcr_lig/lig_model.dic',
     'lpdb':None,
     'out':'hex%s',
     'com':'com/ref_hex.pdb',
     'mac':None,
     'rid':None,
     'lid':None,
     'soln':512,
     }

def syntax(o):
    if len( sys.argv ) < 2:
        print \
"""
Seperately dock several receptor models against several ligand models.

multidock.py -rdic |rec_model.dic| -ldic |lig_model.dic|
            [-rpdb |rec_hex.pdb| -lpdb |lig_hex.pdb| -com |refcomplex_hex.pdb|
             -out |outfolder| -e |excludeHost1 excludeHost2..| mac |1 or 0|
             -rid |A A ..| -lid |B| -soln |int|]

             rdic, ldic  .. dict with PCRModels indexed by 1 .. n (rec, lig)
             rpdb, lpdb  .. HEX-formatted PDB with same models (rec, lig)
             com         .. HEX-formatted PDB with reference complex
             out .. folder for results (created), may contain %s for date
             e   .. dont use these hosts
             mac .. 1|0 force the use of macro docking, if not given, the size
                    of the receptor will decide if macro docking is used.
             rid,lid .. force these chain IDs into HEX PDB file of rec / lig
             soln    .. number of solutions to keep from each docking

Default values:
    """
        for key in o.keys():
            print '\t', key, '\t', o[key]
        print

        sys.exit(0)

running = {}
def nextFreeNode():
    """-> str, next node marked free, wait until one gets available"""
    while 1:
        for h in running:
            if not running[h]:
                return h
        time.sleep( 5 )
    
def callWhenDone( runner ):
    """
    call-back that is passed on to Docker. Docker will call it whenever
    a docking run has been finished and give the RunThread as argument.
    """
    if runner.host in running:
        running[ runner.host ] = 0
    flushPrint("Finished %02i: %s %s\n" %
               (runner.owner.hexDone, runner.getName(),
                stripFilename(runner.finp) ) )

def callWhenFailed( runner ):
    """
    call-back that is passed on to Docker. Docker will call it whenever
    a docking run has failed and give the RunThread as argument.
    Todo, failed job is not re-submitted
    """
    if runner.host in running:
        errWriteln('HEX run failed on %s. Removing host.' % runner.host)
        del running[ runner.host ]

    n_rec, n_lig = runner.nRec, runner.nLig


if __name__ == '__main__':

    syntax( o )

    o = cmdDict( o )

    if o['mac'] is not None:
        o['mac'] = int( o['mac'] )

##    used_nodes = ['daredevil', 'goblin', 'joker']
##    used_nodes = ['penguin', 'spiderman', 'superman']
##    used_nodes = ['twoface', 'yiek', 'fantasio' ]
    used_nodes = [ h for h in nodes_own + nodes_shared if h in dual ]
    
##   + nodes_isd_3 + nodes_isd_2
##    nodes_fast   = len( nodes_dock ) + 1
    nodes_fast = len( used_nodes )
##    nodes_slow   = len( nodes_isd_2 )
    
    ## remove excluded hosts
    if 'e' in o:
        for h in toList( o['e'] ):
            used_nodes.remove( h )
            nodes_fast -= 1

    ## status of nodes
    for h in used_nodes:
        running[ h ] = 0

    flushPrint( "Loading dictionaries..\n" )
    rdic = load( absfile( o['rdic'] ) )
    ldic = load( absfile( o['ldic'] ) )

    flushPrint( "Preparing docking..\n" )
    d = Docker(rdic, ldic, recPdb=o['rpdb'], ligPdb=o['lpdb'], comPdb=o['com'],
               out=o['out'], macDock=o['mac'], recChainId=o['rid'],
               ligChainId=o['lid'], soln=o['soln'] )
    
    d.set_call_when_done( callWhenDone )

    flushPrint( "Start docking .. \n" )

    jobs_total = len( rdic ) * len( ldic )
    jobs_submitted = 0

    for n_rec in rdic:

        for n_lig in ldic:

##             ## don't let very slow nodes fetch late jobs
##             if jobs_submitted == (jobs_total - nodes_fast * 6):
##                 for h in nodes_isd_2:
##                     try:
##                         del running[h]
##                         flushPrint( "Removed %s from host list.\n" % h )
##                     except:
##                         pass

            ## don't let slow nodes fetch jobs that the faster ones can handle
##             if jobs_submitted == (jobs_total - nodes_fast * 4):
##                 for h in nodes_isd_3:
##                     try:
##                         del running[h]
##                         flushPrint( "Removed %s from host list.\n" % h )
##                     except:
##                         pass

            fmac, fout = d.createHexInp( n_rec, n_lig )

            host = nextFreeNode()
            running[host] = 1

            nice = nice_dic.get( host, nice_dic['default'] )

            d.runHex( fmac, ncpu=2, log=1, host=host, nice=0 )

            jobs_submitted += 1

    flushPrint( "\nLast docking submitted.\n" )

    d.waitForLastHex()
    
    flushPrint( "\nDumping result..\n")
    dump( d.result, d.out + '/complexes.cl')

    flushPrint( "\nConsistency check..\n" )

    allOK = 1
    for n_rec in rdic:
        recLst = d.result.filter('model1', n_rec)
        
        for n_lig in ldic:
            ligrecLst = recLst.filter('model2', n_lig )
            ok = len( ligrecLst ) >= 512

            if not ok:
                allOK = 0
                print "Number of solutions for %i : %i -> %i"\
                      % (n_rec, n_lig, len( ligrecLst ) )
            

    if allOK:
        flushPrint( "\nCompressing log files..\n")
        try:
            os.system( 'zip -m %s/log.zip %s/*.log' % (d.out,d.out) )
            os.system( 'zip -m %s/mac.zip %s/*.mac' % (d.out,d.out) )
            os.system( 'zip -m %s/out.zip %s/*.out' % (d.out,d.out) )
        except:
            print lastErrorDetails()
