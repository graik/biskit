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

from Biskit.tools import *
from Biskit.PDBModel import PDBModel
import os

import Biskit.hist
from Numeric import *

import Biskit.settings
import tempfile

## tempDir = settings.tempDirShared
tempDir = Biskit.settings.tempDirLocal

from Biskit.Fold_X import *

def _use():
    print """
Syntax:  a_foldX   -c |complexes.cl| -o |out_folder| -ref |ref.complex|
        
Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def defOptions():
    return {'ref':'com_wet/ref.complex',
            'c':'dock_xray/hex01/complexes_cont.cl',
            'o':'.'}


def populationDist(score, p):
    """
    the score at 'p' percent population
    """
    pos = int( round( len(score) * p ) )
    
    return sort(score)[pos]


def binding_energy( e_rec, e_lig, e_com ):
    """
    """
    e_binding = {}
    e_tot = 0.0

    for key in e_rec.keys():
        e = e_com[key] - (  e_lig[key] +  e_rec[key] )
        e_binding[key] = e
        e_tot += e

    return e_tot, e_binding


def writeRecPdbs( complexes, start=1 ):
    """
    Write all receptor pdb files to disk
    """
    ## collect all different rec models
    models1 = {}
    for c in complexes:
        models1[c['model1']]=c['soln']
    
    ## write to disc
    flushPrint('\nWriting receptor pdbs to disc ')
    for i in models1:
        ## prepare foldX compatible model
        x = Biskit.Fold_X( )
        rec = x.foldX_model( complexes[ models1[i]-1 ].rec(), start=start )
        ## write to temp file
        pdbFile = tempfile.mktemp('.pdb')
        rec.writePdb( pdbFile )
        models1[i] = pdbFile
        flushPrint('.')
    flushPrint('\n')
    return models1


def componentEnergies( complexes ):
    """
    Calculate fold_X energies for all receptor (model1) and ligands (model2)
    """
    ## collect all different models
    rec_models = {}
    lig_models = {}
    for c in complexes:
        rec_models[c['model1']]=c['soln']
        lig_models[c['model2']]=c['soln']
        
    ## get receptor energies
    flushPrint('\nCalculating receptor energies ')
    for i in rec_models:
        r = Biskit.Fold_X( )
        e = r.energy( complexes[ rec_models[i]-1 ].rec() )
        rec_models[i] = e
        flushPrint('.')
        
    ## get ligand energies
    flushPrint('\nCalculating ligand energies ')
    for j in lig_models:
        l = Biskit.Fold_X( )
        e = l.energy( complexes[ lig_models[j]-1 ].lig() )
        lig_models[j] = e
        flushPrint('.')
            
    return rec_models, lig_models


def cleanup( fileDic ):
    for f in fileDic.values():
        try:
            os.unlink(f)
        except:
            print 'could not remove ', f

    
def test():
    options = defOptions()
    #dir = '/home/Bis/raik/data/tb/interfaces/c23/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c15/dock_xray/hex01/'
    dir = '/home/Bis/raik/data/tb/interfaces/c17/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c17/dock_0403_pc2/hex01/'
    ##dir = '/home/Bis/raik/data/tb/interfaces/c23/dock_xray/hex01'
    #dir = '/home/Bis/raik/data/tb/interfaces/c02/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c11/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c02/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c01/dock_xray/hex01/'
    #dir = '/home/Bis/raik/data/tb/interfaces/c16/dock_xray/hex01/'
    #options['l'] = dir + '../../ref_lig.pdb'
    #options['r'] = dir + '../../ref_rec.pdb'
    options['c'] = dir + '/complexes_cont.cl'
    #options['c'] = dir + '/complexes_cluster_cont.cl'
    options['ref'] = dir + '../../com_wet/ref.complex'
    options['o'] = dir

    return options


if __name__ == '__main__':
    if len(sys.argv) < 2:
        _use()
    
#    options = test()
    options = cmdDict( defOptions() )

complexes = load( absfile(options['c']) )
#complexes = complexes[:20]

## create a reference complex, get binding energy
ref_com = load( options['ref'] )
ref_lig = ref_com.lig()
ref_rec = ref_com.rec()
ref_com_model = ref_com.rec().concat( ref_com.lig())

l = Biskit.Fold_X( )
e_ref_lig = l.energy( ref_lig )

r = Biskit.Fold_X( )
e_ref_rec = r.energy( ref_rec )

c = Biskit.Fold_X( )
e_ref_com = c.energy( ref_com_model )

tot_ref, bind_ref = binding_energy( e_ref_rec, e_ref_lig, e_ref_com )


## calculate consservation score and number of native contacts for all
## complexes in complex.dic
l_energy = []
l_tot = []
l_fnc =[]

## energies for free receptor and ligand
rec_energyDic, lig_energyDic = componentEnergies( complexes) 

## receptor pdbs
rec_pdbs = writeRecPdbs( complexes, complexes[0].lig().lenResidues()+1 )

for c in complexes:

    ## energies for free ligand and receptor
    e_lig = lig_energyDic[ c.info['model2'] ]
    e_rec = rec_energyDic[ c.info['model1'] ]

    ## receptor pdb file
    rec_pdb = rec_pdbs[ c.info['model1'] ]
    
    ## prepare foldX compatible ligand model
    x = Biskit.Fold_X( )
    lig = x.foldX_model( c.lig() )
    
    ## write to temp lig
    lig_pdbFile = tempfile.mktemp('.pdb')
    lig.writePdb( lig_pdbFile, ter=0 )
    
    ## create temp complex
#    com_pdbFile = tempfile.mktemp('.pdb')
#    os.system( 'grep -v END ' + lig_pdbFile + ' >> ' + com_pdbFile )
#    os.system( 'cat ' + rec_pdb + ' >> ' + com_pdbFile )
    os.system( 'cat ' + rec_pdb + ' >> ' + lig_pdbFile )

    ## get energy
    c_com = Biskit.Fold_X( )
#    e_com = c_com.runFold_X( com_pdbFile )
    e_com = c_com.runFold_X( lig_pdbFile )

    ## tidy up
    os.unlink(lig_pdbFile)
#    os.unlink(com_pdbFile)
    
    tot, bind = binding_energy( e_rec, e_lig, e_com )

    l_tot += [tot]
    l_energy += [ bind.values() ] 

    # contact difference
    l_fnc += [ c.info['fnac_10'] ]

    print 'Nr: %i   fnac_10: %2.4f  FoldX binding energy: %2.4f'\
          %(c.info['soln'], c.info['fnac_10'], tot )

## delete temp receptor pdbs
cleanup(rec_pdbs)


## Plot results
##
import biggles

## get individual energy terms, remove non calculated terms
for t in bind.keys():
    zero_mask = [ bind[k] == 0.0 for k in bind.keys() ]
    e_val = transpose( compress( logical_not( zero_mask ), l_energy ) )
    e_key = sum( transpose( take( bind.keys(), nonzero( sum( l_energy ) ) ) ) )
    e_key = [ string.strip(key) for key in e_key ]

## initiale FramedArray
biggles.configure('fontsize_min', .5)
p = biggles.FramedArray( len(e_key), 1)

## labels etc.
p.title = '/'.join( string.split(options['o'], '/')[-4:-1] )
p.xlabel = 'fraction of native contacts'
p.ylabel = 'fold_X binding energy'

## Y-data
d0 = biggles.Points( l_fnc, l_tot, type='circle', size=1)

p_energy = []
p_label = []
p_average = []
p_reference = []
p_histogram = []
p_popLines = []

spectrum = colorSpectrum( len(e_key), 'FF0000', '00FF00' )

for i in range( len(e_key) ):
    ## plot data points
    p_energy += [ biggles.Points( l_fnc, e_val[i], color=spectrum[i],\
                                  size=1, type='filled circle') ]
    
    ## frame label
    p_label += [ biggles.PlotLabel( .50, .90, e_key[i], size=6 ) ]
    
    ## average energy line
    p_average += [ biggles.LineY(average(e_val[i]), type='shortdashed') ]
    
    ## reference complex energy line
    p_reference += [  biggles.LineY( bind_ref[str(e_key[i])] ) ]

    ## histogram inset
    inset = biggles.FramedPlot()
    inset.frame.draw_spine = 0
    inset.x2.draw_ticks = 0
    inset.x1.draw_ticks = 0
    inset.x1.draw_ticklabels = 0
    inset.y1.draw_ticklabels = 0
    inset.y1.draw_subticks = 0
    inset.y2.draw_ticks = 0
    dens = Biskit.hist.density( e_val[i], 25 )
    inset.add( biggles.Curve( dens[:,1], dens[:,0], type='solid', width=2, color=spectrum[i] ))
    p_histogram += [ inset ]

    ## population lines
    pop = [ .2, .4, .6, .8 ]
    popName = ['p2', 'p4','p6','p8' ]
    for k in range(0, len(pop)):
        popName[k] = biggles.LineY( populationDist( e_val[i], pop[k] ), type='dotted' )
    p_popLines += [ popName ]
    
## plot and save
for i in range( len(p_energy) ):
    p[i,0].add( p_energy[i] )
    p[i,0].add( p_label[i] )
    p[i,0].add( p_average[i] )
    p[i,0].add( p_reference[i] )
    p[i,0].add( biggles.Inset( (1.0,0.00), (0.8,1.0), p_histogram[i] ) )
    for j in range( shape(p_popLines)[1] ):
        p[i,0].add( p_popLines[i][j] )

p.show()
p.write_eps( options['o'] + '/foldX_plot.eps', width="18cm", height="29cm" )
