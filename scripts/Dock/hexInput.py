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
import sys
import os
from string import *

from Biskit.tools import *
from Biskit import PDBModel, PDBDope, molUtils


def centerSurfDist( model, surf_mask, mask=None ):
    """
    Calculate the longest and shortest distance from
    the center to the surface.
    mask - atoms not to be considerd
    surf_mask - atom surface mask, needed for minimum surface distance
    """
    if mask is None:
            mask = model.maskHeavy()

    ## calculate center of mass
    center, mass = model.centerOfMass(), model.mass()

    ## surface atom coordinates
    surf_xyz = compress( mask*surf_mask, model.getXyz(), 0 )
    
    ## find the atom closest and furthest away from center
    dist = sqrt( sum( (surf_xyz-center)**2 , 1 ) )
    minDist = min(dist)
    maxDist = max(dist)

    return maxDist, minDist
    

##################################################
# Usage and default parameters
##################################################

def _defaultOptions():
    return {
        'r':'receptor.pdb',
        'l':'ligand.pdb',
        'c':'',
        'rm':None,
        'lm':None,
        'sol':512 }

def _use(options):
    print """
hexInput    Create a macro file for hex. 
Syntax      hexInput -r |rec pdb| -l |lig pdb|
                     [-c |com pdb| -rm |rec model| -lm |lig model|]

              r, l   - pdb file in hex format (single or multi model)
              rm, lm - model number to use,
                       if not given perform multi model docking
              c      - a reference complex pdb file (for rmsd output)
              sol    - number of solutions to save
              
Result      Hex macro file

Default values:
    """
    for key in options.keys():
        print '\t', key, '\t', options[key]
    print


##################################################
# main function
##################################################

def main(options):

    ##################
    ## files and names
    receptorName = absfile(options['r'])
    receptorCode = stripFilename(receptorName)[0:4]
    ligandName = absfile(options['l'])          
    ligandCode = stripFilename(ligandName)[0:4]  
    complexName = absfile(options['c'])
    rm = options['rm']
    lm = options['lm']
    
    if rm or lm:
        rn = rm or '1'
        ln = lm or '1'
        baseName = receptorCode + '_' + rn + '-' + ligandCode + '_' + ln
    else:
        baseName = receptorCode + '-' + ligandCode

    ## hex macro name
    macName = baseName + '_hex.mac'

    ## hex rotation matrix output name
    outName_all = baseName + '_hex.out' 
    outName_clust = baseName + '_hex_cluster.out'


    ## load model dictionaries, if they exist in the same directory
    ## as the corresponding pdb-file and start with the same pdb-code
    ## and ends with '_model.dic'
    try:
        rec_dic = stripFilename(options['r'])[:4] + '_model.dic'
        rec_path = options['r'][:options['r'].rfind('/')]

        lig_dic = stripFilename(options['l'])[:4] + '_model.dic'
        lig_path = options['l'][:options['l'].rfind('/')]

        rec = load( absfile(rec_path + '/' + rec_dic) )
        lig = load( absfile(lig_path + '/' + lig_dic) )
        
        print 'Loading dictionaries'
        
        if type(rec) is dict:
            rec = rec[1]
            
        if type(lig) is dict:
            lig = lig[1]
        
    ## could not load dictionaty, will load pdb-file insted
    except:
        print 'Loading pdb files'
        rec = PDBModel(receptorName)
        lig = PDBModel(ligandName)

    #############################
    ## get structural information


    ## add surface profiles if not there
    if not rec.atoms.has_key('relASA'):
        flushPrint('\nCalculating receptor surface profile')
        rec_asa = PDBDope( rec )
        rec_asa.addASA()
    if not lig.atoms.has_key('relASA'):
        flushPrint('\nCalculating ligand surface profile')
        lig_asa = PDBDope( lig )
        lig_asa.addASA()

    ## surface masks, > 95% exposed
    rec_surf_mask = greater( rec.profile('relASA'), 95 )
    lig_surf_mask = greater( lig.profile('relASA'), 95 )

    ## maximun and medisn distance from centre of mass to any surface atom
    recMax, recMin = centerSurfDist( rec, rec_surf_mask )
    ligMax, ligMin = centerSurfDist( lig, lig_surf_mask )

    ## approxinate max and min center to centre distance
    maxDist = recMax + ligMax 
    minDist = recMin + ligMin
    print '\n\nReceptor and ligand max radius are %i and %i A, respectively.'\
          %( recMax, ligMax )
    print 'Receptor and ligand min radius are %i and %i A, respectively.'\
          %( recMin, ligMin )

    ## molecular separation and search range to be used in the docking
    molSep = ( maxDist + minDist ) / 2
    molRange = 2 * ( maxDist - molSep )
    print 'A molecular separation of %i A and a search range of +-%i will be used.'\
          %( molSep, molRange )

    ## determine docking mode to use
    macroDocking = 0
    if recMax > 30 and ligMax > 30:
        print '\nWARNING! Both the receptor and ligand radius is greater than 30 A.\n'     
        
    if recMax > 30:
        print '\nReceptor has a radius that exceeds 30 A -> Macro docking will be used'
        macroDocking = 1


    #####################
    ## write macro file
              
    macOpen= open(macName, 'w')

    macOpen.write('# ------------------- ' + macName + ' -----------------------\n')
    macOpen.write(' \n')
    macOpen.write('open_receptor '+receptorName+'\n')
    macOpen.write('open_ligand '+ligandName+'\n')

    if complexName[-4:] == '.pdb':
        macOpen.write('open_complex '+complexName+'\n')

    macOpen.write('\n')


    head = """
# -------------- general settings ----------------
disc_cache 1                   # disc cache on (0 off)
docking_sort_mode 1            # Sort solutions by cluster (0 by energy)
docking_cluster_mode 1         # Display all clusters (0 display best)
docking_cluster_threshold 2.00
# docking_cluster_bumps  number

# ------------ molecule orientation --------------
molecule_separation %(separation)i
commit_view """%({'separation': round(molSep)} )


    macro ="""
# -------------- macro docking -------------------
macro_min_coverage 25
macro_sphere_radius 15
macro_docking_separation 25
activate_macro_model"""
    
    
    tail = """
# -------------- docking setup -------------------
docking_search_mode 0          # full rotational search

receptor_range_angle 180       # 0, 15, 30, 45, 60, 75, 90, 180
docking_receptor_samples 720   # 362, 492, 642, 720, 980, 1280

ligand_range_angle 180
docking_ligand_samples 720

twist_range_angle 360          # 0, 15, 30, 60, 90, 180, 360
docking_alpha_samples 128      # 64, 128, 256

r12_step 0.500000              # 0.1, 0.2, 0.25, 0.5, 0.75, 1, 1.5, 2
r12_range %(range)i

docking_radial_filter 0        # Radial Envelope Filter - None
# docking_radial_filter 0        # Radial Envelope Filter - Re-entrant
# docking_radial_filter 0        # Radial Envelope Filter - Starlike

grid_size 0.600                # 0.4, 0.5, 0.6, 0.75, 1.0
#  docking_electrostatics 0       # use only surface complimentarity
docking_electrostatics 1      # use electrostatic term for scoring clusters

docking_main_scan 16     # 
docking_main_search 26

max_docking_solutions %(nr_sol)i # number of solutions to save

# -------------- post-processing ----------------
docking_refine 0    # None
#  docking_refine 1    # Backbone Bumps and volumes
#  docking_refine 2    # MM energies
#  docking_refine 3    # MM minimization

# ---------------- run docking ------------------
activate_docking
save_docking %(output_clust)s
#  save_range 1 512 ./ dock .pdb

# ------------ also save all solutions ----------
docking_sort_mode 0            # Sort solutions by energy (1 by cluster)
save_docking %(output_all)s""" \
    %({'range':round(molRange), 'output_all':outName_all,
       'nr_sol':int(options['sol']), 'output_clust':outName_clust} )
    

    macOpen.writelines( head )

    ## select certain models
    if rm or lm:
        macOpen.write('\n')
        macOpen.write('\n# -------------- select models -------------------\n')
    if rm:
        select_rec_model = 'dock_receptor_model %s'%rm
        macOpen.write( select_rec_model + '\n')
    if lm:
        select_lig_model = 'dock_ligand_model %s'%lm
        macOpen.write( select_lig_model + '\n')
        
    ## macro docking will not work with multiple models, if both are added to
    ##  the hex macro file - macrodocking will be skipped during the docking run 
    if macroDocking:
        macOpen.writelines( macro )

    macOpen.writelines( tail )
    
    macOpen.close()

    
if __name__ == '__main__':
    options = get_cmdDict(sys.argv[1:], _defaultOptions() )
    if len(sys.argv) < 2:
        _use(options)
    else:
        main(options)

