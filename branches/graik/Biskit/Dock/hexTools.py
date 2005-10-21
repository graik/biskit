##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## last $Author$
## last $Date$
## $Revision$

import Numeric as N
import os

import Biskit.tools as t
import tempfile
from Biskit import PDBDope, molUtils


def createHexPdb_single( model, fout=None ):
    """
    Write PDB of one structure for hex.
    model - PDBModel
    fout  - str, out file name default is pdbCode + _hex.pdb
    -> str, file name of result PDB
    """
    fout = fout or model.pdbCode + '_hex.pdb'
    fout = t.absfile( fout )
    model.writePdb( fout )
    return fout
    
##     dic = { 1:model }
##     return createHexPdb( dic, fout )


def createHexPdb( modelDic, fout=None ):
    """
    write pdb for hex with models separated by MODEL%i/ENDMODEL
    modelDic - dict {int:PCRModel}
    fout - str, output name, default is pdbCode + _hex.pdb
    -> str, file name of result PDB
    """
    fout = fout or modelDic[1].pdbCode + '_hex.pdb'
    fout = t.absfile( fout )

    ## open new file
    out = open(fout, 'w')
    ## fetch name for temporary file
    tmpFile = tempfile.mktemp('_pcr2hex.pdb')  

    numbers = modelDic.keys()
    numbers.sort()
    for n in numbers:

        ## write temporary pdb and open it
        modelDic[ n ].writePdb( tmpFile )

        pdb_temp = open( tmpFile )

        out.write("MODEL%6i\n" % n)
        lines = pdb_temp.readlines()    # get all lines
        for line in lines:
            if line[:3] <> 'END':       # filter out END
                out.write(line)
        out.write("ENDMDL\n")

        ## remove temporary file
        pdb_temp.close()
        os.remove( tmpFile )

    out.write("END")
    out.close()
    return fout


def centerSurfDist( model, surf_mask, mask=None ):
    """
    Calculate the longest and shortest distance from
    the center to the surface.
    mask - atoms not to be considerd
    surf_mask - atom surface mask, needed for minimum surface distance
    """
    if mask == None:
            mask = model.maskHeavy()

    ## calculate center of mass
    center, mass = model.centerOfMass(), model.mass()

    ## surface atom coordinates
    surf_xyz = N.compress( mask*surf_mask, model.getXyz(), 0 )

    ## find the atom closest and furthest away from center
    dist = N.sqrt( N.sum( (surf_xyz-center)**2 , 1 ) )
    minDist = min(dist)
    maxDist = max(dist)

    return maxDist, minDist


def createHexInp( recPdb, recModel, ligPdb, ligModel, comPdb=None,
                  outFile=None, macDock=None, silent=0, sol=512 ):
    """
    Prepare docking of rec against lig.
    recPdb, ligPdb - str, hex-formatted PDB
    recModel, ligModel - PDBModel, get size from this one
    comPdb - str, reference PDB
    outFile - str, base of file name for mac and out
    macDoc - None -> hex decides, 1 -> force macroDock, 0-> force off
    silent - 0||1, don't print distances and macro warnings (default 0)
    
    -> (str, str, boolean), file name new HEX macro, future HEX out,
                            macro docking status
    """
    ## files and names
    recCode = t.stripFilename( recPdb )[0:4]          
    ligCode = t.stripFilename( ligPdb )[0:4]

    outFile = outFile or recCode + '-' + ligCode

    ## hex macro name
    macName = t.absfile( outFile + '_hex.mac' )

    ## hex rotation matrix output name
    outName_all = t.absfile( outFile + '_hex.out'  )
    outName_clust = t.absfile( outFile + '_hex_cluster.out')

    ## add surface profiles if not there
    if not recModel.aProfiles.has_key('relASA'):
        #t.flushPrint('\nCalculating receptor surface profile')
        rec_asa = PDBDope( recModel )
        rec_asa.addASA()
    if not ligModel.aProfiles.has_key('relASA'):
        #t.flushPrint('\nCalculating ligand surface profile')
        lig_asa = PDBDope( ligModel )
        lig_asa.addASA()

    ## surface masks, > 95% exposed
    rec_surf_mask = N.greater( recModel.profile('relASA'), 95 )
    lig_surf_mask = N.greater( ligModel.profile('relASA'), 95 )

    ## maximun and medisn distance from centre of mass to any surface atom
    recMax, recMin = centerSurfDist( recModel, rec_surf_mask )
    ligMax, ligMin = centerSurfDist( ligModel, lig_surf_mask )

    ## approxinate max and min center to centre distance
    maxDist = recMax + ligMax 
    minDist = recMin + ligMin

    ## molecular separation and search range to be used in the docking
    molSep = ( maxDist + minDist ) / 2
    molRange = 2 * ( maxDist - molSep )

    if not silent:
        print 'Docking setup: %s\nRecMax: %.1f RecMin: %.1f\nLigMax: %.1f LigMin: %.1f\nMaxDist: %.1f MinDist: %.1f\nmolecular_separation: %.1f r12_range: %.1f\n'%(outFile, recMax, recMin, ligMax, ligMin, maxDist, minDist, molSep, molRange)
    
    if recMax > 30 and ligMax > 30 and not silent:
        print '\nWARNING! Both the receptor and ligand radius is ',
        print 'greater than 30A.\n'     

    ## determine docking mode to use
    macroDocking = 0

    if macDock==None:
        if recMax > 35 and not silent:
            print '\nReceptor has a radius that exceeds 35A ',
            print '-> Macro docking will be used'
            macroDocking = 1
    else:
        macroDocking = macDock

    #####################
    ## write macro file

    macOpen= open( macName, 'w')

    macOpen.write('# -- ' + macName + ' --\n')
    macOpen.write(' \n')
    macOpen.write('open_receptor '+ t.absfile(recPdb) +'\n')
    macOpen.write('open_ligand '+ t.absfile(ligPdb) +'\n')

    if comPdb and comPdb[-4:] == '.pdb':
        macOpen.write('open_complex '+comPdb+'\n')

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

receptor_range_angle  180      # 0, 15, 30, 45, 60, 75, 90, 180
docking_receptor_samples 720   # 362, 492, 642, 720, 980, 1280

ligand_range_angle  180
docking_ligand_samples 720

twist_range_angle 360          # 0, 15, 30, 60, 90, 180, 360
docking_alpha_samples 128      # 64, 128, 256

r12_step 0.500000              # 0.1, 0.2, 0.25, 0.5, 0.75, 1, 1.5, 2
r12_range %(range)i

docking_radial_filter 0        # Radial Envelope Filter - None

grid_size 0.600                # 0.4, 0.5, 0.6, 0.75, 1.0
# docking_electrostatics 0       # use only surface complimentarity
docking_electrostatics 1      # use electrostatic term for scoring clusters

docking_main_scan 16     # 
docking_main_search 26

max_docking_solutions %(nr_sol)i # number of solutions to save

# -------------- post-processing ----------------
docking_refine 0    # None
#  docking_refine 1    # Backbone Bumps
#  docking_refine 2    # MM energies
#  docking_refine 3    # MM minimization

# ---------------- run docking ------------------
activate_docking
#  save_docking %(output_clust)s
#  save_range 1 512 ./ dock .pdb

# ------------ also save all solutions ----------
docking_sort_mode 0            # Sort solutions by energy (1 by cluster)
save_docking %(output_all)s""" \
    %({'range':round(molRange), 'output_all':outName_all,
       'nr_sol':int(sol), 'output_clust':outName_clust} )

    macOpen.writelines( head )

    ## macro docking will not work with multiple models, if both are added to
    ## the hex macro file - macrodocking will be skipped during the docking run
    if macroDocking:
        macOpen.writelines( macro )

    macOpen.writelines( tail )

    macOpen.close()

    return macName, outName_all, macroDocking
