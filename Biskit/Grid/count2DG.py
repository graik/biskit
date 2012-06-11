#!/usr/bin/python
## This script will calculate a Free Energy grid from a counts / density grid obtained with ptraj
## or other program from any Molecular dynamics simulation. It averages the counts grid, calulates expected
## densities, corrects the expected density according to the real ratio water/organic solvent, and applies a
## Boltzmann distribution to calculate the free energy.

import numpy as npy
import sys, os, time

import Biskit as bi
import Biskit.Grid.gridManager as gr

#=================================================================================================#
#                    Convert Count Grids to Free Energy                                           #
#=================================================================================================#
# USAGE>> count2DG.py pdbFile countXPLOR_file  SOLVENT  ATOM_NAME(or WAT) number_snapshots Temperature  #
#-------------------------------------------------------------------------------------------------#

_mainBoxes =   {'ANT':{'Volume':7912.290 , 'Atoms':{'N':17, 'C':17, 'WAT':209}, 'num':17, 'WAT':209},
		'WAT':{'Volume':6617.51  , 'Atoms':{'WAT':216}, 'num':216, 'WAT':216},			# TIP3PBOX Water box
                'ISO':{'Volume':8031.000 , 'Atoms':{'O':14, 'C':28, 'WAT':209}, 'num':14, 'WAT':209},
                'ETA':{'Volume':7928.500 , 'Atoms':{'O':17, 'C':17, 'WAT':209}, 'num':17, 'WAT':209},
                'ION_ISO':{'Volume':10091.000, 'Atoms':{'O':12, 'C':24, 'WAT':240}, 'num':12, 'WAT':240},
                'COO':{'Volume':10091.000, 'Atoms':{'O':24, 'WAT':240}, 'num':12, 'WAT':240},
                'COO_new':{'Volume':10893.229, 'Atoms':{'C':12, 'O':24, 'WAT':315}, 'num':12, 'WAT':315},
                'CN3':{'Volume':10091.000, 'Atoms':{'N':12, 'WAT':240}, 'num':12, 'WAT':240},
                'CN3_new':{'Volume':10893.229, 'Atoms':{'C':12,'N':12, 'WAT':315}, 'num':12, 'WAT':315},
                'ACE':{'Volume':7905.800 , 'Atoms':{'O':13, 'C':26, 'WAT':209}, 'num':13, 'WAT':209},
                'TRI':{'Volume':7988.6750, 'Atoms':{'O':13, 'C':13, 'WAT':209}, 'num':13, 'WAT':209},
                'MAM':{'Volume':7969.2880, 'Atoms':{'O':18, 'C':18, 'N':18, 'WAT':209}, 'num':18, 'WAT':209},
                'VYA':{'Volume':7889.7650, 'Atoms':{'N':19, 'C':19, 'WAT':209}, 'num':19, 'WAT':209},
                'PYR':{'Volume':7946.9200, 'Atoms':{'N':12, 'C':36, 'WAT':209}, 'num':12, 'WAT':209} # C considers meta and para C (3 of them)
               }

if len(sys.argv)<4:
    sys.exit("Error: not enough arguments provided.\nUSAGE>> count2DG.py pdbFile countXPLOR_file SOLVENT ATOM number_snapshots Temperature")

# Parse arguments
pdbFile = sys.argv[1]       # PDB file of the system
grid_file = sys.argv[2]     # XPLOR Grid File
solvent = sys.argv[3]       # Box Volume
atom = sys.argv[4]          # Number of atoms in the box
nsnap = int(sys.argv[5])    # Total number of snapshots
T = int(sys.argv[6])        # MD Temperature

# Fetch box information
print "Preparing box information..."
if _mainBoxes.has_key(solvent):
    bvol = _mainBoxes[solvent]['Volume']
    bnres = _mainBoxes[solvent]['num']
    bwat = _mainBoxes[solvent]['WAT']
    bratio = bwat / float(bnres)    # Ratio Water / solvent
    natom = _mainBoxes[solvent]['Atoms'].get(atom) or sys.exit("Invalid atom name for %s box"%solvent)
    
else:
    sys.exit("Solvent not recognized. Check boxes list in this script or add your own there.")

# Check actual ratio water/solvent in the system
print "Calculating corrections from PDB..."
pdb = bi.PDBModel(pdbFile)
pdb.fixNumbering()
pdb['resname'] = pdb.atom2resProfile('residue_name')

if solvent == 'ION_ISO':
	reslist = npy.unique(pdb['resname'])
	if 'ISO' in reslist: solvent = 'ISO'
	elif 'IPR' in reslist: solvent = 'IPR'
	else: sys.exit('ERROR with ionic box isopropanol identification.')
elif solvent == 'CN3_new':
	solvent = 'CN3'
elif solvent == 'COO_new':
	solvent = 'COO'

sys_ratio = pdb['resname'].count('WAT') / float(pdb['resname'].count(solvent))

# Calculate correction
if atom == 'WAT':
    correction = sys_ratio / bratio
else:
    correction = bratio / sys_ratio


# Work on the grid
print "Loading counts grid..."; Xgrid = gr.get_grid(grid_file)
print "Averaging points..."; Xgrid.averageData()   # Will average count grid before conversion
DGgrid = Xgrid.count2DG(bvol, natom, nsnap, T, correction=correction)  # Convert counts into Free Binding Energy
Xgrid.update(DGgrid)                            # Update instance data information

#========================= WRITE Resulting DX File ================================================
out = os.path.splitext(grid_file)[0]+"_DG.dx"
print "Writing File: ",out
Xgrid.writeDX(out)
